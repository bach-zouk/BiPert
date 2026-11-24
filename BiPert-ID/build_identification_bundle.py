"""build_identification_bundle.py  plan | train <variant> | query <dataset> | tables

Builds the data bundle of the BiPert-ID notebook from the rebuilt ExPert benchmark datasets
(one folder per dataset: pseudobulk sums over all valid perturbations, the benchmark cells, ExPert's
held-out test cells per fold). Authors' record: a reader of the notebook does not need this script.

The bundle holds, for ExPert's full benchmark (22 datasets, 12 folds, 34 held-out views):

  benchmark.json                 folds, held-out views, the training contexts of each fold
  candidates.txt                 the 1,157 candidate genes
  expert_input_panel.txt         the 3,096 genes ExPert reads
  train_<variant>.npz            one per training context (cell line; T cells by condition). A context whose
                                 members change between folds (a T-cell donor held out) has one file per
                                 member set. Holds summary statistics only:
       genes        ExPert vocabulary column of each gene
       labels       candidate genes measured in this context (>= 10 cells)
       effect       labels x genes, log-normalised mean expression minus control (float32)
       mad          per gene, robust spread of the effect over ALL knockdowns of the context
       mean_effect  per gene, mean effect over ALL knockdowns of the context
       ctrl_profile per gene, log-normalised mean expression of the control cells
       n_knockdowns number of knockdowns behind mad / mean_effect
  query_<dataset>.npz            one per held-out dataset:
       genes, symbols            stored genes (vocabulary column, symbol)
       counts_*                  CSR raw counts of ExPert's test cells on the stored genes
       total                     per cell, total counts over all genes of the dataset
       label                     per cell, the knocked-down gene
       fold<k>                   rows of the test cells of fold k
       ctrl_mu, ctrl_sd          per stored gene, mean / SD of log-normalised expression over control cells
       ctrl_profile              per stored gene, log-normalised expression of the pooled control cells
       sd_fill                   SD used as placeholder for candidates whose gene is not measured
       n_ctrl
  expert_pooled.csv.gz           ExPert's released pseudobulk rank of the true gene, per fold / dataset / knockdown
  reference_ranks.csv.gz         the authors' cluster run: pseudobulk rank of the true gene, BiPert-ID and its
                                 unweighted version (the notebook checks that it reproduces them)

Stored genes: union over folds of the genes shared by all training datasets of the fold (U), restricted to the
genes of the file; query files additionally hold every candidate's own gene.
"""
import json, os, sys, time
import numpy as np, pandas as pd, scipy.sparse as sps
sys.path.insert(0, os.environ.get('EXPERT_EVAL_DIR', 'eval'))  # folder holding the benchmark helper fbe_common.py
import fbe_common as F

OUT = os.environ.get('BIPERT_OUT', 'identification_data')
os.makedirs(OUT, exist_ok=True)
T0 = time.time()


def log(*a):
    print(time.strftime('%H:%M:%S'), f'[{time.time() - T0:.0f}s]', *a, flush=True)


def vocab(d):
    v = F.load_genes(d).vocab_col.to_numpy(); return set(v[v >= 0].tolist())


def variant_name(c, members):
    """T_Rest with all four donors -> 'T_Rest'; without donor D2 -> 'T_Rest_noD2'."""
    if not c.startswith('T_'):
        return c
    donors = sorted(m.split('_')[1] for m in members)
    missing = [x for x in ('D1', 'D2', 'D3', 'D4') if x not in donors]
    return c if not missing else c + '_no' + ''.join(missing)


def plan():
    """Folds, views, contexts, variants and the gene union U. Pure bookkeeping (gene lists only)."""
    st = F.dataset_status(); FO = F.folds(); V = {d: vocab(d) for d in F.all_datasets()}
    U = set(); variants = {}; folds = {}
    for k in sorted(FO):
        complete, held, train, missing = F.fold_status(k, st); assert complete, (k, missing)
        gt = None
        for d in train:
            gt = V[d] if gt is None else gt & V[d]
        U |= gt
        members = {}
        for d in train:
            members.setdefault(F.context_of(d), []).append(d)
        ctx = {}
        for c in sorted(members):
            vn = variant_name(c, members[c]); ctx[c] = vn
            assert variants.setdefault(vn, sorted(members[c])) == sorted(members[c])
        kind = {'cell_line': 'new cell line', 'condition': 'same donor, new condition'}
        views = [dict(dataset=d, context=F.context_of(d),
                      kind=kind.get(FO[k]['axis'], 'new donor' if d.startswith('Marson2025_') else 'new cell line')) for d in held]
        folds[k] = dict(held_out=FO[k]['held_out'], axis=FO[k]['axis'], views=views, training_contexts=ctx,
                        training_datasets=train, n_genes_shared_by_training=len(gt))
    return dict(folds=folds, variants=variants, U=np.array(sorted(U)))


def build_train(vn):
    P = plan(); members = P['variants'][vn]; U = P['U']
    g = None
    for d in members:
        g = vocab(d) if g is None else g & vocab(d)
    genes = np.array(sorted(set(U.tolist()) & g))
    z = F.load_training_context(members, genes, log=log)
    Lp = np.log1p(1e4 * z['M'] / np.maximum(z['tot'], 1e-12)[:, None])
    Lc = np.log1p(1e4 * z['cM'] / np.maximum(z['ctot'], 1e-12))
    keep = np.asarray(z['n_cells']) >= 10; perts = np.asarray(z['perts']).astype(str)[keep]; Ec = Lp[keep] - Lc
    del Lp
    mad = 1.4826 * np.median(np.abs(Ec - np.median(Ec, 0)), 0)
    mean_effect = Ec.mean(0)
    classes = F.bench_classes()
    first = ~pd.Index(perts).duplicated(); idx = pd.Index(classes).get_indexer(perts); m = first & (idx >= 0)
    o = np.argsort(idx[m]); labels = classes[idx[m]][o]; E = Ec[m][o].astype(np.float32)
    np.savez_compressed(f'{OUT}/train_{vn}.npz', genes=genes.astype(np.int32), labels=labels.astype(str), effect=E,
                        mad=mad, mean_effect=mean_effect, ctrl_profile=Lc, n_knockdowns=int(keep.sum()),
                        members=np.array(members))
    log(f'train_{vn}: members {members}; {len(genes)} genes; {int(keep.sum())} knockdowns with >= 10 cells; '
        f'{len(labels)} candidates; {os.path.getsize(f"{OUT}/train_{vn}.npz") / 1e6:.1f} MB')


def build_query(d):
    P = plan(); U = P['U']; classes = F.bench_classes()
    ks = [k for k, f in P['folds'].items() if d in [v['dataset'] for v in f['views']]]
    g = F.load_genes(d); vc = g.vocab_col.to_numpy(); sym = g.symbol_vocab.astype(str).to_numpy()
    mapped = np.flatnonzero(vc >= 0)
    fi = {sym[j]: j for j in mapped}; tfull = np.array([fi.get(p, -1) for p in classes])   # raw column of each candidate's gene
    inU = mapped[np.isin(vc[mapped], U)]
    raw = np.unique(np.concatenate([inU, tfull[tfull >= 0]]))
    raw = raw[np.argsort(vc[raw], kind='stable')]                    # stored genes, sorted by vocabulary column
    obs = pd.read_csv(f'{F.FB}/out/{d}/bench_obs.csv.gz', usecols=['bench_row', 'label_cleaned', 'is_ctrl'])
    oq = np.load(f'{F.FB}/out/{d}/ood_queries.npz', allow_pickle=True)
    q = {k: np.sort(oq[f'fold{k}']) for k in ks}; uq = np.unique(np.concatenate(list(q.values())))
    ci = np.flatnonzero(obs.is_ctrl.to_numpy())
    X = sps.load_npz(f'{F.FB}/out/{d}/bench_cells.npz').tocsr(); assert X.shape == (len(obs), len(g))
    Xq = X[uq]; Xc = X[ci]; del X
    tq = np.asarray(Xq.sum(1)).ravel().astype(np.float64); tc = np.asarray(Xc.sum(1)).ravel().astype(np.float64)
    Lc = np.log1p(1e4 * Xc[:, raw].toarray().astype(np.float64) / tc[:, None])
    mu = Lc.mean(0); sd = Lc.std(0); del Lc
    prof = np.log1p(1e4 * np.asarray(Xc[:, raw].sum(0), dtype=np.float64).ravel() / max(float(tc.sum()), 1e-12))
    sd_fill = float(np.log1p(1e4 * Xc[:, [0]].toarray().astype(np.float64).ravel() / tc).std())
    C = Xq[:, raw].tocsr(); C.sort_indices()
    mx = C.data.max(); assert (C.data == np.round(C.data)).all() and mx < 2 ** 32
    dt = np.uint16 if mx < 2 ** 16 else np.uint32
    save = dict(genes=vc[raw].astype(np.int32), symbols=sym[raw].astype(str), counts_data=C.data.astype(dt),
                counts_indices=C.indices.astype(np.uint16 if len(raw) < 2 ** 16 else np.int32),
                counts_indptr=C.indptr.astype(np.int64), counts_shape=np.array(C.shape), total=tq,
                label=obs.label_cleaned.astype(str).to_numpy()[uq], ctrl_mu=mu, ctrl_sd=sd, ctrl_profile=prof,
                sd_fill=sd_fill, n_ctrl=int(len(ci)))
    for k in ks:
        save[f'fold{k}'] = np.searchsorted(uq, q[k]).astype(np.int32)
    np.savez_compressed(f'{OUT}/query_{d}.npz', **save)
    log(f'query_{d}: folds {ks}; {C.shape[0]} test cells x {C.shape[1]} genes, nnz {C.nnz}, max count {int(mx)}; '
        f'{len(ci)} control cells; candidates with a measured gene {(tfull >= 0).sum()}; '
        f'{os.path.getsize(f"{OUT}/query_{d}.npz") / 1e6:.1f} MB')


def tables():
    P = plan()
    json.dump(dict(folds={str(k): f for k, f in P['folds'].items()},
                   variants=P['variants'], n_genes_union=int(len(P['U']))), open(f'{OUT}/benchmark.json', 'w'), indent=1)
    open(f'{OUT}/candidates.txt', 'w').write('\n'.join(F.bench_classes()) + '\n')
    _, EG = F.expert_released(); EG = EG[EG['thr'].isna() & (EG['agg'] == 'mean')]
    rows, ref = [], []
    for k, f in P['folds'].items():
        for v in f['views']:
            d = v['dataset']
            e = EG[(EG.fold == k) & (EG.dataset == d)]; assert len(e) and e.label.is_unique
            rows.append(pd.DataFrame(dict(fold=k, dataset=d, label=e.label.to_numpy(), rank=e['rank'].to_numpy(float))))
            w = pd.read_parquet(f'{F.OUT}/fold{k}/W/{d}.groups.parquet'); w = w[w['agg'] == 'mean'].set_index('label')['rank'].astype(float)
            o = pd.read_parquet(f'{F.OUT}/fold{k}/{d}.groups.parquet'); o = o[o['agg'] == 'mean'].set_index('label')['rank'].astype(float)
            assert set(w.index) == set(e.label) == set(o.index)
            ref.append(pd.DataFrame(dict(fold=k, dataset=d, label=w.index, rank_bipert_id=w.to_numpy(),
                                         rank_unweighted=o.reindex(w.index).to_numpy())))
    pd.concat(rows).to_csv(f'{OUT}/expert_pooled.csv.gz', index=False)
    pd.concat(ref).to_csv(f'{OUT}/reference_ranks.csv.gz', index=False)
    panel = np.load(f'{F.P}/scores/expert_HepG2.npz', allow_pickle=True)['panel'].astype(str)
    for line in ('Jurkat', 'K562_essential', 'rpe1'):
        assert set(np.load(f'{F.P}/scores/expert_{line}.npz', allow_pickle=True)['panel'].astype(str)) == set(panel)
    open(f'{OUT}/expert_input_panel.txt', 'w').write('\n'.join(sorted(panel)) + '\n')
    log(f'tables: {len(P["folds"])} folds, {sum(len(f["views"]) for f in P["folds"].values())} views, '
        f'{len(P["variants"])} training files, union genes {len(P["U"])}, panel {len(panel)}')


if __name__ == '__main__':
    a = sys.argv[1]
    if a == 'plan':
        P = plan()
        for vn, m in P['variants'].items():
            print('train', vn, m)
        for d in sorted({v['dataset'] for f in P['folds'].values() for v in f['views']}):
            print('query', d)
        print('U', len(P['U']))
    elif a == 'train':
        build_train(sys.argv[2])
    elif a == 'query':
        build_query(sys.argv[2])
    elif a == 'tables':
        tables()
