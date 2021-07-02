def subsample_by_var(adata, var, subsample=0.1):
    assert var in adata.obs
    idxs_to_keep = \
        sorted(
            np.hstack(
                adata.obs.groupby(var).apply(
                    lambda x: np.random.choice(x.index, 
                                               int(subsample * x.index.shape[0]), 
                                               replace=False)).values), key=lambda x: int(x))
    return adata.loc[idxs_to_keep,:].copy()
