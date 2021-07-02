import numpy as np
import os
import random

from collections import Counter
from functools import partial
from multiprocessing import Process, Value, Array, Queue
from queue import Empty
from scipy.spatial import distance
from sklearn_extra.cluster import KMedoids
from scipy.spatial.distance import cdist
from time import sleep


def run_kmedoids(data, k, sampsize, random_seed, distance_metric='euclidean', mean_inertia=False,
                 km_maxiter=300, verbose=False, **kwargs):
    """
    Identify `k` medoids from a `sampsize`-sized samples of points from `data` using the provided
    `random_seed` to seed the RNG and the provided `distance_metric`. The Global inertia is 
    calculated as either the sum or mean (mean_inertia=True/false) of each point in the global array 
    to its closest medoid in the identified set.

    :param np.ndarray data: 2D Array of points
    :param int k: Number of medoids to calculate
    :param int sampsize: Sample size
    :param int random_seed: random seed for RNG
    :param str distance_metric: Distance metric for the KMedoid algorithm
    :param bool mean_inertia: Use mean for inertia instead of sum?
    :param int km_maxiter: Max iterations for the KMedoid algorithm
    :param bool verbose: Logging verbosity
    :param **kwargs: Futureproofing
    """
    np.random.seed(seed=random_seed)
    _data_idxs = np.random.choice(range(data.shape[0]), size=sampsize, replace=False)
    _data = data[_data_idxs,:]
    km = KMedoids(n_clusters=k,
                  metric=distance_metric,
                  max_iter=km_maxiter,
                  random_state=random_seed).fit(_data)
    medoid_indices = _data_idxs[km.medoid_indices_]
    # It's cheaper to just run the additional `k` distance calculations than it is to filter data to 
    # remove the `k` medoids before calculating pairwise distances.
    # Proof: get_inertia_* is basically the same as the below lines, with and without filtering data
    #        for the k medoids.
    #     >>> data.shape
    #     (20000, 30)
    #     >>> medoid_indices.shape
    #     (10,)
    #     >>> timeit.timeit("get_inertia_filter(data, medoid_indices, distance_metric)",
    #                                           number=20, globals=globals())
    #     1.315036368000051
    #     >>> timeit.timeit("get_inertia_nofilter(data, medoid_indices, distance_metric)",
    #                                             number=20, globals=globals())
    #     0.1129797050000434
    global_inertia = cdist(data,
                           data[medoid_indices,:],
                           metric=distance_metric,
                           **kwargs).min(axis=1).sum()
    if mean_inertia:
        # The denominator needs to be subtracted for the additional distances calculated
        global_inertia /= (data.shape[0] - k)
    return medoid_indices, global_inertia

def run_kmedoids_multiprocessing(tid, medoid_indices, inertia, data, k, sampsize, random_seed,
                                 distance_metric='euclidean', mean_inertia=False, km_maxiter=300,
                                 verbose=False):
    """
    """
    n = 0
    prefix = f'WORKER_{tid}:'
    print(f'{prefix} Spawning worker... Waiting 1 second for queue to populate.')
    sleep(1)
    while True:
        try:
            _random_seed = random_seed.get(block=True, timeout=5)
        except Empty:
            break
        _medoid_indices, _inertia = run_kmedoids(data=data,
                                                 k=k,
                                                 sampsize=sampsize,
                                                 random_seed=_random_seed,
                                                 distance_metric=distance_metric,
                                                 mean_inertia=mean_inertia,
                                                 km_maxiter=km_maxiter,
                                                 verbose=verbose)
        with inertia.get_lock():
            if _inertia < inertia.value:
                print(f'{prefix} identified medoids with lower inertia. '
                      f'Current inertia_min is {_inertia}')
                for i, v in enumerate(_medoid_indices):
                    medoid_indices[i] = v
                inertia.value = _inertia
                if verbose:
                    print(f'{prefix} Medoids found => [', *medoid_indices, ']', sep=' ')
        n += 1
    print(f'{prefix} Processed {n} samples in total. No more samples to process. Signing off....')

def _standardize(x):
    """
    On a per-column-basis, subtract mean and divide by mean absolute deviation
    This is exactly how R clara does it :shrug:
    :param np.ndarray x:
    """
    assert isinstance(x, np.ndarray) and len(x.shape) == 2
    return np.apply_along_axis(lambda y: (y - np.nanmean(y))/np.nanmean(np.abs(y - np.nanmean(y))),
                               axis=1, arr=x)

def clara(data,
          k,
          distance_metric='euclidean',
          samples=5,
          km_maxiter=300,
          sampsize=None,
          random_seed=21212,
          standardize=True,
          mean_inertia=False,
          nworkers=1,
          verbose=True,
          **kwargs):
        assert isinstance(data, np.ndarray)
        assert len(data.shape) == 2 and data.shape[1] < data.shape[0]
        assert k < data.shape[0]
        if sampsize is None:
            sampsize = min(data.shape[0], 40 + 2 * k)
        else:
            assert isinstance(sampsize, int) and sampsize < data.shape[0]
        assert isinstance(km_maxiter, int)
        assert distance_metric in ("euclidean", "manhattan", "jaccard")
        assert isinstance(standardize, bool)
        assert isinstance(mean_inertia, bool)
        assert isinstance(samples, int) and samples > 0
        assert isinstance(random_seed, int)
        assert isinstance(nworkers, int) and nworkers > 0 and samples > 0
        if nworkers > (os.cpu_count()*2):
            print('WARNING: Specified more workers to use than ncpus * 2')
        if verbose:
            print('Running with parameters: ',
                  f'data            : {type(data)} of shape {data.shape}',
                  f'k               : {k}',
                  f'distance_metric : {distance_metric}',
                  f'samples         : {samples}',
                  f'km_maxiter      : {km_maxiter}',
                  f'sampsize        : {sampsize}',
                  f'random_seed     : {random_seed}',
                  f'standardize     : {standardize}',
                  f'mean_inertia    : {mean_inertia}',
                  f'nworkers        : {nworkers}',
                  f'verbose         : {verbose}',
                  sep='\n\t')
        # Seed the random number generator here so it's possible to get similar outputs on
        # subsequent function calls
        np.random.seed(seed=random_seed)
        random_seeds = [np.random.randint(2**32) for _ in range(samples)]

        if nworkers > 1:
            print(f'Running in multiprocessing mode with {nworkers} workers')
            medoid_indices = Array('i', range(k))
            inertia = Value('d', np.inf)
            random_queue = Queue()
            for rs in random_seeds:
                random_queue.put(rs)

            workers = [Process(target=run_kmedoids_multiprocessing,
                               args = (nw+1, # tid
                                       medoid_indices, # medoid_indices
                                       inertia,  # inertia
                                       data,
                                       k,
                                       sampsize,
                                       random_queue,  # , random_seed,
                                       distance_metric,
                                       mean_inertia,
                                       km_maxiter,
                                       verbose))
                            for nw in range(nworkers)]
            for p in workers:
                p.start()
            for p in workers:
                p.join()
            medoid_indices = list(medoid_indices)
        else:
            print(f'Running in single processing mode')
            inertia = np.inf
            medoid_indices = []
            for rs in random_seeds:
                _medoid_indices, _inertia = run_kmedoids(data=data,
                                                         k=k,
                                                         sampsize=sampsize,
                                                         random_seed=rs,
                                                         distance_metric=distance_metric,
                                                         mean_inertia=mean_inertia,
                                                         km_maxiter=km_maxiter,
                                                         verbose=verbose)
                if _inertia < inertia:
                    print('Identified medoids with lower inertia. '
                          f'Current inertia_min is {_inertia}')
                    inertia = _inertia
                    medoid_indices = _medoid_indices
                    if verbose:
                        print(f'Medoids found => [', *medoid_indices, ']', sep=' ')
        
        medoid_indices = list(medoid_indices)
        dm = cdist(data,
                   data[medoid_indices,:],
                   metric=distance_metric,
                   **kwargs)
        clusters = dm.argmin(axis=1)
        if verbose:
            cluster_counts = Counter(clusters)
            rjust = len(str(max(cluster_counts.values())))
            print_string = ('{{k:<3}}: {{counts:>{rjust}}} ({{frac:0.4f}})').format(rjust=rjust)
            print('Cluster counts were as follows:')
            for k in sorted(cluster_counts, key=lambda x: int(x)):
                print(print_string.format(k=k,
                                          counts=cluster_counts[k],
                                          frac=cluster_counts[k]/data.shape[0]))
        return clusters, medoid_indices