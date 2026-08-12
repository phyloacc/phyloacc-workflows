#!/usr/bin/env python3
"""
lib/clustering.py - conserved-region clustering methods for phyloP significant sites.

Each method takes the positions of FDR-significant conserved sites (0-based) and returns
conserved regions as (start, end, n_sites) half-open intervals. Experimental /
comparison-stage - see analyses/phylop-clustering-comparison/. NOT wired into the pipeline
or config templates yet (only useful on deep trees, where phyloP produces conserved sites).

Methods:
- windowed:        fixed-size bins; bins with >= min_sites are conserved; merge adjacent.
- hmm:             2-state online HMM over the per-position conserved/not symbol stream
                   (adapted from workflow/legacy/scripts/conserved_elements_hmm.py).
- hdbscan_cluster: 1-D density clustering of site positions (imports hdbscan lazily).

The existing gap-merge method lives in lib.intervals.cluster_sites (used as a baseline in
the comparison). numpy only, except hdbscan_cluster.
"""

import numpy as np


def as_positions(sites):
    """Accept an iterable of ints or (start, end, ...) tuples; return sorted (starts, ends)
    int arrays (ends default to start+1 for scalar inputs)."""
    starts, ends = [], []
    for s in sites:
        if isinstance(s, (int, np.integer)):
            starts.append(int(s)); ends.append(int(s) + 1)
        else:
            starts.append(int(s[0])); ends.append(int(s[1]))
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    order = np.argsort(starts, kind="mergesort")
    return starts[order], ends[order]


# ---------------------------------------------------------------------------
# 1. windowed - fixed-size bins
# ---------------------------------------------------------------------------

def windowed(sites, length, window_bp=50, min_sites_per_window=5):
    """Tile [0, length) into fixed windows of window_bp; a window is 'conserved' if it
    holds >= min_sites_per_window significant sites; merge runs of adjacent conserved
    windows into regions. Returns [(start, end, n_sites), ...]."""
    starts, _ = as_positions(sites)
    n_windows = (length + window_bp - 1) // window_bp
    counts = np.bincount((starts // window_bp).astype(np.int64), minlength=n_windows)[:n_windows]
    conserved = counts >= min_sites_per_window

    regions = []
    i = 0
    while i < n_windows:
        if not conserved[i]:
            i += 1
            continue
        j = i
        while j + 1 < n_windows and conserved[j + 1]:
            j += 1
        start = i * window_bp
        end = min((j + 1) * window_bp, length)
        regions.append((start, end, int(counts[i:j + 1].sum())))
        i = j + 1
    return regions


# ---------------------------------------------------------------------------
# 2. hmm - 2-state online HMM (adapted from legacy conserved_elements_hmm.py)
# ---------------------------------------------------------------------------

def hmm(sites, length, t0_0=0.9, t1_1=0.8, e0_0=0.8, e1_1=0.5, s0=0.9,
        min_len=20, max_len=10000):
    """2-state online (forward-filtering) HMM. State 0 = outside a conserved element,
    1 = inside; symbol 0 = non-conserved position, 1 = conserved site. Walks every
    position in [0, length), assigns the max-marginal state per position, and returns
    runs of state==1 whose length is strictly between min_len and max_len.

    Probabilities: t0_0 = P(stay outside), t1_1 = P(stay inside), e0_0 = P(emit
    non-conserved | outside), e1_1 = P(emit conserved | inside), s0 = P(start outside).
    Forward filtering is done in the probability domain with per-step renormalization
    (numerically stable), matching the legacy script's `predict`.
    """
    starts, ends = as_positions(sites)
    symbol = np.zeros(length, dtype=np.int8)
    for st, en in zip(starts, ends):
        if st < length:
            symbol[st:min(en, length)] = 1

    def clamp(v):
        return min(max(v, 1e-12), 1.0 - 1e-12)

    t0_0, t1_1, e0_0, e1_1, s0 = map(clamp, (t0_0, t1_1, e0_0, e1_1, s0))
    t01, t10 = 1.0 - t0_0, 1.0 - t1_1
    # emit[state][symbol]
    e0 = (e0_0, 1.0 - e0_0)
    e1 = (1.0 - e1_1, e1_1)

    regions = []
    p0 = p1 = None
    in_elem = False
    elem_start = 0
    sym_arr = symbol  # local ref for speed

    for i in range(length):
        sym = sym_arr[i]
        if p0 is None:
            n0 = s0 * e0[sym]
            n1 = (1.0 - s0) * e1[sym]
        else:
            n0 = (p0 * t0_0 + p1 * t10) * e0[sym]
            n1 = (p0 * t01 + p1 * t1_1) * e1[sym]
        z = n0 + n1
        p0 = n0 / z
        p1 = n1 / z
        state1 = p1 > p0

        if state1 and not in_elem:
            in_elem = True
            elem_start = i
        elif not state1 and in_elem:
            in_elem = False
            elem_len = i - elem_start
            if min_len < elem_len < max_len:
                regions.append((elem_start, i, int(symbol[elem_start:i].sum())))

    if in_elem:
        elem_len = length - elem_start
        if min_len < elem_len < max_len:
            regions.append((elem_start, length, int(symbol[elem_start:length].sum())))
    return regions


# ---------------------------------------------------------------------------
# 3. hdbscan_cluster - 1-D density clustering
# ---------------------------------------------------------------------------

def hdbscan_cluster(sites, min_cluster_size=10, min_samples=None, cluster_selection_method="eom",
                    cluster_selection_epsilon=0.0, min_len=1, min_sites=1):
    """Density-cluster site positions in 1-D with HDBSCAN; each cluster becomes a region
    spanning its members (noise, label -1, is dropped). Returns [(start, end, n_sites), ...]
    sorted by start. Imports hdbscan lazily so the rest of the module works without it.

    cluster_selection_method: 'eom' (default, most-stable) or 'leaf' (finer, smaller clusters
    - less gap-bridging). cluster_selection_epsilon: distance below which clusters are merged
    (0 = off; larger = more merging)."""
    import hdbscan

    starts, ends = as_positions(sites)
    if len(starts) == 0:
        return []
    X = starts.reshape(-1, 1).astype(float)
    labels = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples,
                             cluster_selection_method=cluster_selection_method,
                             cluster_selection_epsilon=cluster_selection_epsilon).fit_predict(X)

    regions = []
    for lab in np.unique(labels):
        if lab == -1:
            continue
        mask = labels == lab
        start = int(starts[mask].min())
        end = int(ends[mask].max())
        n_sites = int(mask.sum())
        if end - start >= min_len and n_sites >= min_sites:
            regions.append((start, end, n_sites))
    regions.sort()
    return regions
