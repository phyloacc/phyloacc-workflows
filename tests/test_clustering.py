#############################################################################
# Unit tests for lib/clustering.py - the conserved-region clustering methods
# (windowed, hmm, hdbscan_cluster). Pure-Python/numpy; hdbscan tests skip if the
# hdbscan package isn't installed (it's not in the base test/CI env).
#############################################################################

import importlib.util

import pytest

import lib.clustering as C


def test_as_positions_ints_and_tuples_sorted():
    starts, ends = C.as_positions([(7, 8), 3, (0, 2)])
    assert list(starts) == [0, 3, 7]
    assert list(ends) == [2, 4, 8]  # scalar 3 -> end 4; scalar 7 -> 8


# --- windowed -------------------------------------------------------------

def test_windowed_threshold_and_merge():
    # window_bp=10, min 3 sites/window, length 70 -> 7 windows.
    # window 0 (0-9): 3 sites -> conserved; window 1: 2 sites -> not;
    # windows 5 and 6: 3 sites each -> conserved and adjacent -> merged.
    sites = [0, 3, 7, 12, 15, 50, 51, 52, 60, 61, 63]
    regions = C.windowed(sites, length=70, window_bp=10, min_sites_per_window=3)
    assert regions == [(0, 10, 3), (50, 70, 6)]


def test_windowed_clips_last_window_to_length():
    # length not a multiple of window_bp -> final region end clipped to length.
    sites = [50, 51, 52]
    regions = C.windowed(sites, length=55, window_bp=10, min_sites_per_window=3)
    assert regions == [(50, 55, 3)]


def test_windowed_empty():
    assert C.windowed([], length=100, window_bp=10, min_sites_per_window=3) == []


# --- hmm ------------------------------------------------------------------

def test_hmm_empty_sites_no_regions():
    assert C.hmm([], length=1000) == []


def test_hmm_recovers_a_dense_block():
    # A contiguous run of conserved sites should be recovered as (roughly) one region,
    # with essentially no coverage outside it. Use permissive-but-sane probs.
    block = list(range(100, 160))  # 60 contiguous conserved positions
    regions = C.hmm(block, length=300, t0_0=0.9, t1_1=0.95, e0_0=0.9, e1_1=0.9,
                    s0=0.9, min_len=5, max_len=1000)
    assert len(regions) >= 1
    covered = set()
    for s, e, n in regions:
        covered.update(range(s, e))
        assert e - s > 5              # min_len filter (strict) respected
        assert n == sum(1 for p in range(s, e) if 100 <= p < 160)  # n_sites counts real sites
    # most of the block is covered, and coverage stays near the block
    assert len(covered & set(range(100, 160))) >= 40
    assert not (covered & set(range(0, 90)))
    assert not (covered & set(range(180, 300)))


# --- hdbscan (optional dependency) ----------------------------------------

_HAS_HDBSCAN = importlib.util.find_spec("hdbscan") is not None


@pytest.mark.skipif(not _HAS_HDBSCAN, reason="hdbscan not installed")
def test_hdbscan_two_separated_clusters():
    sites = list(range(0, 30)) + list(range(1000, 1030))
    regions = C.hdbscan_cluster(sites, min_cluster_size=5)
    assert len(regions) == 2
    (s0, e0, n0), (s1, e1, n1) = regions           # returned sorted by start
    assert s0 == 0 and e0 == 30 and n0 == 30
    assert s1 == 1000 and e1 == 1030 and n1 == 30


@pytest.mark.skipif(not _HAS_HDBSCAN, reason="hdbscan not installed")
def test_hdbscan_empty():
    assert C.hdbscan_cluster([], min_cluster_size=5) == []
