"""
Tests that chunked metric functions produce identical results to the original
non-chunked implementations.
"""
import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import UnitMatchPy.metric_functions as mf
import UnitMatchPy.default_params as default_params


def _make_test_data(n_units=50, spike_width=82, n_flips=2, seed=42):
    """Create synthetic avg_waveform_per_tp_flip and avg_centroid for testing."""
    rng = np.random.RandomState(seed)
    avg_waveform_per_tp_flip = rng.randn(3, n_units, spike_width, 2, n_flips)
    # Sprinkle some NaNs to test NaN handling
    nan_mask = rng.rand(3, n_units, spike_width, 2, n_flips) < 0.02
    avg_waveform_per_tp_flip[nan_mask] = np.nan
    avg_centroid = rng.randn(3, n_units, 2)
    param = default_params.get_default_param()
    param['n_units'] = n_units
    param['spike_width'] = spike_width
    param['n_channels'] = 384
    return avg_waveform_per_tp_flip, avg_centroid, param


def test_euclidean_metrics_match():
    """Chunked euclidean metrics must match the original three-step pipeline."""
    avg_wf_flip, _, param = _make_test_data()

    # Original pipeline
    euclid_dist_orig = mf.get_Euclidean_dist(avg_wf_flip, param)
    centroid_dist_orig, centroid_var_orig = mf.centroid_metrics(euclid_dist_orig, param)
    euclid_squeezed_orig = np.nanmin(
        euclid_dist_orig[:, param['peak_loc'] - param['waveidx'] == 0, :, :].squeeze(),
        axis=1)

    # Chunked pipeline (use small chunk to exercise multiple iterations)
    param['chunk_size'] = 17
    centroid_dist_new, centroid_var_new, euclid_squeezed_new = \
        mf.get_euclidean_metrics_chunked(avg_wf_flip, param)

    np.testing.assert_allclose(euclid_squeezed_new, euclid_squeezed_orig,
                               rtol=1e-10, equal_nan=True)
    np.testing.assert_allclose(centroid_dist_new, centroid_dist_orig,
                               rtol=1e-10, equal_nan=True)
    np.testing.assert_allclose(centroid_var_new, centroid_var_orig,
                               rtol=1e-10, equal_nan=True)


def test_recentered_metrics_match():
    """Chunked recentered metrics must match the original two-step pipeline."""
    avg_wf_flip, avg_centroid, param = _make_test_data()

    # Original pipeline
    euclid_dist_rc = mf.get_recentered_euclidean_dist(avg_wf_flip, avg_centroid, param)
    centroid_rc_orig = mf.recentered_metrics(euclid_dist_rc)

    # Chunked pipeline
    param['chunk_size'] = 17
    centroid_rc_new = mf.get_recentered_metrics_chunked(avg_wf_flip, avg_centroid, param)

    np.testing.assert_allclose(centroid_rc_new, centroid_rc_orig,
                               rtol=1e-10, equal_nan=True)


def test_chunk_size_edge_cases():
    """Chunked functions work with chunk_size larger than n_units and chunk_size=1."""
    avg_wf_flip, avg_centroid, param = _make_test_data(n_units=10)

    # chunk_size > n_units (single chunk)
    param['chunk_size'] = 9999
    cd1, cv1, ed1 = mf.get_euclidean_metrics_chunked(avg_wf_flip, param)
    rc1 = mf.get_recentered_metrics_chunked(avg_wf_flip, avg_centroid, param)

    # chunk_size = 1 (one unit per chunk)
    param['chunk_size'] = 1
    cd2, cv2, ed2 = mf.get_euclidean_metrics_chunked(avg_wf_flip, param)
    rc2 = mf.get_recentered_metrics_chunked(avg_wf_flip, avg_centroid, param)

    np.testing.assert_allclose(cd1, cd2, rtol=1e-10, equal_nan=True)
    np.testing.assert_allclose(cv1, cv2, rtol=1e-10, equal_nan=True)
    np.testing.assert_allclose(ed1, ed2, rtol=1e-10, equal_nan=True)
    np.testing.assert_allclose(rc1, rc2, rtol=1e-10, equal_nan=True)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
