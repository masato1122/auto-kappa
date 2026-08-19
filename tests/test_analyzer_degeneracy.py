"""
Tests for auto_kappa/io/analyzer.py: get_average_at_degenerate_point().

The averaging must be applied once per degenerate group, in the same way as
ANPHON's Conductivity::average_self_energy_at_degenerate_point(), so that a
group is never re-processed from each of its members (which used to overwrite
the averaged values with partial averages).
"""
import numpy as np

from auto_kappa.io.analyzer import get_average_at_degenerate_point


def test_triply_degenerate_group_is_averaged_once():
    omega = np.array([[0.0, 0.0, 0.0, 120.0, 120.0, 300.0]])
    tau = np.array([[3.0, 5.0, 10.0, 2.0, 8.0, 1.0]])
    expected = np.array([[6.0, 6.0, 6.0, 5.0, 5.0, 1.0]])
    ave = get_average_at_degenerate_point(omega, tau)
    assert np.allclose(ave, expected)


def test_non_degenerate_values_are_kept():
    omega = np.array([[10.0, 20.0, 30.0]])
    tau = np.array([[1.0, 2.0, 3.0]])
    ave = get_average_at_degenerate_point(omega, tau)
    assert np.allclose(ave, tau)


def test_last_band_of_a_group_is_averaged():
    """The highest band must not be left unaveraged."""
    omega = np.array([[0.0, 100.0, 100.0]])
    tau = np.array([[1.0, 2.0, 4.0]])
    expected = np.array([[1.0, 3.0, 3.0]])
    ave = get_average_at_degenerate_point(omega, tau)
    assert np.allclose(ave, expected)


def test_all_bands_degenerate():
    omega = np.array([[50.0, 50.0, 50.0, 50.0]])
    tau = np.array([[1.0, 2.0, 3.0, 4.0]])
    ave = get_average_at_degenerate_point(omega, tau)
    assert np.allclose(ave, 2.5)


def test_multiple_kpoints_are_treated_independently():
    omega = np.array([[0.0, 0.0, 200.0],
                      [10.0, 20.0, 30.0]])
    tau = np.array([[1.0, 3.0, 7.0],
                    [1.0, 2.0, 3.0]])
    expected = np.array([[2.0, 2.0, 7.0],
                         [1.0, 2.0, 3.0]])
    ave = get_average_at_degenerate_point(omega, tau)
    assert np.allclose(ave, expected)


def test_sum_is_conserved():
    """Averaging redistributes values within a group without changing the sum."""
    rng = np.random.default_rng(0)
    omega = np.array([[0.0, 0.0, 0.0, 120.0, 120.0, 300.0]])
    tau = rng.random((1, 6))
    ave = get_average_at_degenerate_point(omega, tau)
    assert np.isclose(ave.sum(), tau.sum())


def _make_scattering(frequencies, gammas):
    """Build a minimal Scattering object without touching the file system."""
    from auto_kappa.io.scattering import Scattering
    scat = object.__new__(Scattering)
    scat.verbose = 0
    scat._temperature = 300.0
    scat.result = {'frequencies': np.array(frequencies),
                   'temperatures': np.array([300.0]),
                   'gammas': np.array([gammas])}
    scat.scattering_rates = {'phph': None, 'isotope': None, 'boundary': None}
    scat._total_scattering_rate = None
    scat._lifetime = None
    return scat


def test_phph_scattering_rate_is_averaged_over_degenerate_modes():
    """ANPHON averages the self-energy Gamma, not the lifetime."""
    from auto_kappa.io.scattering import convert_gamma2scatrate

    frequencies = [[0.0, 0.0, 0.0, 120.0, 120.0, 300.0]]
    gammas = [[1.0, 2.0, 3.0, 4.0, 6.0, 5.0]]  # Kayser
    scat = _make_scattering(frequencies, gammas)
    scat.set_scattering_rate_phph()

    expected = convert_gamma2scatrate(np.array([[2.0, 2.0, 2.0, 5.0, 5.0, 5.0]]))
    assert np.allclose(scat.scattering_rates['phph'], expected)


def test_lifetime_is_the_inverse_of_the_averaged_rate():
    """tau = 1/<rate>, not <1/rate>: the two differ for non-identical rates."""
    frequencies = [[10.0, 100.0, 100.0]]
    gammas = [[1.0, 1.0, 3.0]]
    scat = _make_scattering(frequencies, gammas)
    scat.set_scattering_rate_phph()
    scat._total_scattering_rate = scat.scattering_rates['phph']
    scat.set_lifetime()

    rate = scat.scattering_rates['phph'][0]
    assert np.allclose(scat.lifetime[0], 1.0 / rate)
    # the degenerate pair shares one value, different from the average lifetime
    assert np.isclose(scat.lifetime[0, 1], scat.lifetime[0, 2])
    raw_rate = rate[1] * np.array([0.5, 1.5])  # gammas 1 and 3 -> rates
    assert not np.isclose(scat.lifetime[0, 1], np.average(1.0 / raw_rate))


def test_group_split_by_the_eps_threshold():
    """Bands farther apart than eps (1e-3 cm^-1) belong to separate groups."""
    omega = np.array([[100.0, 100.0005, 100.1]])
    tau = np.array([[2.0, 4.0, 9.0]])
    ave = get_average_at_degenerate_point(omega, tau)
    assert np.allclose(ave[0, :2], 3.0)
    assert np.isclose(ave[0, 2], 9.0)


def test_top_band_is_not_zeroed():
    """The band loop must cover the highest band (it used to run to nb-1)."""
    omega = np.array([[1.0, 2.0, 3.0, 4.0]])
    tau = np.array([[10.0, 20.0, 30.0, 40.0]])
    ave = get_average_at_degenerate_point(omega, tau)
    assert ave[0, -1] != 0.0
    assert np.isclose(ave[0, -1], 40.0)


def test_matches_anphon_arithmetic_mean_of_gamma():
    """Reproduce ANPHON's average_over_degenerate_modes() (degeneracy_utils.h).

    ANPHON takes the plain arithmetic mean of the damping Gamma within each
    degenerate group, over every temperature column at once, and applies it to
    the ph-ph damping only -- before the isotope/boundary terms are added and
    before the conversion to a lifetime (conductivity.cpp:compute_kappa).
    """
    freqs = np.array([[0.0, 0.0, 0.0, 120.0, 120.0, 300.0],
                      [50.0, 50.0, 50.0, 50.0, 210.0, 400.0]])
    gammas = np.array([[1.7, 0.8, 2.4, 1.1, 2.9, 0.6],
                       [1.3, 2.2, 0.9, 1.8, 0.7, 2.6]])

    ref = gammas.copy()
    nk, ns = freqs.shape
    for i in range(nk):
        deg, omega_prev, ideg = [], freqs[i][0], 1
        for j in range(1, ns):
            if abs(freqs[i][j] - omega_prev) < 1e-3:
                ideg += 1
            else:
                deg.append(ideg)
                ideg = 1
                omega_prev = freqs[i][j]
        deg.append(ideg)
        im = 0
        for d in deg:
            if d > 1:
                data_sum = sum(ref[i][k] for k in range(im, im + d))
                for k in range(im, im + d):
                    ref[i][k] = data_sum / d
            im += d

    assert np.allclose(ref, get_average_at_degenerate_point(freqs, gammas))
