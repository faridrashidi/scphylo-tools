"""Validate the bundled discrete distributions against SciPy."""

import math

import pytest
from scipy import stats

from scphylo.external import _betabinom as distribution


@pytest.mark.parametrize(
    ("i", "r", "shape1", "shape2"),
    [
        (0, 1, 0.5, 0.5),
        (1, 1, 0.5, 0.5),
        (5, 1, 0.5, 0.5),
        (0, 2, 2, 3),
        (20, 2, 2, 3),
        (5, 100, 0.1, 0.0001),
        (5, 100, 2, 3),
        (100, 5, 0.2, 10),
        (10, 2, 0.5, 30),
        (10, 30, 0.5, 2),
        (100, 30, 0.5, 40),
    ],
)
def test_beta_negative_binomial_matches_scipy(i, r, shape1, shape2):
    """Exercise stable central and tail algorithms against a trusted reference."""
    actual = (
        distribution.pmf_BetaNegativeBinomial(i, r, shape1, shape2),
        distribution.cdf_BetaNegativeBinomial(i, r, shape1, shape2),
        distribution.sf_BetaNegativeBinomial(i, r, shape1, shape2),
    )
    expected = (
        stats.betanbinom.pmf(i, r, shape1, shape2),
        stats.betanbinom.cdf(i, r, shape1, shape2),
        stats.betanbinom.sf(i, r, shape1, shape2),
    )

    assert actual == pytest.approx(expected, rel=1e-10, abs=1e-14)
    assert actual[1] + actual[2] == pytest.approx(1)


@pytest.mark.parametrize(
    ("i", "size", "shape1", "shape2"),
    [
        (-1, 10, 2, 3),
        (0, 10, 0.5, 0.5),
        (5, 10, 2, 3),
        (10, 10, 2, 3),
        (11, 10, 2, 3),
        (25, 50, 0.01, 40),
    ],
)
def test_beta_binomial_matches_scipy(i, size, shape1, shape2):
    """Check mass and both tails, including values outside the support."""
    assert distribution.pmf_BetaBinomial(i, size, shape1, shape2) == pytest.approx(
        stats.betabinom.pmf(i, size, shape1, shape2), rel=1e-10
    )
    assert distribution.cdf_BetaBinomial(i, size, shape1, shape2) == pytest.approx(
        stats.betabinom.cdf(i, size, shape1, shape2), rel=1e-10
    )
    assert distribution.sf_BetaBinomial(i, size, shape1, shape2) == pytest.approx(
        stats.betabinom.sf(i, size, shape1, shape2), rel=1e-10
    )


@pytest.mark.parametrize("observed", [-1, 0, 2, 5, 9])
def test_hypergeometric_matches_scipy(observed):
    """Validate support edges and the central hypergeometric calculation."""
    size, successes, population = 5, 6, 20
    assert distribution.pmf_hypergeometric(
        observed, size, successes, population
    ) == pytest.approx(stats.hypergeom.pmf(observed, population, successes, size))
    assert distribution.cdf_hypergeometric(
        observed, size, successes, population
    ) == pytest.approx(stats.hypergeom.cdf(observed, population, successes, size))
    assert distribution.sf_hypergeometric(
        observed, size, successes, population
    ) == pytest.approx(stats.hypergeom.sf(observed, population, successes, size))


@pytest.mark.parametrize("failures", [-1, 0, 2, 14, 15])
def test_negative_hypergeometric_matches_scipy(failures):
    """Validate the negative-hypergeometric support and cumulative functions."""
    required, successes, population = 2, 6, 20
    scipy_failures = population - successes
    assert distribution.pmf_neghypergeometric(
        failures, required, successes, population
    ) == pytest.approx(
        stats.nhypergeom.pmf(failures, population, scipy_failures, required)
    )
    assert distribution.cdf_neghypergeometric(
        failures, required, successes, population
    ) == pytest.approx(
        stats.nhypergeom.cdf(failures, population, scipy_failures, required)
    )
    assert distribution.sf_neghypergeometric(
        failures, required, successes, population
    ) == pytest.approx(
        stats.nhypergeom.sf(failures, population, scipy_failures, required)
    )


@pytest.mark.parametrize(
    ("function", "arguments"),
    [
        (distribution.pmf_BetaNegativeBinomial, (1, 0, 1, 1)),
        (distribution.cdf_BetaNegativeBinomial, (1, 1, 0, 1)),
        (distribution.sf_BetaNegativeBinomial, (1, 1, 1, 0)),
        (distribution.pmf_BetaBinomial, (1, -1, 1, 1)),
        (distribution.cdf_BetaBinomial, (1, 1, 0, 1)),
        (distribution.sf_BetaBinomial, (1, 1, 1, 0)),
        (distribution.pmf_hypergeometric, (1, 5, 21, 20)),
        (distribution.cdf_hypergeometric, (1, 21, 5, 20)),
        (distribution.sf_hypergeometric, (1, -1, 5, 20)),
        (distribution.pmf_neghypergeometric, (1, 0, 5, 20)),
        (distribution.cdf_neghypergeometric, (1, 6, 5, 20)),
        (distribution.sf_neghypergeometric, (1, 1, 21, 20)),
    ],
)
def test_distribution_parameters_are_validated(function, arguments):
    """Reject invalid shapes and population counts at public boundaries."""
    with pytest.raises(ValueError):
        function(*arguments)


def test_integral_policy_can_floor_values(monkeypatch):
    """Cover the documented opt-in policy for non-integral observations."""
    with pytest.raises(ValueError):
        distribution.pmf_BetaBinomial(1.5, 4, 2, 3)

    monkeypatch.setattr(distribution, "NonIntegralValuesAllowed_Others", True)

    assert distribution.pmf_BetaBinomial(1.5, 4.5, 2, 3) == pytest.approx(
        distribution.pmf_BetaBinomial(1, 4, 2, 3)
    )


@pytest.mark.parametrize(
    "value",
    [
        -0.9,
        -0.7,
        -0.65,
        -0.5,
        -0.2,
        0,
        0.5,
        1,
        1.2,
        1.56,
        2,
        2.2,
        2.5,
        3,
        3.5,
        4,
        5,
        6,
        10,
        20,
        50,
        999,
        1000,
    ],
)
def test_stirling_error_helpers_match_lgamma(value):
    """Exercise all approximation ranges against the defining expression."""
    expected = (
        math.lgamma(value + 1)
        - 0.5 * math.log(2 * math.pi)
        + value
        + 1
        - (value + 0.5) * math.log(value + 1)
    )

    assert distribution.logfbit(value) == pytest.approx(expected, abs=1e-12)
    assert distribution.logfbita(value) == pytest.approx(expected, abs=1e-12)


def test_stirling_derivatives_cover_ranges_and_domain_errors():
    """Cover recursive, asymptotic, and invalid derivative inputs."""
    for value in (-0.5, 0, 6.5, 7, 1e10):
        assert math.isfinite(distribution.logfbit2(value))
        assert math.isfinite(distribution.logfbit4(value))
    assert distribution.logfbit2dif(2) > 0
    assert distribution.logfbit4dif(2) > 0
    assert distribution.logfbit(1e8) > 0
    assert distribution.logfbita(1e8) > 0
    for function in (
        distribution.logfbit,
        distribution.logfbita,
        distribution.logfbit2,
        distribution.logfbit4,
    ):
        with pytest.raises(ValueError):
            function(-1)


@pytest.mark.parametrize("value", [1e-16, -0.005, 0.5, -0.9, 2, 5])
def test_accurate_log_helpers(value):
    """Check small-value and ordinary paths in the logarithm helpers."""
    assert distribution.log0(value) == pytest.approx(math.log1p(value))
    assert distribution.log1(value) == pytest.approx(math.log1p(value) - value)


@pytest.mark.parametrize(
    ("a", "b"),
    [
        (-1, 10),
        (2, 10),
        (2, 3),
        (2, 2),
        (2, -0.7),
        (1, -0.7),
        (0.5, -0.7),
        (1e-6, -0.2),
        (2, -0.2),
        (0.1, 0),
    ],
)
def test_accurate_stirling_difference(a, b):
    """Exercise every stable difference strategy against direct evaluation."""
    assert distribution.lfbaccdif1(a, b) == pytest.approx(
        distribution.logfbit(b) - distribution.logfbit(a + b), abs=1e-14
    )


@pytest.mark.parametrize(
    ("a", "b", "c", "d"),
    [
        (0, 2, 3, 4),
        (4, 5, 2, 3),
        (2, 3, 4, 5),
        (1e100, 1e100, 1, 1),
        (1, 1, 1e100, 1e100),
    ],
)
def test_accurate_product_difference(a, b, c, d):
    """Cover both exponent orderings, including widely separated products."""
    assert distribution.Generalabminuscd(a, b, c, d) == pytest.approx(a * b - c * d)
