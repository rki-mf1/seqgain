"""Accessible-fraction gamma-Poisson scenario model for depth breadth."""
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import nbinom


def effort_grid():
    # Fine resolution near measured effort; bounded report size at long horizons.
    return np.concatenate((np.linspace(1, 2, 101), np.linspace(2.25, 100, 392)))


def breadth(effort, threshold, accessible, mean, shape):
    effort = np.asarray(effort, dtype=float)
    return accessible * nbinom.sf(threshold - 1, shape, shape / (shape + mean * effort))


def _fit_parameters(x, depth, y):
    def residual(p):
        return breadth(x, depth, *p) - y
    candidates = [least_squares(residual, [max(.1, float(max(y))), 10, start],
                                bounds=([.000001, .000001, .03], [1, 1e6, 1e4]), max_nfev=1500)
                  for start in (.2, 2, 100)]
    return min(candidates, key=lambda result: np.sum(result.fun**2))


def _conditional_prediction(histogram, thresholds, efforts, accessible, mean, shape):
    """Posterior-predict future depth conditional on the measured depth histogram."""
    total = sum(histogram.values())
    threshold_array = np.asarray(thresholds, dtype=int)
    beta = shape / mean
    zero_likelihood = (beta / (beta + 1)) ** shape
    zero_denominator = (1 - accessible) + accessible * zero_likelihood
    zero_accessible = accessible * zero_likelihood / zero_denominator if zero_denominator else 0
    observed = np.array([sum(count for depth, count in histogram.items() if depth >= threshold) / total
                         for threshold in threshold_array])
    matrix = np.empty((len(threshold_array), len(efforts)))
    for index, effort in enumerate(efforts):
        if effort <= 1:
            matrix[:, index] = observed
            continue
        probability = (beta + 1) / (beta + effort)
        predicted = np.zeros(len(threshold_array))
        for depth, count in histogram.items():
            needed = threshold_array - depth
            chance = np.where(needed <= 0, 1.0,
                              nbinom.sf(needed - 1, shape + depth, probability))
            if depth == 0:
                chance = np.where(needed <= 0, 1.0, chance * zero_accessible)
            predicted += count * chance
        matrix[:, index] = np.clip(np.maximum(observed, predicted / total), 0, 1)
    return {str(threshold): [float(value) for value in matrix[index]]
            for index, threshold in enumerate(threshold_array)}


def _retrospective_backtests(rows, thresholds, histograms):
    fractions = sorted({r["fraction"] for r in rows})
    available = {round(f, 12) for f in fractions}
    tests = []
    for anchor in fractions:
        validation = round(2 * anchor, 12)
        training = [f for f in fractions if f <= anchor]
        if validation not in available or len(training) < 4:
            continue
        x = np.array([f/anchor for f in training for _ in thresholds])
        depth = np.array(thresholds * len(training))
        y = np.array([np.mean([r["breadth"][str(d)] for r in rows if r["fraction"] == f])
                      for f in training for d in thresholds])
        if not np.any(y):
            continue
        result = _fit_parameters(x, depth, y)
        effort = validation / anchor
        replicate_predictions = [
            _conditional_prediction(histograms[(replicate, anchor)], thresholds, [effort], *result.x)
            for replicate in sorted({r["replicate"] for r in rows})
        ]
        predicted = np.array([
            np.mean([prediction[str(threshold)][0] for prediction in replicate_predictions])
            for threshold in thresholds
        ])
        observed = np.array([np.mean([r["breadth"][str(d)] for r in rows
                                     if round(r["fraction"], 12) == validation]) for d in thresholds])
        tests.append({"training_fraction": float(anchor), "validation_fraction": float(validation),
                      "max_absolute_error": float(np.max(np.abs(predicted-observed)))})
    return tests


def fit_breadth(rows, thresholds, full_depth_histogram,
                histograms, forecast_thresholds):
    """Fit a depth prior, then forecast conditionally from the measured histogram."""
    fractions = sorted({r["fraction"] for r in rows})
    x = np.array([f for f in fractions for _ in thresholds])
    depth = np.array(thresholds * len(fractions))
    y = np.array([np.mean([r["breadth"][str(d)] for r in rows if r["fraction"] == f])
                  for f in fractions for d in thresholds])
    if not np.any(y):
        return {"status": "no_coverage", "forecast_efforts": [],
                "forecast_breadth": {},
                "retrospective_backtests": [],
                "warnings": ["No usable target coverage; additional depth cannot be estimated."]}

    result = _fit_parameters(x, depth, y)
    a, mean, shape = map(float, result.x)
    rmse = float(np.sqrt(np.mean(result.fun**2)))
    holdout = _fit_parameters(x[x < 1], depth[x < 1], y[x < 1])
    holdout_error = float(np.max(np.abs(breadth(1, np.array(thresholds), *holdout.x) - y[x == 1])))
    backtests = _retrospective_backtests(rows, thresholds, histograms)
    warnings = []
    if not result.success or rmse > .03 or holdout_error > .05:
        warnings.append("Poor fit or >5 percentage point holdout error: forecast is unreliable.")
    if np.linalg.cond(result.jac) > 1e6 or a < .00001 or shape < .031 or shape > 9999:
        warnings.append("Model parameters are weakly identified or at a bound; asymptote is uncertain.")
    if backtests and max(test["max_absolute_error"] for test in backtests) > .05:
        warnings.append("A retrospective 2× breadth forecast missed an observed threshold by more than 5 percentage points.")
    grid = effort_grid()
    forecasts = _conditional_prediction(full_depth_histogram, forecast_thresholds, grid,
                                        a, mean, shape)
    return {"status": "caution" if warnings else "fit", "accessible_fraction": a,
            "accessible_mean_depth": mean, "dispersion": shape, "rmse": rmse,
            "holdout_max_absolute_error": holdout_error,
            "forecast_efforts": [float(t) for t in grid],
            "forecast_breadth": forecasts,
            "retrospective_backtests": backtests,
            "warnings": warnings}
