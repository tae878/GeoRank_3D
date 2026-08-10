"""Statistical characterisation of the detachable-block inventory.

Reproduces manuscript Section 4.4 and Figures 16 to 18 from the block
inventory of Table 5:

*   descriptive statistics and sphere-equivalent diameters,
*   candidate distribution fitting (log-normal, Weibull, gamma, exponential,
    normal) ranked by AIC and BIC, with a Kolmogorov-Smirnov test and a
    Shapiro-Wilk test on ln V and the normal quantile-quantile R-squared,
*   the magnitude-frequency relation N(>=V) and the volume-weighted
    exceedance curve,
*   the Gini coefficient of the volume distribution, and
*   sector-wise volume shares with bootstrap confidence intervals and a
    Kruskal-Wallis test on median block size.

None of this analysis existed in the submitted code, although every number it
produces is reported in the manuscript.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
from scipy import stats

__all__ = ["read_block_inventory", "describe", "fit_distributions",
           "magnitude_frequency", "gini_coefficient", "sector_analysis",
           "analyse_blocks"]

_CANDIDATES = {
    "lognormal": (stats.lognorm, dict(floc=0), 2),
    "weibull": (stats.weibull_min, dict(floc=0), 2),
    "gamma": (stats.gamma, dict(floc=0), 2),
    "exponential": (stats.expon, dict(floc=0), 1),
    "normal": (stats.norm, {}, 2),
}


def read_block_inventory(path) -> List[Dict]:
    """Read a block inventory CSV with at least ``volume_m3``.

    Recognised optional columns: ``block_id`` and ``sector``.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"block inventory not found: {p}")
    with p.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"{p} contains no data rows")
    headers = list(rows[0].keys())
    lowered = {h: h.strip().lower() for h in headers}
    vol_col = next((h for h in headers if "vol" in lowered[h]), None)
    if vol_col is None:
        raise ValueError(f"{p}: no volume column found (looked for a header containing 'vol'); "
                         f"columns present: {headers}")
    sec_col = next((h for h in headers if "sector" in lowered[h]), None)
    id_col = next((h for h in headers if lowered[h] in ("block_id", "blockid", "id")), None)

    out = []
    for i, row in enumerate(rows):
        try:
            v = float(row[vol_col])
        except (TypeError, ValueError):
            continue
        if not np.isfinite(v) or v <= 0:
            raise ValueError(f"{p} row {i + 2}: volume must be finite and positive (got {row[vol_col]!r})")
        out.append(dict(block_id=(row[id_col] if id_col else str(i + 1)),
                        sector=(row[sec_col].strip() if sec_col and row[sec_col] else "unassigned"),
                        volume_m3=v))
    if not out:
        raise ValueError(f"{p}: no valid volumes parsed")
    return out


def describe(volumes: Sequence[float]) -> Dict:
    V = np.asarray(volumes, dtype=float)
    return dict(count=int(V.size), sum_m3=float(V.sum()), mean_m3=float(V.mean()),
                median_m3=float(np.median(V)), min_m3=float(V.min()), max_m3=float(V.max()),
                std_m3=float(V.std(ddof=1)),
                coefficient_of_variation=float(V.std(ddof=1) / V.mean()),
                equivalent_diameter_m=[float(x) for x in (6.0 * V / np.pi) ** (1.0 / 3.0)])


def fit_distributions(volumes: Sequence[float]) -> Dict:
    """Fit the five candidate models and rank them by AIC and BIC."""
    V = np.asarray(volumes, dtype=float)
    n = V.size
    models = []
    for name, (dist, kwargs, k) in _CANDIDATES.items():
        params = dist.fit(V, **kwargs)
        loglik = float(np.sum(dist.logpdf(V, *params)))
        models.append(dict(model=name, params=[float(x) for x in params], k=k,
                           loglik=loglik, aic=float(2 * k - 2 * loglik),
                           bic=float(k * np.log(n) - 2 * loglik)))
    models.sort(key=lambda m: m["aic"])

    lnV = np.log(V)
    mu, sigma = float(lnV.mean()), float(lnV.std(ddof=0))
    ks = stats.kstest(lnV, "norm", args=(mu, sigma))
    sw = stats.shapiro(lnV)
    osm, osr = stats.probplot(lnV, dist="norm")[0]
    return dict(models=models, best_by_aic=models[0]["model"],
                best_by_bic=min(models, key=lambda m: m["bic"])["model"],
                lognormal=dict(ln_mean=mu, ln_std=sigma, geometric_mean_m3=float(np.exp(mu)),
                               ks_statistic=float(ks.statistic), ks_pvalue=float(ks.pvalue),
                               shapiro_w=float(sw.statistic), shapiro_pvalue=float(sw.pvalue),
                               qq_r_squared=float(np.corrcoef(osm, osr)[0, 1] ** 2)))


def magnitude_frequency(volumes: Sequence[float], thresholds: Sequence[float] = (1.0, 5.0, 20.0)) -> Dict:
    V = np.sort(np.asarray(volumes, dtype=float))[::-1]
    total = V.sum()
    return dict(
        exceedance=[dict(volume_m3=float(v), n_ge=int(i + 1),
                         volume_share_ge=float(V[:i + 1].sum() / total))
                    for i, v in enumerate(V)],
        thresholds=[dict(threshold_m3=float(t), n_ge=int((V >= t).sum()),
                         volume_share_percent=float(100.0 * V[V >= t].sum() / total))
                    for t in thresholds])


def gini_coefficient(volumes: Sequence[float]) -> float:
    V = np.sort(np.asarray(volumes, dtype=float))
    n = V.size
    idx = np.arange(1, n + 1)
    return float((2.0 * (idx * V).sum()) / (n * V.sum()) - (n + 1) / n)


def sector_analysis(records: List[Dict], n_bootstrap: int = 10_000, seed: int = 42) -> Dict:
    """Sector volume shares with bootstrap CIs, plus a Kruskal-Wallis test."""
    rng = np.random.default_rng(seed)
    V = np.array([r["volume_m3"] for r in records], dtype=float)
    sectors = np.array([r["sector"] for r in records])
    total = V.sum()
    names = sorted(set(sectors))

    out = []
    for s in names:
        m = sectors == s
        boot = np.empty(n_bootstrap)
        n = V.size
        for b in range(n_bootstrap):
            pick = rng.integers(0, n, n)
            vb, sb = V[pick], sectors[pick]
            boot[b] = 100.0 * vb[sb == s].sum() / vb.sum()
        out.append(dict(sector=s, n_blocks=int(m.sum()), total_volume_m3=float(V[m].sum()),
                        volume_share_percent=float(100.0 * V[m].sum() / total),
                        share_ci95_percent=[float(np.percentile(boot, 2.5)),
                                            float(np.percentile(boot, 97.5))],
                        median_volume_m3=float(np.median(V[m])),
                        min_volume_m3=float(V[m].min()), max_volume_m3=float(V[m].max()),
                        n_blocks_ge_20m3=int((V[m] >= 20).sum())))
    out.sort(key=lambda r: r["volume_share_percent"], reverse=True)

    groups = [V[sectors == s] for s in names if (sectors == s).sum() > 0]
    kw = stats.kruskal(*groups) if len(groups) > 1 else None
    return dict(sectors=out,
                kruskal_h=(float(kw.statistic) if kw else None),
                kruskal_pvalue=(float(kw.pvalue) if kw else None))


def classify_volume(v: float) -> str:
    if v < 1.0:
        return "small"
    if v < 5.0:
        return "medium"
    if v < 20.0:
        return "large"
    return "extra-large"


def analyse_blocks(records: List[Dict], n_bootstrap: int = 10_000, seed: int = 42) -> Dict:
    V = [r["volume_m3"] for r in records]
    for r, d in zip(records, describe(V)["equivalent_diameter_m"]):
        r["equivalent_diameter_m"] = d
        r["volume_class"] = classify_volume(r["volume_m3"])
    return dict(inventory=records, descriptive=describe(V), distributions=fit_distributions(V),
                magnitude_frequency=magnitude_frequency(V), gini=gini_coefficient(V),
                sectors=sector_analysis(records, n_bootstrap=n_bootstrap, seed=seed))
