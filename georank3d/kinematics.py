"""Markland kinematic screening, mean-based and dispersion-aware.

Manuscript Section 3.3.1.  A set is compatible with planar sliding when its
dip direction lies within ``planar_window_deg`` of the slope dip direction,
its dip exceeds the friction angle and its dip is less than the slope-face dip
(so that it daylights).  It is compatible with direct toppling when its dip
direction lies within ``toppling_window_deg`` of the anti-slope direction and
its dip exceeds ``90 - slope_dip + phi/2``.  Wedge sliding is evaluated
pairwise from the line of intersection (Eq. 11): the plunge must exceed the
friction angle, be less than the slope dip, and the trend must lie within
``wedge_window_deg`` of the slope dip direction.

Both the set means and the *full orientation population* of every set are
screened against identical thresholds; the population fractions are what the
manuscript calls dispersion-aware screening.  Oblique and flexural toppling
are not implemented, consistent with the scope limitation stated in the
manuscript.
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np

from .geometry import dip_from_normal, line_of_intersection, wrap180

__all__ = ["planar_toppling_flags", "kinematic_score", "wedge_pairs",
           "dispersion_screening", "wedge_dispersion"]


def planar_toppling_flags(dip, dipdir, cfg):
    """Vectorised planar-sliding and direct-toppling compatibility flags."""
    dip = np.atleast_1d(np.asarray(dip, dtype=float))
    dipdir = np.atleast_1d(np.asarray(dipdir, dtype=float))
    azi = np.abs(wrap180(dipdir - cfg.slope_dipdir_deg))
    opp = np.abs(wrap180(dipdir - (cfg.slope_dipdir_deg + 180.0) % 360.0))
    planar = (azi <= cfg.planar_window_deg) & (dip > cfg.friction_angle_deg) & (dip < cfg.slope_dip_deg)
    toppling = (opp <= cfg.toppling_window_deg) & (dip > (90.0 - cfg.slope_dip_deg + cfg.friction_angle_deg / 2.0))
    return planar, toppling, azi, opp


def kinematic_score(planar, toppling, cfg) -> np.ndarray:
    """K_s of manuscript Eq. (5): weighted, capped at 1."""
    score = cfg.planar_kin_weight * np.asarray(planar, dtype=float) \
        + cfg.toppling_kin_weight * np.asarray(toppling, dtype=float)
    return np.minimum(score, 1.0)


def wedge_pairs(set_records: List[Dict], cfg) -> List[Dict]:
    """Mean-orientation wedge screening for every pair of retained sets."""
    rows = []
    for i in range(len(set_records) - 1):
        for j in range(i + 1, len(set_records)):
            a, b = set_records[i], set_records[j]
            plunge, trend = line_of_intersection(a["axis"], b["axis"])
            if not np.isfinite(plunge):
                continue
            azi = abs(float(wrap180(trend - cfg.slope_dipdir_deg)))
            feasible = bool(plunge > cfg.friction_angle_deg
                            and plunge < cfg.slope_dip_deg
                            and azi <= cfg.wedge_window_deg)
            score = (1.0 - azi / cfg.wedge_window_deg) * min(plunge / cfg.slope_dip_deg, 1.0) if feasible else 0.0
            rows.append(dict(set_a=a["set_id"], set_b=b["set_id"],
                             intersection_plunge_deg=float(plunge),
                             intersection_trend_deg=float(trend),
                             slope_trend_difference_deg=azi,
                             wedge_feasible_mean=int(feasible),
                             wedge_score=float(max(score, 0.0))))
    rows.sort(key=lambda r: r["wedge_score"], reverse=True)
    return rows


def dispersion_screening(normals_by_set: Dict[int, np.ndarray], set_records: List[Dict], cfg) -> List[Dict]:
    """Population planar / toppling feasibility fractions per set."""
    rows = []
    for rec in set_records:
        N = normals_by_set[rec["set_id"]]
        if len(N) == 0:
            continue
        dip, dipdir = dip_from_normal(N)
        planar, toppling, _, _ = planar_toppling_flags(dip, dipdir, cfg)
        rows.append(dict(set_id=rec["set_id"], n_sampled=int(len(N)),
                         planar_feasible_fraction=float(planar.mean()),
                         toppling_feasible_fraction=float(toppling.mean()),
                         planar_feasible_mean=int(rec["planar_feasible_mean"]),
                         toppling_feasible_mean=int(rec["toppling_feasible_mean"]),
                         dip_p05=float(np.percentile(dip, 5)),
                         dip_p95=float(np.percentile(dip, 95))))
    return rows


def wedge_dispersion(normals_by_set: Dict[int, np.ndarray], mean_rows: List[Dict], cfg, rng) -> List[Dict]:
    """Population wedge feasibility fractions for every pair of sets.

    Intersections are formed from randomly paired members of the two
    populations, which samples the joint distribution of the two orientation
    populations under the assumption that they are independent.
    """
    lookup = {(r["set_a"], r["set_b"]): r["wedge_feasible_mean"] for r in mean_rows}
    ids = sorted(normals_by_set)
    rows = []
    for i in range(len(ids) - 1):
        for j in range(i + 1, len(ids)):
            A, B = normals_by_set[ids[i]], normals_by_set[ids[j]]
            if len(A) == 0 or len(B) == 0:
                continue
            m = min(len(A), len(B), cfg.dispersion_sample_cap)
            Aa = A[rng.integers(0, len(A), m)]
            Bb = B[rng.integers(0, len(B), m)]
            plunge, trend = line_of_intersection(Aa, Bb)
            ok = np.isfinite(plunge)
            plunge, trend = plunge[ok], trend[ok]
            if plunge.size == 0:
                continue
            azi = np.abs(wrap180(trend - cfg.slope_dipdir_deg))
            feasible = (plunge > cfg.friction_angle_deg) & (plunge < cfg.slope_dip_deg) \
                & (azi <= cfg.wedge_window_deg)
            rows.append(dict(set_a=ids[i], set_b=ids[j], n_sampled=int(plunge.size),
                             wedge_feasible_fraction=float(feasible.mean()),
                             wedge_feasible_mean=int(lookup.get((ids[i], ids[j]),
                                                                lookup.get((ids[j], ids[i]), 0))),
                             median_plunge_deg=float(np.median(plunge)),
                             p95_plunge_deg=float(np.percentile(plunge, 95))))
    return rows
