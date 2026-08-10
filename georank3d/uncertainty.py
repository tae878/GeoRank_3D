"""Monte Carlo cluster-tracking orientation uncertainty (Section 3.3.2).

Each run draws a random subsample of the retained normals, rotates every
normal by a random angle of order ``mc_rotation_sigma_deg`` in a random
tangent direction, jitters the angular tolerance and the minimum cluster size,
re-runs the orientation clustering and matches the resulting clusters back to
the reference sets by minimum axial deviation.  Matches beyond
``mc_assignment_threshold_deg`` are rejected.

Two changes relative to the reference MATLAB release:

*   the assignment threshold is read from the configuration instead of being a
    literal, and
*   the reported interval is a percentile interval of the *run-to-run
    orientation*, together with the Fisher (1993) cone of 95 per cent
    confidence about the mean axis.  A standard error of the mean across runs,
    as previously reported, shrinks with the number of Monte Carlo runs and
    therefore describes the precision of the averaging rather than the
    stability of the orientation.
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np

from .clustering import cluster_normals
from .geometry import dip_from_normal, mean_axial_direction, normal_from_dip

__all__ = ["monte_carlo_uncertainty", "perturb_normals", "fisher_cone_95"]


def perturb_normals(normals: np.ndarray, sigma_deg: float, rng) -> np.ndarray:
    """Rotate each normal by a N(0, sigma) angle in a random tangent direction."""
    n = len(normals)
    v = rng.standard_normal((n, 3))
    v -= (v * normals).sum(axis=1, keepdims=True) * normals
    norms = np.linalg.norm(v, axis=1, keepdims=True)
    ok = norms[:, 0] > 1e-12
    v[ok] /= norms[ok]
    angles = np.radians(rng.standard_normal(n) * sigma_deg)[:, None]
    out = np.cos(angles) * normals + np.sin(angles) * v
    out[~ok] = normals[~ok]
    return out / np.maximum(np.linalg.norm(out, axis=1, keepdims=True), 1e-300)


def fisher_cone_95(axes: np.ndarray) -> float:
    """Semi-apical angle of the 95 per cent confidence cone (Fisher 1993).

    Uses the standard Fisher estimate ``alpha_95`` for axial data reduced to a
    common hemisphere.  Returns NaN when fewer than three axes are supplied or
    when the resultant is degenerate.
    """
    A = np.asarray(axes, dtype=float)
    if len(A) < 3:
        return float("nan")
    axis, _ = mean_axial_direction(A)
    A = np.where((A @ axis < 0)[:, None], -A, A)
    n = len(A)
    R = float(np.linalg.norm(A.sum(axis=0)))
    if R <= 0 or n - R <= 0:
        return float("nan")
    ratio = (n - R) / R
    value = 1.0 - ratio * ((1.0 / 0.05) ** (1.0 / (n - 1)) - 1.0)
    if not -1.0 <= value <= 1.0:
        return float("nan")
    return float(np.degrees(np.arccos(value)))


def monte_carlo_uncertainty(normals: np.ndarray, reference_sets: List[Dict], cfg, rng) -> Dict:
    """Run the Monte Carlo and summarise per-set orientation stability."""
    runs: List[Dict] = []
    per_set: Dict[int, List[np.ndarray]] = {rec["set_id"]: [] for rec in reference_sets}
    if cfg.monte_carlo_runs <= 0 or len(normals) < cfg.min_cluster_points or not reference_sets:
        return dict(runs=runs, per_set=[], n_runs=0)

    ref_axes = {rec["set_id"]: normal_from_dip(rec["dip_deg"], rec["dipdir_deg"])
                for rec in reference_sets}
    sub_n = min(cfg.monte_carlo_subsample, len(normals))

    for r in range(int(cfg.monte_carlo_runs)):
        sel = rng.choice(len(normals), sub_n, replace=False)
        Ns = perturb_normals(normals[sel], cfg.mc_rotation_sigma_deg, rng)

        run_cfg = _jittered(cfg, rng)
        labels, poles = cluster_normals(Ns, run_cfg, rng)
        found = []
        for lab in np.unique(labels[labels > 0]):
            axis, _ = mean_axial_direction(Ns[labels == lab])
            found.append(axis)

        matched = 0
        for sid, ref_axis in ref_axes.items():
            if not found:
                continue
            devs = [np.degrees(np.arccos(min(abs(float(ref_axis @ a)), 1.0))) for a in found]
            best = int(np.argmin(devs))
            if devs[best] <= cfg.mc_assignment_threshold_deg:
                per_set[sid].append(found[best])
                matched += 1
        dips = [dip_from_normal(a)[0] for a in found]
        runs.append(dict(run=r + 1,
                         density_bandwidth_deg=float(run_cfg.density_bandwidth_deg),
                         min_cluster_points=int(run_cfg.min_cluster_points),
                         n_clusters=len(found),
                         mean_cluster_dip_deg=float(np.mean(dips)) if dips else float("nan"),
                         sets_matched=matched))

    summary = []
    for rec in reference_sets:
        axes = np.array(per_set[rec["set_id"]]) if per_set[rec["set_id"]] else np.zeros((0, 3))
        row = dict(set_id=rec["set_id"], reference_dip_deg=rec["dip_deg"],
                   reference_dipdir_deg=rec["dipdir_deg"],
                   runs_matched=len(axes), runs_total=int(cfg.monte_carlo_runs))
        if len(axes) >= 3:
            dips, dipdirs = dip_from_normal(axes)
            mean_axis, _ = mean_axial_direction(axes)
            mdip, mdir = dip_from_normal(mean_axis)
            centred = (dipdirs - mdir + 180.0) % 360.0 - 180.0
            row.update(mc_mean_dip_deg=float(mdip), mc_mean_dipdir_deg=float(mdir),
                       dip_p2_5_deg=float(np.percentile(dips, 2.5)),
                       dip_p97_5_deg=float(np.percentile(dips, 97.5)),
                       dipdir_p2_5_deg=float((mdir + np.percentile(centred, 2.5)) % 360.0),
                       dipdir_p97_5_deg=float((mdir + np.percentile(centred, 97.5)) % 360.0),
                       fisher_alpha95_deg=fisher_cone_95(axes))
        else:
            row.update(mc_mean_dip_deg=float("nan"), mc_mean_dipdir_deg=float("nan"),
                       dip_p2_5_deg=float("nan"), dip_p97_5_deg=float("nan"),
                       dipdir_p2_5_deg=float("nan"), dipdir_p97_5_deg=float("nan"),
                       fisher_alpha95_deg=float("nan"))
        summary.append(row)
    return dict(runs=runs, per_set=summary, n_runs=int(cfg.monte_carlo_runs))


def _jittered(cfg, rng):
    from copy import copy
    run_cfg = copy(cfg)
    run_cfg.density_bandwidth_deg = max(
        4.0, cfg.density_bandwidth_deg + cfg.mc_angular_tolerance_sigma_deg * rng.standard_normal())
    jitter = 1.0 + cfg.mc_min_cluster_jitter * (2.0 * rng.random() - 1.0)
    run_cfg.min_cluster_points = max(30, int(round(cfg.min_cluster_points * jitter)))
    run_cfg.cluster_sample_cap = min(cfg.cluster_sample_cap, cfg.monte_carlo_subsample)
    return run_cfg
