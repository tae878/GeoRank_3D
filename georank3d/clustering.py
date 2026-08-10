"""Two-stage clustering: axial density-peak pole finding, then spatial patches.

Manuscript Section 3.3.1 defines the discontinuity sets as the principal
orientation poles located as peaks of the normal-vector density on the sphere,
with each normal assigned to its nearest pole; a second, independent
density-based pass on the 3-D coordinates then splits each orientation cluster
into spatially coherent patches.

Two departures from the reference MATLAB release are deliberate:

*   Pole membership uses the axial angle ``arccos(|n.p|)`` rather than the
    Euclidean chord between hemisphere-flipped vectors.  The chord test has a
    seam at ``n_z = 0``; because the pole of an 80-degree plane plunges only
    10 degrees, that seam cuts through the middle of every steeply dipping
    set and discards the members that dip past vertical.
*   Set mean orientations come from the orientation-tensor principal
    eigenvector (:func:`georank3d.geometry.mean_axial_direction`) rather than
    the arithmetic mean of flipped unit vectors, which is not invariant to
    the arbitrary sign of each normal.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.spatial import cKDTree

from .geometry import dip_from_normal, mean_axial_direction

__all__ = ["find_orientation_poles", "assign_to_poles", "cluster_normals", "spatial_patches"]


def _axial_pairs(points: np.ndarray) -> np.ndarray:
    """Duplicate an axial sample with its antipodes for seam-free KD queries."""
    return np.vstack([points, -points])


def find_orientation_poles(normals: np.ndarray, cfg, rng) -> np.ndarray:
    """Greedy density-peak (mean-shift) search for principal orientation poles.

    Densities, mean-shift updates and support queries are all evaluated on the
    antipodally symmetrised sample, so a pole close to the equator collects its
    whole population instead of half of it.
    """
    n = len(normals)
    if n == 0:
        return np.zeros((0, 3))
    cap = min(n, int(cfg.cluster_sample_cap))
    sub = normals if cap >= n else normals[rng.choice(n, cap, replace=False)]

    h_chord = 2.0 * np.sin(np.radians(cfg.density_bandwidth_deg) / 2.0)
    doubled = _axial_pairs(sub)
    tree_all = cKDTree(doubled)
    m = len(sub)

    assigned = np.zeros(m, dtype=bool)
    poles: list[np.ndarray] = []
    first_peak = None
    min_support = max(30, int(round(cfg.min_cluster_points * m / n)))
    seed_cap = min(m, 20_000)

    for _ in range(int(cfg.max_sets) * 3):
        available = np.flatnonzero(~assigned)
        if len(available) < min_support:
            break
        seeds = available if len(available) <= seed_cap else \
            available[rng.choice(len(available), seed_cap, replace=False)]
        tree_avail = cKDTree(_axial_pairs(sub[available]))
        density = tree_avail.query_ball_point(sub[seeds], h_chord, workers=-1, return_length=True)
        p = sub[seeds[int(np.argmax(density))]].copy()

        for _ in range(60):                                # mean-shift to the mode
            nb = tree_all.query_ball_point(p, h_chord)
            if not nb:
                break
            block = doubled[nb]
            new = block.mean(axis=0)
            norm = np.linalg.norm(new)
            if norm < 1e-12:
                break
            new /= norm
            shift = np.degrees(np.arccos(min(abs(float(new @ p)), 1.0)))
            p = new
            if shift < 0.05:
                break

        peak = float(tree_all.query_ball_point(p, h_chord, return_length=True)) / 2.0
        if first_peak is None:
            first_peak = peak
        elif peak < cfg.pole_prominence_frac * first_peak:
            break

        within = np.degrees(np.arccos(np.clip(np.abs(sub @ p), 0.0, 1.0))) <= cfg.set_cone_deg
        support = np.flatnonzero(within & ~assigned)
        if len(support) < min_support:
            near = np.degrees(np.arccos(np.clip(np.abs(sub @ p), 0.0, 1.0))) <= cfg.density_bandwidth_deg
            assigned |= near
            continue

        merged = any(np.degrees(np.arccos(min(abs(float(q @ p)), 1.0))) < cfg.merge_deg for q in poles)
        if not merged:
            poles.append(p if p[2] >= 0 else -p)
        assigned[support] = True

    return np.array(poles) if poles else np.zeros((0, 3))


def assign_to_poles(normals: np.ndarray, poles: np.ndarray, cfg) -> np.ndarray:
    """Assign every normal to its nearest pole within the set cone (0 = noise)."""
    if len(poles) == 0:
        return np.zeros(len(normals), dtype=int)
    ang = np.degrees(np.arccos(np.clip(np.abs(normals @ poles.T), 0.0, 1.0)))
    nearest = np.argmin(ang, axis=1)
    labels = nearest + 1
    labels[ang[np.arange(len(normals)), nearest] > cfg.set_cone_deg] = 0
    return labels


def _drop_small_and_compact(labels: np.ndarray, min_points: int) -> np.ndarray:
    out = labels.copy()
    for lab in np.unique(out[out > 0]):
        if int((out == lab).sum()) < min_points:
            out[out == lab] = 0
    remap = {lab: i + 1 for i, lab in enumerate(np.unique(out[out > 0]))}
    return np.array([remap.get(v, 0) for v in out], dtype=int)


def cluster_normals(normals: np.ndarray, cfg, rng) -> Tuple[np.ndarray, np.ndarray]:
    """Full orientation-clustering stage.

    Returns ``(labels, poles)`` where ``labels`` is 1-based with 0 for
    unassigned points, and ``poles`` are the retained, re-estimated set axes
    ordered to match the labels.
    """
    poles = find_orientation_poles(normals, cfg, rng)
    if len(poles) == 0:
        return np.zeros(len(normals), dtype=int), poles
    labels = assign_to_poles(normals, poles, cfg)
    labels = _drop_small_and_compact(labels, cfg.min_cluster_points)

    refined = []
    for lab in np.unique(labels[labels > 0]):
        axis, _ = mean_axial_direction(normals[labels == lab])
        refined.append(axis)
    return labels, (np.array(refined) if refined else np.zeros((0, 3)))


def _dbscan(points: np.ndarray, eps: float, min_points: int) -> np.ndarray:
    """Minimal DBSCAN (0 = noise) so the package needs no clustering library."""
    n = len(points)
    tree = cKDTree(points)
    neighbours = tree.query_ball_point(points, eps, workers=-1)
    is_core = np.array([len(nb) >= min_points for nb in neighbours])
    labels = np.zeros(n, dtype=int)
    current = 0
    for i in range(n):
        if labels[i] != 0 or not is_core[i]:
            continue
        current += 1
        stack = [i]
        labels[i] = current
        while stack:
            j = stack.pop()
            if not is_core[j]:
                continue
            for nb in neighbours[j]:
                if labels[nb] == 0:
                    labels[nb] = current
                    stack.append(nb)
    return labels


def spatial_patches(points: np.ndarray, cfg, rng) -> np.ndarray:
    """Split one orientation cluster into spatially coherent patches.

    Above ``patch_sample_cap`` the patches are defined on a spatial subsample
    and propagated to every point by nearest core point, which keeps the patch
    extent (and therefore area and persistence) while staying tractable.
    """
    n = len(points)
    if n < cfg.min_patch_points:
        return np.zeros(n, dtype=int)
    cap = int(cfg.patch_sample_cap)
    if n <= cap:
        return _dbscan(points, cfg.patch_radius_m, cfg.min_patch_points)

    sel = rng.choice(n, cap, replace=False)
    min_sub = max(10, int(round(cfg.min_patch_points * cap / n)))
    sub_labels = _dbscan(points[sel], cfg.patch_radius_m, min_sub)
    if not (sub_labels > 0).any():
        return np.zeros(n, dtype=int)
    core = sel[sub_labels > 0]
    core_lab = sub_labels[sub_labels > 0]
    dist, nearest = cKDTree(points[core]).query(points, k=1, workers=-1)
    labels = core_lab[nearest]
    labels[dist > 3.0 * cfg.patch_radius_m] = 0
    remap = {lab: i + 1 for i, lab in enumerate(np.unique(labels[labels > 0]))}
    return np.array([remap.get(v, 0) for v in labels], dtype=int)


def set_orientation(normals: np.ndarray) -> tuple[float, float, np.ndarray]:
    """Convenience wrapper: dip, dip direction and axis of a set of normals."""
    axis, _ = mean_axial_direction(normals)
    dip, dipdir = dip_from_normal(axis)
    return dip, dipdir, axis
