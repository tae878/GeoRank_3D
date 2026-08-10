"""Orientation conventions and plane geometry for GeoRank_3D.

Coordinate frame: X = east, Y = north, Z = up (right-handed).

A discontinuity plane is described by its dip ``delta`` and dip direction
``alpha``.  Its *upward* unit normal is

    n = (sin(delta) sin(alpha), sin(delta) cos(alpha), cos(delta))

so that

    delta = arccos(|n_z|)                      (manuscript Eq. 8)
    alpha = mod(atan2(n_x, n_y), 360)

The horizontal projection of the **upward** normal points *down dip*: for a
surface descending eastward, ``z = -tan(delta) x``, the upward normal is
``(sin delta, 0, cos delta)``, which points east.  It is the **downward**
(lower-hemisphere) pole whose trend is ``alpha + 180``.

Manuscript Eq. (9) gives ``alpha = mod(atan2(-a, -b), 360)`` for a normal that
Section 3.3.3 explicitly orients "into the upper hemisphere", and
``GeoRank_3D_orientation.m`` and ``GeoRank_3D_dipdir2normal.m`` implement it
that way.  That combination returns the dip direction rotated by 180 degrees;
the minus signs belong with a downward-oriented normal.  See
``AUDIT_REPORT.md`` for the consequences for Table 8 of the manuscript.

Normals are treated as *axial* data throughout: n and -n denote the same
plane.  Angular separations use arccos(|n_i . n_j|) (manuscript Eq. 4) and
mean orientations come from the principal eigenvector of the orientation
tensor rather than from an arithmetic mean of hemisphere-flipped vectors.
Hemisphere flipping introduces a seam at n_z = 0 that truncates precisely the
steeply dipping sets this study is concerned with.

Normals are treated as *axial* data throughout: n and -n denote the same
plane.  Angular separations use arccos(|n_i . n_j|) (manuscript Eq. 4) and
mean orientations come from the principal eigenvector of the orientation
tensor rather than from an arithmetic mean of hemisphere-flipped vectors.
Hemisphere flipping introduces a seam at n_z = 0 that truncates precisely the
steeply dipping sets this study is concerned with.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "normal_from_dip",
    "dip_from_normal",
    "angular_difference_deg",
    "mean_axial_direction",
    "orientation_tensor",
    "fit_plane",
    "projected_hull",
    "line_of_intersection",
    "wrap180",
]


def wrap180(angles_deg):
    """Wrap angular differences to the interval (-180, 180]."""
    a = np.asarray(angles_deg, dtype=float)
    return np.degrees(np.arctan2(np.sin(np.radians(a)), np.cos(np.radians(a))))


def normal_from_dip(dip_deg, dipdir_deg):
    """Upward unit normal(s) from dip / dip direction in degrees."""
    d = np.radians(np.atleast_1d(np.asarray(dip_deg, dtype=float)))
    a = np.radians(np.atleast_1d(np.asarray(dipdir_deg, dtype=float)))
    n = np.stack([np.sin(d) * np.sin(a), np.sin(d) * np.cos(a), np.cos(d)], axis=-1)
    n /= np.maximum(np.linalg.norm(n, axis=-1, keepdims=True), 1e-300)
    return n[0] if np.ndim(dip_deg) == 0 else n


def dip_from_normal(normals):
    """Dip and dip direction in degrees from plane normal(s).

    Antipodal input is accepted: the normal is interpreted axially, so n and
    -n give the same answer.  The dip follows manuscript Eq. (8); the dip
    direction is ``atan2(n_x, n_y)`` of the upward normal, which is the
    corrected form of manuscript Eq. (9) (see the module docstring).
    """
    n = np.atleast_2d(np.asarray(normals, dtype=float))
    n = n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-300)
    n = np.where((n[:, 2] < 0)[:, None], -n, n)          # axial -> upward rep.
    dip = np.degrees(np.arccos(np.clip(n[:, 2], -1.0, 1.0)))
    dipdir = np.degrees(np.arctan2(n[:, 0], n[:, 1])) % 360.0
    if np.ndim(normals) == 1:
        return float(dip[0]), float(dipdir[0])
    return dip, dipdir


def angular_difference_deg(a, b):
    """Axial angle between unit vectors, manuscript Eq. (4): arccos(|a.b|)."""
    a = np.atleast_2d(np.asarray(a, dtype=float))
    b = np.atleast_2d(np.asarray(b, dtype=float))
    a = a / np.maximum(np.linalg.norm(a, axis=1, keepdims=True), 1e-300)
    b = b / np.maximum(np.linalg.norm(b, axis=1, keepdims=True), 1e-300)
    dots = np.clip(np.abs(a @ b.T), 0.0, 1.0)
    out = np.degrees(np.arccos(dots))
    return out


def orientation_tensor(normals):
    """Normalised orientation tensor T = (1/n) sum n_i n_i^T (Fisher 1993)."""
    n = np.asarray(normals, dtype=float)
    n = n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-300)
    return (n.T @ n) / max(len(n), 1)


def mean_axial_direction(normals):
    """Principal axis of a set of axial unit vectors, plus its eigenvalues.

    Returns ``(axis, eigenvalues_descending)``.  The axis is the eigenvector
    of the orientation tensor with the largest eigenvalue, oriented into the
    upper hemisphere.  Unlike the arithmetic mean of hemisphere-flipped
    normals, this is invariant to the arbitrary sign of each input vector.
    """
    n = np.asarray(normals, dtype=float)
    if len(n) == 0:
        return np.array([0.0, 0.0, 1.0]), np.zeros(3)
    T = orientation_tensor(n)
    vals, vecs = np.linalg.eigh(T)
    order = np.argsort(vals)[::-1]
    axis = vecs[:, order[0]]
    if axis[2] < 0:
        axis = -axis
    return axis, vals[order]


def fit_plane(points):
    """Total-least-squares plane fit by SVD (manuscript Eqs. 6 and 7).

    Returns a dict with the unit normal, centroid, offset ``d``, the RMS
    perpendicular residual and the dip / dip direction.
    """
    P = np.asarray(points, dtype=float)
    if len(P) < 3:
        return dict(normal=np.array([0.0, 0.0, 1.0]), centroid=P.mean(axis=0) if len(P) else np.zeros(3),
                    d=0.0, rmse=0.0, dip=0.0, dipdir=0.0, n_points=len(P))
    c = P.mean(axis=0)
    A = P - c
    _, _, Vt = np.linalg.svd(A, full_matrices=False)
    n = Vt[-1]
    if n[2] < 0:
        n = -n
    n = n / np.maximum(np.linalg.norm(n), 1e-300)
    dist = A @ n
    dip, dipdir = dip_from_normal(n)
    return dict(normal=n, centroid=c, d=float(-n @ c), rmse=float(np.sqrt(np.mean(dist ** 2))),
                dip=dip, dipdir=dipdir, n_points=len(P))


def projected_hull(points, normal):
    """Exposed area and persistence of a patch (manuscript Eq. 10).

    The points are projected onto their best-fit plane, the 2-D convex hull is
    built, the enclosed polygon area is the exposed area and the maximum
    pairwise distance between hull vertices is the persistence (maximum
    exposed trace length).
    """
    from scipy.spatial import ConvexHull, QhullError

    P = np.asarray(points, dtype=float)
    if len(P) < 3:
        return 0.0, 0.0
    n = np.asarray(normal, dtype=float)
    n = n / np.maximum(np.linalg.norm(n), 1e-300)
    # orthonormal basis of the plane
    seed = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(n, seed)
    e1 /= np.maximum(np.linalg.norm(e1), 1e-300)
    e2 = np.cross(n, e1)
    XY = (P - P.mean(axis=0)) @ np.stack([e1, e2], axis=1)
    try:
        hull = ConvexHull(XY)
    except (QhullError, ValueError):
        return 0.0, 0.0
    verts = XY[hull.vertices]
    area = float(hull.volume)                      # 2-D ConvexHull.volume == area
    diff = verts[:, None, :] - verts[None, :, :]
    persistence = float(np.sqrt((diff ** 2).sum(-1)).max()) if len(verts) > 1 else 0.0
    return area, persistence


def line_of_intersection(n1, n2):
    """Downward-oriented line of intersection of two planes (Eq. 11).

    Returns ``(plunge_deg, trend_deg)``.  The trend is the azimuth of the
    horizontal projection of the downward-pointing intersection vector.
    """
    n1 = np.atleast_2d(np.asarray(n1, dtype=float))
    n2 = np.atleast_2d(np.asarray(n2, dtype=float))
    L = np.cross(n1, n2)
    norm = np.linalg.norm(L, axis=1)
    ok = norm > 1e-12
    out_plunge = np.full(len(L), np.nan)
    out_trend = np.full(len(L), np.nan)
    Lv = L[ok] / norm[ok, None]
    Lv = np.where((Lv[:, 2] > 0)[:, None], -Lv, Lv)      # force downward
    horiz = np.hypot(Lv[:, 0], Lv[:, 1])
    out_plunge[ok] = np.degrees(np.arctan2(np.abs(Lv[:, 2]), horiz))
    out_trend[ok] = np.degrees(np.arctan2(Lv[:, 0], Lv[:, 1])) % 360.0
    if len(L) == 1:
        return float(out_plunge[0]), float(out_trend[0])
    return out_plunge, out_trend
