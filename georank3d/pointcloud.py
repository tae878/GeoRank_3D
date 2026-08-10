"""Point-cloud import, voxel reduction and PCA normal estimation.

Implements manuscript Section 3.3.2 stages 1-3: import and preprocessing,
local normal and feature estimation (Eqs. 1-3), and planarity / roughness
filtering.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
from scipy.spatial import cKDTree

__all__ = ["read_point_cloud", "voxel_downsample", "reduce_to_budget",
           "estimate_normals", "PointFeatures"]

_ASCII_SUFFIXES = {".xyz", ".txt", ".csv", ".asc", ".pts"}


def read_point_cloud(path) -> np.ndarray:
    """Read an ASCII (XYZ/TXT/CSV/ASC/PTS) or binary/ASCII PLY point cloud.

    Only the first three numeric columns are used.  LAS/LAZ must be converted
    beforehand, as stated in manuscript Section 3.3.2.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(
            f"point cloud not found: {p}\n"
            "Set 'input_cloud' in the configuration to an existing "
            "XYZ/TXT/CSV/PLY file (LAS/LAZ must be converted first)."
        )
    suffix = p.suffix.lower()
    if suffix == ".ply":
        P = _read_ply(p)
    elif suffix in _ASCII_SUFFIXES:
        P = _read_ascii(p)
    else:
        raise ValueError(f"unsupported point-cloud extension {suffix!r}; "
                         f"expected one of {sorted(_ASCII_SUFFIXES | {'.ply'})}")
    if P.ndim != 2 or P.shape[1] < 3:
        raise ValueError(f"{p} must contain at least three columns (X Y Z); got shape {P.shape}")
    P = np.asarray(P[:, :3], dtype=float)
    finite = np.isfinite(P).all(axis=1)
    if not finite.all():
        P = P[finite]
    if len(P) < 10:
        raise ValueError(f"{p} contains only {len(P)} finite points; nothing to analyse")
    return P


def _read_ascii(p: Path) -> np.ndarray:
    for delim in (None, ",", ";", "\t"):
        try:
            arr = np.loadtxt(p, delimiter=delim, comments=("#", "//"), ndmin=2)
        except ValueError:
            continue
        if arr.size and arr.shape[1] >= 3:
            return arr
    # last resort: skip a header line
    for delim in (None, ",", ";", "\t"):
        try:
            arr = np.loadtxt(p, delimiter=delim, comments=("#", "//"), skiprows=1, ndmin=2)
        except ValueError:
            continue
        if arr.size and arr.shape[1] >= 3:
            return arr
    raise ValueError(f"could not parse {p} as an ASCII point cloud with >= 3 numeric columns")


def _read_ply(p: Path) -> np.ndarray:
    _NP = {"float": "f4", "float32": "f4", "double": "f8", "float64": "f8",
           "char": "i1", "int8": "i1", "uchar": "u1", "uint8": "u1",
           "short": "i2", "int16": "i2", "ushort": "u2", "uint16": "u2",
           "int": "i4", "int32": "i4", "uint": "u4", "uint32": "u4"}
    with open(p, "rb") as fh:
        if fh.readline().strip() != b"ply":
            raise ValueError(f"{p} is not a PLY file")
        fmt, props, count, in_vertex = None, [], 0, False
        while True:
            line = fh.readline()
            if not line:
                raise ValueError(f"{p}: unexpected end of PLY header")
            tok = line.decode("ascii", "replace").split()
            if not tok:
                continue
            if tok[0] == "format":
                fmt = tok[1]
            elif tok[0] == "element":
                in_vertex = tok[1] == "vertex"
                if in_vertex:
                    count = int(tok[2])
            elif tok[0] == "property" and in_vertex:
                if tok[1] == "list":
                    raise ValueError(f"{p}: list properties on vertices are not supported")
                props.append((tok[2], _NP[tok[1]]))
            elif tok[0] == "end_header":
                break
        names = [n for n, _ in props]
        for axis in ("x", "y", "z"):
            if axis not in names:
                raise ValueError(f"{p}: PLY vertex element lacks property {axis!r}")
        if fmt == "ascii":
            arr = np.loadtxt(fh, max_rows=count, ndmin=2)
            cols = [names.index(a) for a in ("x", "y", "z")]
            return arr[:, cols]
        dtype = np.dtype([(n, ("<" if "little" in fmt else ">") + t) for n, t in props])
        raw = np.frombuffer(fh.read(count * dtype.itemsize), dtype=dtype, count=count)
        return np.stack([raw["x"], raw["y"], raw["z"]], axis=1).astype(float)


def voxel_downsample(P: np.ndarray, voxel_size: float) -> np.ndarray:
    """Centroid-based voxel-grid down-sampling."""
    if voxel_size <= 0 or len(P) == 0:
        return P
    origin = P.min(axis=0)
    keys = np.floor((P - origin) / voxel_size).astype(np.int64)
    _, inverse, counts = np.unique(keys, axis=0, return_inverse=True, return_counts=True)
    out = np.zeros((len(counts), 3))
    for c in range(3):
        out[:, c] = np.bincount(inverse, weights=P[:, c], minlength=len(counts)) / counts
    return out


def reduce_to_budget(P: np.ndarray, cfg) -> Tuple[np.ndarray, Dict]:
    """Voxel filter, then coarsen uniformly until the analysis budget is met.

    Uniform coarsening (rather than random thinning) preserves the even
    spatial coverage on which the orientation statistics depend.
    """
    info = {"input_count": int(len(P)), "requested_voxel_size_m": float(cfg.voxel_size_m),
            "budget_method": "none", "applied_voxel_size_m": float(cfg.voxel_size_m)}
    out = voxel_downsample(P, cfg.voxel_size_m) if cfg.voxel_size_m > 0 else P
    info["after_voxel_count"] = int(len(out))

    if len(out) > cfg.max_points:
        if cfg.uniform_reduce_to_budget and cfg.voxel_size_m > 0:
            lo = cfg.voxel_size_m
            hi = max(cfg.voxel_size_m * np.sqrt(len(out) / cfg.max_points), cfg.voxel_size_m * 1.01)
            cand = voxel_downsample(out, hi)
            grow = 0
            while len(cand) > cfg.max_points and grow < 30:
                lo, hi = hi, hi * 1.4
                cand = voxel_downsample(out, hi)
                grow += 1
            for _ in range(18):
                mid = 0.5 * (lo + hi)
                trial = voxel_downsample(out, mid)
                if len(trial) > cfg.max_points:
                    lo = mid
                else:
                    hi, cand = mid, trial
            out = cand
            info["budget_method"] = "uniform_voxel"
            info["applied_voxel_size_m"] = float(hi)
        if len(out) > cfg.max_points:                      # safety net
            rng = np.random.default_rng(cfg.random_seed)
            out = out[rng.choice(len(out), cfg.max_points, replace=False)]
            info["budget_method"] = "random_subsample"
    info["working_count"] = int(len(out))
    return out, info


@dataclass
class PointFeatures:
    """Per-point PCA descriptors (manuscript Eqs. 2 and 3)."""
    normals: np.ndarray        # upward-oriented unit normals, (n, 3)
    planarity: np.ndarray      # P_i = (l2 - l1) / l3
    linearity: np.ndarray      # diagnostic only, per Section 3.3.1
    surface_variation: np.ndarray   # S_i = l1 / (l1 + l2 + l3)
    roughness: np.ndarray      # sqrt(l1): RMS distance to the local plane

    def subset(self, mask) -> "PointFeatures":
        return PointFeatures(self.normals[mask], self.planarity[mask], self.linearity[mask],
                             self.surface_variation[mask], self.roughness[mask])


def estimate_normals(P: np.ndarray, cfg, block_size: int = 200_000) -> PointFeatures:
    """Blocked PCA normal estimation with eigenvalue descriptors.

    For every point the covariance of its ``normal_k`` nearest neighbours is
    formed (Eq. 2); the eigenvector of the smallest eigenvalue is the surface
    normal and the eigenvalues give planarity, linearity and surface variation
    (Eq. 3).  Local roughness is ``sqrt(lambda_min)``, the RMS distance of the
    neighbourhood to its own best-fit plane.

    Queries run in blocks so that peak memory is ``block_size * k * 3``
    doubles rather than ``n * k * 3``.
    """
    n = len(P)
    k = int(cfg.normal_k)
    if n <= k:
        raise ValueError(f"the working cloud has {n} points but normal_k = {k}; "
                         "reduce normal_k or supply a denser cloud")
    tree = cKDTree(P)
    normals = np.empty((n, 3))
    planarity = np.empty(n)
    linearity = np.empty(n)
    surfvar = np.empty(n)
    roughness = np.empty(n)

    for start in range(0, n, block_size):
        stop = min(start + block_size, n)
        _, idx = tree.query(P[start:stop], k=k + 1, workers=-1)
        nbr = P[idx[:, 1:]]                              # (b, k, 3), self dropped
        centred = nbr - nbr.mean(axis=1, keepdims=True)
        cov = np.einsum("bki,bkj->bij", centred, centred) / k
        vals, vecs = np.linalg.eigh(cov)                 # ascending eigenvalues
        l1, l2, l3 = vals[:, 0], vals[:, 1], vals[:, 2]
        nvec = vecs[:, :, 0]
        nvec = np.where((nvec[:, 2] < 0)[:, None], -nvec, nvec)
        normals[start:stop] = nvec / np.maximum(np.linalg.norm(nvec, axis=1, keepdims=True), 1e-300)
        denom3 = np.maximum(l3, 1e-300)
        total = np.maximum(l1 + l2 + l3, 1e-300)
        planarity[start:stop] = (l2 - l1) / denom3
        linearity[start:stop] = (l3 - l2) / denom3
        surfvar[start:stop] = l1 / total
        roughness[start:stop] = np.sqrt(np.maximum(l1, 0.0))

    return PointFeatures(normals, planarity, linearity, surfvar, roughness)
