#!/usr/bin/env python3
"""Generate a synthetic rock slope with known discontinuity sets.

The cloud contains a slope face plus three joint sets whose orientations are
supplied on the command line (defaults: the Chilas manual-mapping attitudes
JS-1 82/011, JS-2 67/029 and JS-3 78/346 of manuscript Table 4).  Because the
ground truth is known exactly, this is the end-to-end test that the pipeline
recovers the orientations it is given.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import _bootstrap  # noqa: F401  (adds the repository root to sys.path)
from georank3d.geometry import normal_from_dip


def plane_patch(dip, dipdir, centre, extent, n, roughness, rng):
    """Points on a rough planar patch of a given attitude, with waviness."""
    n_vec = normal_from_dip(dip, dipdir)
    seed = np.array([1.0, 0.0, 0.0]) if abs(n_vec[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(n_vec, seed)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n_vec, e1)
    uv = (rng.random((n, 2)) - 0.5) * np.asarray(extent)
    waviness = roughness * np.sin(uv[:, 0] / max(extent[0] / 6.0, 1e-6)) \
        + roughness * rng.standard_normal(n) * 0.3
    return np.asarray(centre) + uv[:, :1] * e1 + uv[:, 1:2] * e2 + waviness[:, None] * n_vec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output", required=True, type=Path, help="output .xyz path")
    ap.add_argument("--sets", nargs="*", default=["82/11", "67/29", "78/346"],
                    help="joint-set attitudes as dip/dipdir (default: Chilas Table 4 manual sets)")
    ap.add_argument("--slope", default="80/135", help="slope-face attitude dip/dipdir")
    ap.add_argument("--points-per-set", type=int, default=9000)
    ap.add_argument("--slope-points", type=int, default=14000)
    ap.add_argument("--noise-points", type=int, default=1500)
    ap.add_argument("--roughness", type=float, default=0.02, help="surface roughness amplitude (m)")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    clouds = []

    s_dip, s_dir = (float(x) for x in args.slope.split("/"))
    clouds.append(plane_patch(s_dip, s_dir, [0, 0, 18], (60, 34), args.slope_points,
                              args.roughness * 2.5, rng))

    centres = [[-14, 4, 26], [12, -6, 14], [2, 8, 6], [-6, -10, 22], [16, 6, 24]]
    for i, spec in enumerate(args.sets):
        dip, dipdir = (float(x) for x in spec.split("/"))
        per_patch = max(args.points_per_set // 3, 200)
        for j in range(3):
            centre = centres[(i * 3 + j) % len(centres)]
            centre = [centre[0] + 3 * j, centre[1] - 2 * i, centre[2] - 2 * j]
            clouds.append(plane_patch(dip, dipdir, centre, (11, 8), per_patch,
                                      args.roughness, rng))

    clouds.append((rng.random((args.noise_points, 3)) - 0.5) * np.array([60.0, 40.0, 34.0])
                  + np.array([0.0, 0.0, 17.0]))

    P = np.vstack(clouds)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output, P, fmt="%.4f")
    print(f"wrote {len(P)} points to {args.output}")
    print("ground-truth sets: " + ", ".join(args.sets) + f"   slope face: {args.slope}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
