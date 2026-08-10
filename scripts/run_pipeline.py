#!/usr/bin/env python3
"""Command-line driver for the GeoRank_3D discontinuity pipeline.

    python scripts/run_pipeline.py --config configs/chilas.json
    python scripts/run_pipeline.py --config configs/chilas.json \
        --set input_cloud=/data/chilas.xyz --set slope_dipdir_deg=140

Every parameter can be overridden from the command line with ``--set key=value``
so that sensitivity analyses need no edits to source or configuration files.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import _bootstrap  # noqa: F401

from georank3d.config import Config, ConfigError
from georank3d.pipeline import run_pipeline


def _coerce(raw: str):
    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        return raw


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, help="JSON configuration file")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="override a configuration entry (repeatable)")
    ap.add_argument("--print-config", action="store_true",
                    help="print the resolved configuration and exit without running")
    args = ap.parse_args(argv)

    try:
        data = json.loads(args.config.read_text()) if args.config else {}
        for item in args.set:
            if "=" not in item:
                raise ConfigError(f"--set expects KEY=VALUE, got {item!r}")
            key, value = item.split("=", 1)
            data[key.strip()] = _coerce(value.strip())
        cfg = Config.from_dict(data)
    except (ConfigError, FileNotFoundError, json.JSONDecodeError) as exc:
        print(f"configuration error: {exc}", file=sys.stderr)
        return 2

    if args.print_config:
        print(json.dumps(cfg.to_dict(), indent=2, sort_keys=True))
        return 0
    if not cfg.input_cloud:
        print("configuration error: 'input_cloud' is empty; point it at an "
              "XYZ/TXT/CSV/PLY point cloud", file=sys.stderr)
        return 2

    try:
        results = run_pipeline(cfg)
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"run failed: {exc}", file=sys.stderr)
        return 1

    print("\n=== discontinuity sets, ranked by Geomechanical Reliability Index ===")
    header = f"{'set':>4} {'N':>9} {'dip':>6} {'dipdir':>7} {'GRI':>6} {'accept':>7} " \
             f"{'planar':>7} {'topple':>7}"
    print(header)
    print("-" * len(header))
    for r in results["sets"]:
        print(f"{r['set_id']:>4} {r['n_points']:>9} {r['dip_deg']:>6.1f} {r['dipdir_deg']:>7.1f} "
              f"{r['gri']:>6.3f} {str(r['accepted']):>7} "
              f"{r['planar_feasible_mean']:>7} {r['toppling_feasible_mean']:>7}")

    if results["set_dispersion"]:
        print("\n=== dispersion-aware kinematic screening (population vs set mean) ===")
        for r in results["set_dispersion"]:
            print(f"  set {r['set_id']}: planar {100 * r['planar_feasible_fraction']:5.1f}% / "
                  f"toppling {100 * r['toppling_feasible_fraction']:5.1f}% of population feasible "
                  f"(mean flags {r['planar_feasible_mean']}/{r['toppling_feasible_mean']})")
    for r in results["wedge_dispersion"]:
        print(f"  wedge {r['set_a']}-{r['set_b']}: "
              f"{100 * r['wedge_feasible_fraction']:5.1f}% of sampled intersections feasible "
              f"(mean flag {r['wedge_feasible_mean']}, P95 plunge {r['p95_plunge_deg']:.1f} deg)")

    print(f"\noutputs written to {cfg.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
