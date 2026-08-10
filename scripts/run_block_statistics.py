#!/usr/bin/env python3
"""Reproduce manuscript Section 4.4 and Figures 16 to 18 from a block inventory.

    python scripts/run_block_statistics.py \
        --inventory data/chilas_block_inventory.csv \
        --output-dir output/block_statistics
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import _bootstrap  # noqa: F401

from georank3d.blockstats import analyse_blocks, read_block_inventory
from georank3d.pipeline import write_csv


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inventory", required=True, type=Path,
                    help="CSV with a volume column, optionally block_id and sector")
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--bootstrap", type=int, default=10_000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)

    try:
        records = read_block_inventory(args.inventory)
    except (FileNotFoundError, ValueError) as exc:
        print(f"input error: {exc}")
        return 2

    args.output_dir.mkdir(parents=True, exist_ok=True)
    analysis = analyse_blocks(records, n_bootstrap=args.bootstrap, seed=args.seed)

    d = analysis["descriptive"]
    ln = analysis["distributions"]["lognormal"]
    print(f"blocks               : {d['count']}")
    print(f"total volume         : {d['sum_m3']:.2f} m3")
    print(f"mean / median        : {d['mean_m3']:.2f} / {d['median_m3']:.2f} m3")
    print(f"range                : {d['min_m3']:.2f} - {d['max_m3']:.2f} m3")
    print(f"coefficient of var.  : {d['coefficient_of_variation']:.2f}")
    print(f"ln-mean / ln-sd      : {ln['ln_mean']:.2f} / {ln['ln_std']:.2f}"
          f"  (geometric mean {ln['geometric_mean_m3']:.2f} m3)")
    print(f"KS on ln V           : D = {ln['ks_statistic']:.3f}, p = {ln['ks_pvalue']:.3f}")
    print(f"Shapiro-Wilk on ln V : W = {ln['shapiro_w']:.3f}, p = {ln['shapiro_pvalue']:.3f}")
    print(f"normal Q-Q R^2       : {ln['qq_r_squared']:.3f}")
    print(f"Gini coefficient     : {analysis['gini']:.2f}")
    print("\nmodel ranking by AIC")
    for m in analysis["distributions"]["models"]:
        print(f"  {m['model']:<12} logL = {m['loglik']:9.3f}   AIC = {m['aic']:8.2f}   BIC = {m['bic']:8.2f}")
    print("\nvolume concentration")
    for t in analysis["magnitude_frequency"]["thresholds"]:
        print(f"  blocks >= {t['threshold_m3']:>5.1f} m3 : n = {t['n_ge']:>3}, "
              f"{t['volume_share_percent']:.2f}% of mapped volume")
    print("\nsource sectors")
    for s in analysis["sectors"]["sectors"]:
        lo, hi = s["share_ci95_percent"]
        print(f"  {s['sector']:<10} n = {s['n_blocks']:>2}  V = {s['total_volume_m3']:>7.2f} m3  "
              f"{s['volume_share_percent']:>6.2f}%  (bootstrap 95% CI {lo:.2f}-{hi:.2f})  "
              f"median {s['median_volume_m3']:.2f} m3  blocks >= 20 m3: {s['n_blocks_ge_20m3']}")
    print(f"  Kruskal-Wallis on median block size: H = {analysis['sectors']['kruskal_h']:.3f}, "
          f"p = {analysis['sectors']['kruskal_pvalue']:.3f}")

    write_csv(args.output_dir / "block_inventory_classified.csv", analysis["inventory"])
    write_csv(args.output_dir / "distribution_model_comparison.csv",
              analysis["distributions"]["models"])
    write_csv(args.output_dir / "sector_summary.csv", analysis["sectors"]["sectors"])
    write_csv(args.output_dir / "volume_exceedance.csv",
              analysis["magnitude_frequency"]["exceedance"])
    (args.output_dir / "block_statistics.json").write_text(json.dumps(analysis, indent=2, default=str))

    if not args.no_figures:
        from georank3d.plotting import block_statistics_figures
        for p in block_statistics_figures(args.output_dir, analysis):
            print(f"figure: {p}")
    print(f"\noutputs written to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
