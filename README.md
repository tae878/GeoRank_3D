# GeoRank_3D — corrected reference implementation and reproducibility audit

This repository contains an independent audit of the manuscript *Discontinuity-
Controlled Rockfall Source Identification and Trajectory Analysis for Hazard
Zonation at Chilas Slope, Karakoram Highway* and its accompanying MATLAB
archive, together with a corrected, executable implementation of the workflow
the manuscript describes.

**Start with [`AUDIT_REPORT.md`](AUDIT_REPORT.md)** — findings, evidence and the
list of corrections required before publication.

## What is here

| Path | Contents |
|---|---|
| `AUDIT_REPORT.md` | The audit: inventory, claim-by-claim correspondence, defects, verdict |
| `georank3d/` | Corrected Python implementation of manuscript Section 3.3 and Section 4.4 |
| `scripts/` | Two command-line entry points and a synthetic-slope generator |
| `configs/` | JSON run configurations (Chilas parameters, synthetic demo) |
| `data/` | Table 4 and Table 5 transcribed to CSV — see `data/README.md` |
| `matlab_patches/` | Drop-in corrections for the authors' own MATLAB release |
| `tests/` | 62 tests, including regression guards for every defect found |

## Install

Python 3.10 or newer.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt        # numpy, scipy, matplotlib
python -m pytest                       # 62 tests, ~6 s
```

No compiled extensions, no clustering library, no GPU. The package reads
ASCII (`.xyz`, `.txt`, `.csv`, `.asc`, `.pts`) and PLY (ASCII or binary) clouds.

## Run

### Discontinuity pipeline

```bash
# build a synthetic slope with four known orientations
python scripts/make_synthetic_slope.py --output output/synthetic/synthetic_slope.xyz

# run the full workflow
python scripts/run_pipeline.py --config configs/synthetic_demo.json
```

On a real cloud, point `configs/chilas.json` at it:

```bash
python scripts/run_pipeline.py --config configs/chilas.json \
    --set input_cloud=/path/to/chilas.xyz \
    --set output_dir=output/chilas
```

Any configuration entry can be overridden with `--set key=value`, so a
sensitivity analysis needs no edits to source or config:

```bash
for w in '[0.2,0.25,0.15,0.2,0.2]' '[0.1,0.3,0.2,0.2,0.2]' '[0.0,0.3,0.2,0.3,0.2]'; do
  python scripts/run_pipeline.py --config configs/chilas.json \
      --set input_cloud=/path/to/chilas.xyz \
      --set "reliability_weights=$w" --set output_dir="output/weights_$w"
done
```

`--print-config` shows the resolved configuration without running.

Outputs (all written to `output_dir`): `set_metrics.csv`, `patch_metrics.csv`,
`wedge_pairs_mean.csv`, `kinematic_dispersion_sets.csv`,
`kinematic_dispersion_wedges.csv`, `uncertainty_runs.csv`,
`uncertainty_sets.csv`, `validation_*.csv`, `classified_cloud.xyz`,
`pole_plot_schmidt.png`, `reliability_ranking.png`, `results.json`,
`resolved_config.json` and `run_provenance.json`.

### Block-volume statistics (manuscript Section 4.4, Figs. 16–18)

```bash
python scripts/run_block_statistics.py \
    --inventory data/chilas_block_inventory.csv \
    --output-dir output/block_statistics
```

Reproduces every published value: total 222.13 m³, CV 1.49, ln-mean 0.94,
ln-sd 1.53, KS D = 0.120 p = 0.722, Shapiro–Wilk p = 0.637, Q-Q R² = 0.976,
log-normal AIC 176.66 / BIC 179.53, Gini 0.66, sector shares
74.59 / 19.04 / 6.37 %, Kruskal–Wallis p = 0.424. These are asserted in
`tests/test_blockstats.py` against the manuscript's numbers.

## How this implementation differs from the MATLAB archive

Each difference is a correction of a defect documented in `AUDIT_REPORT.md`,
not a change of method.

| Change | Reason |
|---|---|
| Dip direction is `atan2(n_x, n_y)` of the upward normal | Manuscript Eq. (9) and `GeoRank_3D_orientation.m` return it 180° out (§5.1) |
| One orientation convention, used everywhere | The archive contains two, with opposite signs (§5.2) |
| Axial cone membership `arccos(\|n·p\|)`; set means from the orientation tensor | The chord test on hemisphere-flipped normals truncates steeply dipping sets at the `n_z = 0` seam (§5.3) |
| GRI weights read from the configuration and validated | `opts.reliabilityWeights` was declared but never read (§5.4) |
| GRI persistence fed from maximum patch persistence, and all five terms exported | The set-wide convex hull saturates the term for every set (§5.5) |
| Monte Carlo reports percentile intervals and the Fisher α₉₅ cone | A standard error of the mean across runs shrinks with the run count (§5.8) |
| Two-way validation: reference → extracted as well as extracted → reference | One-way matching hides reference sets that no extracted set recovers (§5.1) |
| Block-volume ceiling sensitivity and KS p-value reported | The ceiling is derived from the census being tested (§5.17) |
| Validation and block comparison available from the CLI | They were reachable only from the GUI (§5.15) |
| Section 4.4 statistics implemented | No code for them existed (§8) |
| Provenance, input validation, informative errors, tests | None existed (§6) |

The Schmidt pole plot is unchanged in behaviour: `GeoRank_3D_plot_poles.m`
contains a second 180° error that cancels the first, so the published stereonet
is correct. Fixing only the orientation convention would silently rotate it —
see `matlab_patches/GeoRank_3D_plot_poles_patch_notes.m`.

## Scope

Implemented: manuscript Section 3.3 (Eqs. 1–12) and Section 4.4.

Not implemented, deliberately: RocFall3 trajectory simulation (commercial
third-party software; Section 4.5), the CloudCompare manual cleaning and 2.5D
volume computation (interactive, and no reproducible parameterisation is given),
and DSE / qFacet extraction (external published tools). Nothing has been added
that the manuscript does not describe.

## Coordinate and orientation conventions

X = east, Y = north, Z = up. A plane dipping `δ` towards `α` has upward unit
normal `(sin δ sin α, sin δ cos α, cos δ)`; its horizontal projection points
**down dip**. The lower-hemisphere pole plunges `90 − δ` towards `α + 180` and
plots at radius `√2 sin(δ/2)`. Normals are axial: `n` and `−n` denote the same
plane, angles use `arccos(|n₁·n₂|)`, and mean orientations come from the
orientation-tensor principal eigenvector.

If your cloud is not in an east/north/up frame, transform it before import — the
package has no way to detect this, and an incorrect frame silently rotates every
reported azimuth.
