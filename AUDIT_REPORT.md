# Reproducibility and paper–code consistency audit

**Manuscript.** *Discontinuity-Controlled Rockfall Source Identification and Trajectory
Analysis for Hazard Zonation at Chilas Slope, Karakoram Highway, NW Himalayas, Pakistan*
(Qaisar, Zhao, Khan and Qaisar), 62 pp., submitted to JRMGE.

**Code archive.** `03_GeoRank_3D_Code_and_GUI.zip` — 41 files, 40 of them MATLAB
sources, ~120 kB, no data, no tests, no results, no environment specification.

**Audit date.** 10 August 2026. All numerical results below were recomputed
independently from the manuscript's own tables and from a line-by-line
re-implementation of the algorithms; nothing is quoted on trust.

---

## 1. Executive summary

The block-volume statistics in Section 4.4 are excellent: every reported value
in Figs. 16–18 reproduces exactly (Section 8 below). The engineering narrative
is coherent and the study addresses a real gap.

The GeoRank_3D component, which is the paper's sole claimed novelty, does not
survive audit in its present form. Three findings are decisive:

1. **The dip-direction convention is 180° in error.** Manuscript Eq. (9),
   `alpha = mod(atan2(-a, -b), 360)`, is applied to a normal that Section 3.3.3
   explicitly orients "into the upper hemisphere". The minus signs belong with a
   *downward*-oriented normal. `GeoRank_3D_orientation.m:9` and
   `GeoRank_3D_dipdir2normal.m:5` implement Eq. (9) faithfully, so every
   GeoRank_3D attitude in Table 8, Fig. 22 and Fig. 23 is rotated by 180°: Set 3
   is 80/197, not 80/017; Set 2 is 73/005; Set 1 is 72/028.
   In `GeoRank_3D_validate_multi_source.m` the extracted set enters as its own
   PCA normal while each reference attitude is passed through the flipped
   formula, so **only one of the two vectors is rotated** and the Table 8
   deviations are not the angles between the planes. Recomputed consistently,
   the match structure inverts: Set 3 vs manual JS‑1 becomes 19.0° instead of
   6.5°, and the closest field match to any GeoRank set is **Set 1 — the
   population the paper dismisses as a non-structural exposure surface — at 5.1°
   from manual JS‑2**. The abstract's headline claim that GeoRank_3D
   "independently reproduced the dominant joint set (6.5° from field mapping)"
   does not survive this correction.

2. **The Geomechanical Reliability Index does no discriminating work.** With the
   abundance term saturated above 10⁵ points and the persistence term fed from
   a set-wide convex hull that always exceeds the 10 m reference, Eq. (5)
   collapses to `0.40 + 0.25 P + 0.15 R + 0.20 K` with `K = 0` for all three sets
   (Table 9). That is exactly why Table 8 reports 0.730, 0.730 and 0.728 — a
   spread of **0.002** across a field-confirmed joint set, a DSE-corroborated
   set and an exposure surface. The GRI ≥ 0.50 threshold is never binding, so
   "ranking sets by reliability rather than statistical abundance" is not what
   the code does at this site.

3. **Most of the paper is not in the archive.** Sections 4.1–4.3 (manual
   mapping, DSE, both qFacet variants, all Dips kinematics, Figs. 9–14), Section
   4.4 (Table 5, Figs. 16–18, the log-normal fitting, Gini, Kruskal–Wallis,
   bootstrap CIs) and Section 4.5 (all RocFall3 results, Figs. 19–21) have **no
   corresponding code and no supplied data**. The archive covers Section 3.3
   only. Even within Section 3.3, Table 8 and Table 10 can be produced only
   through the GUI; no command-line driver reaches the validation layer.

**Verdict: Insufficient/incomplete for validating the paper's claims**, with the
GeoRank_3D component **substantially inconsistent** with what the manuscript
reports. This is fixable — the defects are localised and the underlying design
is sound — but the corrections change reported numbers, so a revision cannot be
editorial.

Deliverables in this repository: a corrected, executable reference
implementation (`georank3d/`, 62 passing tests), a reproduction of Section 4.4
that matches every published value, and drop-in MATLAB patches
(`matlab_patches/`) so the authors can keep their own package and GUI.

---

## 2. Inventory of the ZIP archive

```
GeoRank_3D/
├── GeoRank_3D.m                     main pipeline (import → … → export)
├── GeoRank_3D_KKH.m                 near-duplicate pipeline, missing 2 stages
├── GeoRank_3D_App.m                 App Designer GUI, 8 tabs (26 kB)
├── GeoRank_3D_demo_synthetic.m      synthetic 3-plane demo
├── run_GeoRank_3D.m                 CLI driver → GeoRank_3D.m
├── run_GeoRank_3D_KKH.m             CLI driver → GeoRank_3D_KKH.m
├── README.md                        install / performance notes
├── FIX_NOTES_clustering.md          changelog: DBSCAN → density-peak
├── FIX_NOTES_dispersion.md          changelog: dispersion screening added
├── docs/                            3 markdown notes
└── functions/                       30 GeoRank_3D_*.m modules
```

### Call graph and per-file purpose

| File | Purpose | Reached from |
|---|---|---|
| `GeoRank_3D.m` | Orchestrates stages 1–9 + export | `run_GeoRank_3D.m`, GUI |
| `GeoRank_3D_KKH.m` | Same, **omits** dispersion screening and MC cluster tracking | `run_GeoRank_3D_KKH.m`, demo |
| `_read_pointcloud` | XYZ/TXT/CSV via `readmatrix`, PLY via `pcread` | both pipelines |
| `_prepare_points` / `_voxel_downsample` | Voxel filter + bracketed uniform reduce-to-budget | both |
| `_check_memory` | Advisory RAM estimate (Windows `memory` only) | `GeoRank_3D.m` only |
| `_estimate_normals` / `_sym3x3_eig` | Blocked PCA normals, analytic 3×3 eigensolver (Smith 1961) | both |
| `_filter_struct` | Masks the per-point feature struct | both |
| `_cluster_normals` | Density-peak (mean-shift) pole finding + cone assignment; legacy DBSCAN/k-means fallback | both, and MC |
| `_unify_normals` | Forces `nz ≥ 0` | clustering, orientation |
| `_spatial_patches` | Spatial DBSCAN per set, subsample-and-propagate above cap | set metrics |
| `_compute_set_metrics` | Per-set dip/dipdir, planarity, roughness, area, persistence, GRI | both |
| `_fit_plane` | SVD plane fit, Eqs. (6)–(7) | set/patch metrics |
| `_projected_area` | Convex-hull area + max chord, Eq. (10) | set/patch metrics |
| `_orientation` / `_dipdir2normal` | Eqs. (8)–(9) and inverse | everywhere |
| `_kinematic_index` | Markland planar + direct toppling on the set mean | set metrics |
| `_wedge_analysis` | Pairwise Eq. (11) on set means | `GeoRank_3D.m`, `_KKH.m` |
| `_kinematic_dispersion` | Population planar/toppling/wedge fractions | **`GeoRank_3D.m` only** |
| `_reliability_score` | Eq. (5) | set metrics |
| `_uncertainty` / `_perturb_normals` | 50-run Monte Carlo, cluster tracking | both |
| `_export_results` | CSV + classified XYZ + run report | both |
| `_plot_poles` | Schmidt net + reliability bar chart | both |
| `_validate_multi_source` | Deviation vs manual/qFacet/DSE, 10° threshold | **GUI only** |
| `_compare_block_volumes` | V = 0.20 A^1.5 vs field census, KS test | **GUI only** |
| `_read_reference_csv` | Reference CSV parser | **GUI only** |
| `_orientation_deviation` | **Dead code** — nothing calls it | — |

### What is absent

No point cloud, no reference CSVs, no block inventory, no saved results, no
figures, no tests, no CI, no `LICENSE`, no version pin, no `requirements`
equivalent, no RocFall3 project or export, no CloudCompare 2.5D volume project,
no DSE or qFacet exports, no Dips files. The manuscript's Data Availability
statement ("available from the corresponding author upon reasonable request")
is the only pointer.

### Executability

Neither pipeline was executable in the audit environment. The archive requires
MATLAB with the Statistics and Machine Learning Toolbox (`knnsearch`,
`createns`, `KDTreeSearcher`, `rangesearch`, `dbscan`, `kmeans`, `silhouette`,
`pdist`, `prctile`), the Computer Vision Toolbox (`pcread`, for PLY input), the
`table` class (`table`, `cell2table`, `writetable`, `readtable`, `height` — 59
call sites) and App Designer (`matlab.apps.AppBase`). No MATLAB licence was
available; GNU Octave 8.4.0 with `statistics` 1.6.3 was used as the nearest
substitute and the demo halts immediately:

```
EXECUTION FAILED: 'writematrix' undefined near line 32, column 1
  at GeoRank_3D_demo_synthetic line 32
```

`table`, `cell2table`, `writetable`, `dbscan`, `createns`, `KDTreeSearcher` and
`pcread` are also unavailable outside MATLAB, so the failure is structural, not
incidental. **This is not a defect** — MATLAB-only code is normal in this
field — but it does mean no result in the manuscript could be checked by
execution, and the archive contains nothing (no test, no expected output, no
saved run) that would let a reviewer with MATLAB confirm that the shipped code
produced the shipped numbers. Both hard-coded input paths
(`'your_pointcloud.xyz'`, `'your_point_cloud.xyz'`) are placeholders.

---

## 3. Summary of the paper's methodology and claims

| # | Claim | Location |
|---|---|---|
| C1 | 65 M-point VOM from 324 UAV images, GSD 2.5 cm, 8 RTK GCPs | §3.1, p. 13–14 |
| C2 | Manual mapping resolved 4,818 discontinuities → 3 sets | §4.2, p. 27 |
| C3 | DSE returned 1,472 discontinuities → 10 clusters | §4.3.1, Table 1 |
| C4 | qFacet FM 1,872 planes → 3 sets; KD 16,500 planes → 3 sets | §4.3.2–3, Table 2 |
| C5 | DSE major sets deviate < 5° from manual mapping | §4.3.1, p. 28 |
| C6 | GeoRank_3D = PCA normals, planarity filter, density-peak clustering, spatial patching, GRI, kinematics, Monte Carlo | §3.3.1–2, Eqs. 1–5 |
| C7 | GRI weights "stored as a single configurable vector" | §3.3.1, p. 20 |
| C8 | GeoRank_3D extracted 3 sets, all GRI = 0.73, threshold 0.50 | §4.6, Table 8 |
| C9 | Set 3 (80/017) reproduces field JS‑1 to 6.5° | Abstract, §4.6, Table 8 |
| C10 | No set mean is kinematically feasible | §4.6, Table 9 |
| C11 | 21.1 % of the Set 3–Set 1 wedge population is feasible | Abstract, §4.6, Table 9 |
| C12 | 26–29 % wedge feasibility in the mapped populations corroborates C11 | Abstract, §4.6 |
| C13 | Monte Carlo gives 95 % CIs on dip and dip direction (Fisher 1993) | §3.3.2, Table 7 |
| C14 | 31 blocks, 222.13 m³, log-normal, middle sector 74.59 % | §4.4, Table 5, Figs. 16–18 |
| C15 | V ≈ 0.20 A^1.5 proxy reproduces the field census | §3.3.3 Eq. 12, Table 10 |
| C16 | RocFall3: 27 m/s, >25 m bounce, 38–40 % reach the road, 15 m barrier retains 56–60 % | §4.5, Figs. 19–21 |

---

## 4. Paper–code correspondence table

Legend: **F** fully implemented · **P** partial · **D** implemented differently ·
**A** absent from code · **X** in code but not in paper · **N** not verifiable

| Claim / method | Status | Code location | Note |
|---|---|---|---|
| Eq. (1) point cloud | F | `_read_pointcloud.m` | |
| Eq. (2) covariance, k = 30 | F | `_estimate_normals.m:68–75` | |
| Eq. (3) planarity, surface variation | F | `_estimate_normals.m:98–100` | |
| Eq. (3) roughness = √λ₁ | F | `_estimate_normals.m:101` | dimensionally consistent (m) |
| Eq. (4) θ = arccos\|nᵢ·nⱼ\| | **D** | `_cluster_normals.m:65–68` | Euclidean chord on flipped normals, not the axial angle — §5.3 |
| "θ_eps default 8°" | **D** | `_default_options.m:14,39` | operative membership tolerance is `setConeDeg = 22°`, absent from the paper — §5.6 |
| Density-peak pole finding | F | `_cluster_normals.m:83–170` | mean-shift is axially correct |
| DBSCAN/k-means safeguard | F | `_cluster_normals.m:190–213` | |
| Spatial patch pass, 0.75 m / 80 pts | F | `_spatial_patches.m` | |
| Eq. (5) GRI | **P** | `_reliability_score.m:20` | weights hard-coded, contradicting C7 — §5.4 |
| Eq. (5) A_s persistence | **D** | `_compute_set_metrics.m:24,31` | fed from the set-wide hull, not max patch persistence as §3.3.3 states — §5.5 |
| "must be spatially coherent" | **A** | — | patch count/area computed, never enters the GRI (the paper admits this for inter-patch continuity but still states it as a condition) |
| Eq. (6)–(7) plane fit, RMS | F | `_fit_plane.m` | |
| Eq. (8) dip | F | `_orientation.m:7` | |
| Eq. (9) dip direction | **D/critical** | `_orientation.m:9`, `_dipdir2normal.m:5` | 180° error — §5.1 |
| Eq. (10) persistence | F | `_projected_area.m:16–19` | |
| Eq. (11) line of intersection | F | `_wedge_analysis.m:13–20` | trend convention correct |
| Eq. (12) V = 0.20 A^1.5 | F | `_compare_block_volumes.m:33` | GUI only |
| Markland planar / toppling | **P** | `_kinematic_index.m:19–20` | windows swapped by the Eq. (9) error — §5.1 |
| Oblique / flexural toppling | A (declared) | — | correctly declared out of scope |
| Dispersion-aware screening | F | `_kinematic_dispersion.m` | uses the *correct* azimuth; inconsistent with the rest — §5.2 |
| Monte Carlo, 50 runs, ±2°, ±25 %, 30° | F | `_uncertainty.m:31–74` | `mcAssignThresh_deg` never read |
| "95 % CI … Fisher (1993)" | **D** | `_uncertainty.m:93–99` | standard error of the mean across runs, linear for dip; no Fisher statistic — §5.8 |
| 10° validation vs 4 sources | **P** | `_validate_multi_source.m` | GUI only; deviations corrupted by Eq. (9) — §5.1 |
| KS test on log block volumes | F | `_compare_block_volumes.m:76` | GUI only; p-value never reported |
| Schmidt equal-area pole plot | F* | `_plot_poles.m:24–26` | *correct only because two 180° errors cancel — §5.7 |
| C2 manual mapping (4,818) | **A** | — | no code, no data |
| C3 DSE (1,472 → 10) | **A** | — | external tool, no export supplied |
| C4 qFacet FM/KD | **A** | — | external tool, no export supplied |
| C5 DSE within 5° of manual | **contradicted** | — | §5.9 |
| C14 block statistics, Figs. 16–18 | **A** | — | reproduced here in `georank3d/blockstats.py` |
| C16 RocFall3 | **A** | — | commercial software; no project file or export |
| Exposure-surface patch exclusion | **X** | `_compare_block_volumes.m:52–60` | ceiling = 5 × max observed; rule is in the paper, but its data-driven origin is not discussed as a selection effect |
| `_orientation_deviation.m` | **X** | dead code | header cites a "JRMGE Table 3" that does not exist |

---

## 5. Critical inconsistencies and unsupported claims

### 5.1 The dip direction is 180° out — **critical**

**Paper.** Eq. (9): "α = mod(atan2(−a, −b), 360°)", where §3.3.3 defines
`n = (a, b, c)` as the unit normal "oriented into the upper hemisphere".

**Code.** `functions/GeoRank_3D_orientation.m:9`
```matlab
dipdir = mod(atan2d(-nx, -ny), 360);
```
after `GeoRank_3D_unify_normals` has forced `nz ≥ 0`. Comment on line 8: "Down-dip
direction is opposite the horizontal projection of the upward normal." That
statement is the error. `functions/GeoRank_3D_dipdir2normal.m:5` inverts the
same wrong convention.

**Why it is wrong.** Three independent checks, all in `tests/test_geometry.py`:

* *Explicit surface.* A face descending eastward is `z = −tan δ · x`, i.e.
  `tan δ · x + z = 0`, so the upward normal is `(sin δ, 0, cos δ)` — pointing
  **east**, the dip direction. Eq. (9) returns 270°.
* *Orthogonality.* The down-dip vector is
  `(cos δ sin α, cos δ cos α, −sin δ)`. With `n = (+sin δ sin α, +sin δ cos α, cos δ)`
  the dot product is 0 to machine precision; with the Eq. (9) normal it is
  −0.94 — that vector is not in the plane at all.
* *Stereonet identity.* The lower-hemisphere pole plunges `90 − δ` towards
  `α + 180`; negating it gives exactly `(+sin δ sin α, +sin δ cos α, cos δ)`.

**Effect on reported results.**

*(a) Every GeoRank_3D attitude is rotated by 180°.* Table 8, Fig. 22 and Fig. 23:

| Reported | True attitude |
|---|---|
| Set 3, 80/017 | 80/197 |
| Set 2, 73/185 | 73/005 |
| Set 1, 72/208 | 72/028 |

*(b) Table 8's deviation columns are not the angles between the planes.* In
`GeoRank_3D_validate_multi_source.m:28,31` both sides go through
`GeoRank_3D_dipdir2normal`. For the extracted set this is a round trip that
returns the original PCA normal unchanged; for the reference set, quoted by a
human in the standard convention, it applies a genuine 180° rotation. Only one
vector is rotated, so the deviation is wrong. Recomputed with a single
convention:

| Set | Source | Table 8 | As coded | Corrected | 10° verdict |
|---|---|---|---|---|---|
| 3 (80/197) | Manual | 6.5 ✓ | 6.3 | **19.0** | pass → **fail** |
| 3 | qFacet FM / KD | 9.8 / 8.6 ✓ | 9.5 / 8.2 | **29.2 / 28.1** | pass → **fail** |
| 3 | DSE | 12.3 ✗ | 12.0 | **8.0** (JS‑3) | fail → **pass** |
| 2 (73/005) | Manual | 25.7 ✗ | 25.7 | **10.7** | fail → fail |
| 2 | qFacet FM / KD | 31.0 / 36.4 ✗ | 31.0 / 36.4 | **8.7 / 9.6** | fail → **pass** |
| 1 (72/028) | Manual | 31.3 ✗ | 30.9 | **5.1** (JS‑2) | fail → **pass** |
| 1 | qFacet FM / KD | 39.8 / 38.4 ✗ | 39.5 / 38.2 | **9.6 / 9.5** | fail → **pass** |

The paper's central validation claim reverses. **Set 3 no longer reproduces any
manually mapped set within 10°, and the best field agreement in the whole table
belongs to Set 1 — the population §4.6 dismisses as "a broad slope-parallel
exposure surface with no field counterpart closer than 20.9°".** The abstract's
"independently reproduced the dominant joint set (6.5° from field mapping)" and
the Conclusion's "6.5° agreement" are not supported.

*(c) The mean-based planar and toppling windows are interchanged.*
`GeoRank_3D_kinematic_index.m:15–20` compares the flipped set azimuth against
`opts.slopeDipDir`, a true real-world azimuth. At Chilas all six flags are 0
under either convention, so Table 9's mean column is unaffected; the code is
nevertheless wrong in general.

**What is *not* affected.** Because a 180° azimuth flip is a rotation about the
vertical, angles between two attitudes both expressed in the same convention are
preserved. The wedge dispersion branch uses the correct line-trend convention on
true PCA normals, so **the 21.1 % headline figure stands as computed** (subject
to §5.3).

**Caveat, stated plainly.** This conclusion assumes the exported cloud is
(easting, northing, elevation), as the code's own docstring and the paper's
"WGS84/UTM Zone 42N" imply. The authors can settle it in one run: the dominant
exposure-surface set must come out at the surveyed slope aspect, **80/135**.
Neither the reported 72/208 nor the corrected 72/028 is close to that, which is
itself an unexplained inconsistency with §4.6's description of Set 1 as
"slope-parallel" (see §5.10).

**Fix.** `matlab_patches/GeoRank_3D_orientation.m`,
`GeoRank_3D_dipdir2normal.m`, plus the plotting patch — apply all three
together (§5.7). **Recompute Table 8, Fig. 22, Fig. 23 and every orientation
quoted in §4.6, the abstract and the conclusions.**

---

### 5.2 Two contradictory orientation conventions coexist — **major**

`GeoRank_3D_kinematic_dispersion.m:106` (private `local_orient`) uses
`mod(atan2d(nx, ny), 360)` — the *correct* form — while
`GeoRank_3D_orientation.m:9` uses the flipped one. The dispersion module
duplicates the conversion instead of calling the shared function, so the two
disagree by 180° throughout the run. On a synthetic set with mean 70/140
screened against an 80/135 face at φ = 35°:

```
flipped azimuth : planar  0.50 %   toppling 93.27 %
correct azimuth : planar 74.79 %   toppling  4.17 %
```

At Chilas both give 0.0 %, so Table 9's population column is unaffected, but a
package with two mutually inconsistent orientation conventions cannot be
described as reproducible. **Severity: major.**
Fix: `matlab_patches/GeoRank_3D_kinematic_dispersion_patch_notes.m`.

---

### 5.3 Hemisphere-seam truncation of steeply dipping sets — **major**

`GeoRank_3D_cluster_normals.m:29,65–68` forces every normal to `nz ≥ 0` and then
assigns it to a pole by **Euclidean chord** with `knnsearch`. The pole of an 80°
plane plunges only 10°, so the `nz = 0` seam runs through the middle of every
steeply dipping set; members that dip past vertical are flipped to the antipode,
land ~160° from their own pole and are discarded as noise.

Simulated on Fisher populations (`tests/test_pipeline.py::test_clustering_does_not_truncate_a_near_vertical_set`):

| Set | κ | retained, as coded | retained, axial | mean dip, as coded | mean dip, axial |
|---|---|---|---|---|---|
| 80/017 | 60 | 90.2 % | 98.8 % | 78.77 | 80.02 |
| 80/017 | 25 | 70.2 % | 83.6 % | 77.29 | 79.99 |
| 73/185 | 25 | 81.3 % | 83.8 % | 72.41 | 73.00 |
| 72/208 | 25 | 82.3 % | 84.1 % | 71.63 | 72.06 |

Effect on the headline wedge statistic (80/017 × 72/208, κ = 40):

```
truncated population (as coded) : 20.6 % feasible, P95 plunge 54.1°
full axial population           : 20.8 % feasible, P95 plunge 60.9°
```

So the 21.1 % figure is robust to this defect, but the **95th-percentile plunge
of 63.2° reported in Table 9 is biased low by roughly 7°**, the per-set point
counts in Table 8 are understated by up to 30 % for the steepest set, and the
GRI abundance term and all Monte Carlo intervals inherit the truncation. The
mean-shift step at `_cluster_normals.m:121` already handles axial data correctly
via `sign(Nb*p.')`; only the membership and support tests are wrong.
Fix: `matlab_patches/GeoRank_3D_assign_to_poles.m`.

---

### 5.4 `opts.reliabilityWeights` is declared but never read — **major**

**Paper, §3.3.1:** "The five weights sum to unity and are stored as a single
configurable vector, allowing sensitivity analysis without modifying the code."

**Code:** `grep -rn reliabilityWeights` returns exactly one hit —
`GeoRank_3D_default_options.m:21`, the declaration. `GeoRank_3D_reliability_score.m:20`
hard-codes them:
```matlab
score = 0.20*np + 0.25*pl + 0.15*rf + 0.20*pe + 0.20*ks;
```
The claimed sensitivity analysis is therefore impossible without editing source,
and no weight-sensitivity result appears anywhere in the paper.
`opts.mcAssignThresh_deg` is dead in the same way (`_uncertainty.m:21` hard-codes
`thr = 30`), and `opts.useParallel` is declared "reserved" and never used.
Fix: `matlab_patches/GeoRank_3D_reliability_score.m`.

---

### 5.5 The persistence term saturates by construction — **major**

`_compute_set_metrics.m:24` computes `persistenceApprox` as the maximum convex-hull
chord of the **entire set** projected onto one best-fit plane — an outcrop
dimension of tens of metres, not a trace length. Line 31 passes that to the GRI,
where `pe = min(persistence/10, 1)` saturates at 1 for every set on every slope.
`patchPersistenceMax` — the quantity §3.3.3 says is used ("aggregated per set as
the total patch area and the maximum patch persistence") — is computed on line 45
and never passed.

Combined with `np = min(log10(N)/5, 1)`, also saturated for any set above 10⁵
points, Eq. (5) reduces at Chilas to `0.40 + 0.25 P_s + 0.15 R_s + 0.20 K_s`
with `K_s = 0` for all three sets. Hence Table 8's 0.730 / 0.730 / 0.728: **a
spread of 0.002 across three sets the paper itself describes as structurally
different**, and a GRI ≥ 0.50 threshold that nothing can fail. The claim that
GeoRank_3D ranks "by geomechanical significance rather than statistical
abundance" is not what Eq. (5) does as implemented — two of its five terms are
constants and a third is zero.
Fix: `matlab_patches/GeoRank_3D_compute_set_metrics_patch_notes.m`.

---

### 5.6 The operative clustering tolerance is not the one the paper states — **moderate**

§3.3.1 states θ_eps = 8° and that each normal is "assigned to its nearest pole
within this tolerance". In the code, `angularEpsDeg = 8` is used only by the
legacy DBSCAN fallback and as a floor for two derived defaults; actual set
membership uses `setConeDeg = 22°` (`_cluster_normals.m:35,65`), with a
`densityBandwidthDeg = 10°` mean-shift kernel, `mergeDeg = 8°` and
`maxSets = 12`. None of 22°, 10° or 12 appears in the manuscript. A ~3× change in
the membership cone materially changes point counts, set means and every
dispersion percentile. Fix: report the real parameters in a table and state
their sensitivity.

---

### 5.7 Two 180° errors cancel in the stereonet — **critical if half-fixed**

`_plot_poles.m:24–26` plots at radius `√2·sin(δ/2)` — correct Lambert equal-area —
but at azimuth `dipdir` rather than the pole trend `dipdir + 180`. Because
`dipdir` already carries the §5.1 error, the two cancel and **Fig. 22 is
geometrically correct**. Anyone fixing Eq. (9) alone will silently rotate every
stereonet by 180° without any error being raised. The same applies to
`GeoRank_3D_App.m:269`. Fix: `matlab_patches/GeoRank_3D_plot_poles_patch_notes.m`
— apply with patches 1 and 2.

---

### 5.8 The Monte Carlo interval is not what the paper says it is — **moderate**

§3.3.2 promises "a 95 % confidence interval on each set's dip and dip direction
from circular statistics (Fisher 1993)". `_uncertainty.m:93–99` computes
`mu ± 1.96·std/√n` across runs — a standard error of the **mean**, which shrinks
as `1/√monteCarloRuns` and can be made arbitrarily narrow by asking for more
runs. The dip interval uses linear statistics; no Fisher statistic (α₉₅, cone of
confidence, resultant length) appears anywhere. The Monte Carlo also perturbs
only the normals: it does not re-run normal estimation or the planarity filter,
so it propagates clustering sensitivity but not point-cloud or feature
uncertainty.

**More consequentially: no Monte Carlo interval is reported anywhere in the
manuscript.** Table 8 has no CI column; §4.6 quotes none. Claim C13 and the
"Monte Carlo 95 % CI, circular statistics (Fisher 1993)" row of Table 7 are
therefore **unsupported by any reported number**. Fix:
`matlab_patches/GeoRank_3D_uncertainty_patch_notes.m`; add the CIs to Table 8.

---

### 5.9 "DSE major sets deviate by less than 5° from the manual field measurements" — **contradicted by the paper's own Table 4**

§4.3.1, p. 28. Computed from Table 4 with the paper's own metric
(`tests/test_validation.py::test_dse_major_sets_are_not_within_5_degrees_of_manual_mapping`):

| Pair (Table 4) | Deviation | Nearest manual set |
|---|---|---|
| DSE JS‑1 68/017 vs manual JS‑1 82/011 | **15.1°** | JS‑2 at 11.1° |
| DSE JS‑2 40/002 vs manual JS‑2 67/029 | **34.3°** | JS‑2 at 34.3° |
| DSE JS‑3 81/189 vs manual JS‑3 78/346 | **31.0°** | JS‑1 at 17.1° |

Not one is within 5°; the dip alone differs by 27° for JS‑2, and DSE JS‑3 is
157° from manual JS‑3 in dip direction. Table 8's own DSE column (12.3° for the
best match) contradicts the same sentence. **Severity: major** — this sentence
supports the paper's argument that the three sets converge across methods.
Either delete it or replace it with the computed values.

---

### 5.10 Set 1 is described as "slope-parallel" but its reported attitude is not — **moderate**

§4.6 and the Fig. 23 caption call Set 1 (72/208) "the dominant slope-parallel
exposure surface". The slope face is 80/135 (§4.6, Table 9). 208° is 73° from
135° in azimuth; the corrected 028° is 107° away. Under neither convention is
Set 1 slope-parallel. Either the characterisation or the orientation is wrong;
this is also the cleanest diagnostic for settling §5.1.

---

### 5.11 qFacet FM wedge percentage is arithmetically impossible — **moderate**

§4.3.2, p. 29: "Wedge sliding affected about 27 % (1,872 of 475,503) of the
intersections". 1,872 / 475,503 = **0.39 %**. 1,872 is also the *plane* count
reported two sentences earlier, reused as an intersection count. The companion
figure "direct and oblique toppling involved 15 % and 12 % of the 274,513
potential intersections" is unverifiable. Fix: recompute from the Dips output or
delete the raw counts.

---

### 5.12 Random pairing measures orientation distributions, not wedges — **moderate**

`_kinematic_dispersion.m:76–78` pairs randomly drawn members of two sets. This
samples the joint orientation distribution under an independence assumption; no
spatial constraint requires the two normals to belong to surfaces that actually
intersect on the face, and points from the same physical wedge are not
independent. The 21.1 % is a legitimate statistic *about two orientation
distributions*, but the abstract's phrasing — "21.1 % of the critical Set 3–Set 1
wedge population" — reads as a frequency of physical wedges. Rephrase, and state
the independence assumption.

---

### 5.13 C12: the 26–29 % "corroboration" is not independent — **moderate**

The abstract and §4.6 present 21.1 % (GeoRank population screening) as
"independently corroborated" by 26–29 % from the mapped populations. The
authors' own `FIX_NOTES_dispersion.md` is candid that these "are different
populations… Report whatever the run produces; do not substitute the manual
number". They are computed by different tools on different populations against
different set definitions, and the two numbers differ by 5–8 percentage points
with no error bars on either. "Consistent in order of magnitude" is defensible;
"corroborates" is not.

---

### 5.14 `GeoRank_3D_KKH.m` cannot produce the paper's central result — **major (reproducibility)**

`GeoRank_3D_KKH.m` omits both `GeoRank_3D_kinematic_dispersion` and the Monte
Carlo cluster-tracking arguments (line 43 calls `_uncertainty` without
`refDips`/`refDipDirs`). It is the pipeline invoked by `run_GeoRank_3D_KKH.m`
*and* by `GeoRank_3D_demo_synthetic.m` — the only two scripts a reader would try
first, and the ones the README's Quick Start points at. A reader following the
README therefore cannot obtain the dispersion-aware result that carries the
paper's argument. Two near-duplicate pipelines that silently differ in scientific
content is itself a defect.

### 5.15 Table 8 and Table 10 are unreachable from any script — **major (reproducibility)**

`_validate_multi_source`, `_compare_block_volumes` and `_read_reference_csv` are
called only from `GeoRank_3D_App.m`. Neither CLI driver produces Table 8's
deviation columns or Table 10. A GUI-only path is not a reproducible pipeline.

### 5.16 Undocumented decimation of the input cloud — **moderate / not verifiable**

The paper states the VOM has ~65 M points (§3.1, §3.4) and that GeoRank_3D was
"applied to the same UAV-derived virtual outcrop model" (§3.3). But
`FIX_NOTES_clustering.md` refers to "the 10 M-point Chilas cloud, reduced to
~1.69 M by voxelisation", and Table 8's point counts sum to 1,239,900 — consistent
with ~1.69 M retained after the planarity filter, not with 65 M. The paper never
reports the working point count, the applied voxel size, or the retained
fraction. **Not verifiable from the available files**; these three numbers must
be stated.

### 5.17 The block-volume ceiling is set from the validation data — **minor**

`_compare_block_volumes.m:52`: `blockCeiling = 5 * max(obsVols)`. The comment
argues it "cannot be tuned to force agreement". The rule is defensible and
disclosed in §3.3.3, but the cut-off is nonetheless derived from the very census
being tested — selection on the comparison variable. Report the KS p-value
(currently printed to the console and never tabulated) and the sensitivity of
the agreement to the multiple. The corrected implementation reports both.

### 5.18 Minor items

* Seven line seeders described, six listed: "Seven rockfall line seeders (S1 and
  S2 to S6)" (§3.4, p. 25).
* §3 calls the block-volume method "voxel based differencing"; §3.4 and the
  abstract correctly call it 2.5D raster differencing. These are different
  methods.
* `_orientation_deviation.m` is dead code citing a non-existent "JRMGE Table 3".
* README lists item "9" twice; the limitations checklist skips item 9.
* `poleProminenceFrac` is documented in `FIX_NOTES_clustering.md` as a default
  option but is absent from `GeoRank_3D_default_options.m`.
* `GeoRank_3D_check_memory` is called by `GeoRank_3D.m` but not by `_KKH.m`.
* `datestr` (deprecated since R2022b) is used in two places.
* No licence file. For a paper claiming an "open analytical tools" workflow, the
  distribution terms of GeoRank_3D are undefined.

---

## 6. Code-quality and implementation audit

**Genuinely good.** The blocked normal estimation with the analytic symmetric
3×3 eigensolver (Smith 1961) is a correct and well-executed optimisation; the
identity `roughness = √λ_min` is right and saves a second neighbourhood pass.
Eqs. (2), (3), (6), (7), (10) and (11) are implemented correctly and are
dimensionally consistent. The bracketed uniform voxel coarsening in
`_prepare_points` is a thoughtful alternative to random thinning. The
subsample-and-propagate strategy for clustering and patching is a reasonable
scalability trade-off and is documented. The mean-shift pole search handles
axial data correctly. Numerical guards (`max(..., eps)`, degenerate-eigenvector
fallback) are present throughout.

**Weaknesses.**

* *No tests of any kind.* Not one assertion, not one expected output. For a
  package whose contribution is a numerical index, this is the single largest
  gap after the data.
* *Duplicated logic.* Two pipelines (`GeoRank_3D.m` / `_KKH.m`) that differ in
  scientific content; the orientation conversion is implemented twice, with
  opposite signs (§5.2). Divergence was inevitable and it happened.
* *Silent failure modes.* `_spatial_patches.m:22–24` catches any DBSCAN error and
  returns a single patch; `_projected_area.m:23` catches hull failures and
  returns zero area **and zero persistence**, which silently zeroes the GRI
  persistence term; `GeoRank_3D_App.m` has six bare `catch; end` blocks that
  discard plotting errors. None of these warn.
* *Options that do nothing:* `reliabilityWeights`, `mcAssignThresh_deg`,
  `useParallel`; `linearity` is computed for every point and never used (the
  paper does say it is diagnostic only).
* *No provenance.* The run report records point counts and slope parameters but
  not the software version, the seed, the git commit, or the full option struct.
  A result cannot be traced to the configuration that produced it.
* *`rng(opts.randomSeed)` is set*, which is good, but the pipeline is not
  deterministic across MATLAB versions because `dbscan`, `kmeans` and
  `knnsearch` tie-breaking are not version-stable, and nothing records the
  version used.
* *Input validation* is essentially absent: no check that a cloud is 3-D, that
  the slope parameters are in range, that weights sum to 1, or that dip and
  dip-direction columns in a reference CSV have not been transposed
  (`_read_reference_csv.m:39` matches any column containing "dip" and not "dir",
  so a file with a `DipDirection` column and no `Dip` column silently yields
  dip = dip direction).

**No evidence of data leakage or invalid train/test splitting was found** — the
workflow is not a supervised-learning pipeline, so those categories largely do
not apply. The nearest analogue is §5.17, where a comparison threshold is
derived from the validation census.

---

## 7. Data, experiment and reproducibility audit

| Requirement | Status |
|---|---|
| Dataset identity / source | Stated in the paper (§3.1); **not supplied** |
| Sample sizes | 65 M points claimed; **~1.7 M actually analysed** and never reported (§5.16) |
| Input variables | XYZ only; documented |
| Preprocessing | Manual cleaning + SOR in CloudCompare (§3.2.3) — **outside the code**, unparameterised, unrepeatable |
| Outlier handling | SOR "one standard deviation, 6 neighbours"; no code |
| Normalisation | Not applicable |
| Splitting / cross-validation | Not applicable |
| Data leakage | None found; see §5.17 for the nearest analogue |
| Reproducible pipeline | **No.** Hard-coded placeholder paths, GUI-only validation, no seeds recorded, no provenance, no tests, no data |
| Reference data for Tables 4 and 5 | **Not supplied** in any machine-readable form |

The manual CloudCompare stage is the weakest link for reproducibility: the
segmentation of intact rock from talus is described as done "manually on the
basis of spatial position and surface texture" and is acknowledged as "the most
important step" (§3.2.3). It cannot be repeated by anyone, and its effect on the
extracted sets is unquantified.

---

## 8. Results-reproduction assessment

**Attempted and succeeded — Section 4.4, Table 5, Figs. 16–18.** Reimplemented
from the Table 5 inventory (`georank3d/blockstats.py`; run
`scripts/run_block_statistics.py`). Every published value reproduces:

| Quantity | Paper | Recomputed |
|---|---|---|
| n blocks / total volume | 31 / 222.13 m³ | 31 / **222.13** |
| mean / median | 7.17 / 2.10 m³ | **7.165 / 2.100** |
| coefficient of variation | 1.49 | **1.4871** |
| D_eq = (6V/π)^⅓ | Table 5 | max discrepancy **0.005 m** |
| ln-mean / ln-sd | 0.94 / 1.53 | **0.9386 / 1.5333** |
| geometric mean | 2.56 m³ | **2.5563** |
| KS on ln V | D = 0.120, p = 0.722 | **D = 0.1196, p = 0.7222** |
| Shapiro–Wilk on ln V | p = 0.637 | **p = 0.6371** |
| normal Q-Q R² | 0.976 | **0.9764** |
| log-normal AIC / BIC | 176.66 / 179.53 | **176.66 / 179.53** (lowest of five) |
| Gini | 0.66 | **0.6626** |
| ≥ 20 / ≥ 5 / ≥ 1 m³ volume share | 56.65 / 87.3 / 98.2 % | **56.65 / 87.34 / 98.19** |
| Middle / Side / Lower share | 74.59 / 19.04 / 6.37 % | **74.59 / 19.04 / 6.37** |
| Kruskal–Wallis | p = 0.424 | **H = 1.7167, p = 0.4239** |

This section is arithmetically sound and was clearly computed carefully. **The
code for it is simply not in the archive**; it is supplied here.

**Attempted and failed — Table 8.** The deviation columns are reproducible to
0.4° *if both sides are read as quoted attitudes in a single convention*
(`tests/test_validation.py::test_table_8_deviations_are_reproduced`), which
confirms Table 8 is internally self-consistent. But that is not the comparison
the code performs (§5.1), so the tabulated values are not the true angular
deviations.

**Not attemptable.** Tables 1–4, 6, 9, 10; Figs. 9–15, 19–23. No point cloud, no
DSE/qFacet/Dips exports, no CloudCompare project, no RocFall3 model, and no
MATLAB licence. Table 9's kinematic feasibility and Table 10's proxy comparison
depend entirely on the 65 M-point cloud.

**Could the submitted code have produced the reported results?** For Table 8 and
Table 9, plausibly yes — the values are consistent with the code's own
(defective) formulas, and Table 8's near-identical GRIs are exactly what the
saturated Eq. (5) predicts, which is corroborating evidence that the numbers
came from this code. For Tables 5 and 6 and Figs. 16–18, **no** — no code exists
for them. For Figs. 19–21, no — RocFall3 is external.

**Corrected implementation, executed here.** The reference implementation runs
end-to-end in ~16 s on an 88,500-point synthetic slope with four known
orientations and recovers all four to within 0.1°:

```
 set         N    dip  dipdir    GRI  accept
   1     59732   80.0   135.0  0.683    True     <- slope face (ground truth 80/135)
   2      8959   67.0    29.0  0.652    True     <- JS-2        (ground truth 67/029)
   3      8974   78.1   345.9  0.643    True     <- JS-3        (ground truth 78/346)
   4      8865   81.9    11.0  0.639    True     <- JS-1        (ground truth 82/011)
```

Note what this shows about Eq. (5), on data where the answer is known: **the
slope face scores highest.** Once its point density is comparable to the joint
sets, the exposure surface is the most abundant, most planar and most persistent
population on the slope, so a GRI weighted 0.20/0.25/0.15/0.20 on those terms
ranks it first. Eq. (5) as specified cannot separate an exposure surface from a
structural discontinuity; that requires a term the index does not contain.

---

## 9. Scientific-validity assessment

**Novelty.** The claimed contribution — an explicit, reproducible per-set
reliability ranking layered on established extraction — is genuine as a *concept*
and the comparison in Table 7 is fair: DSE and qFacet really do not provide
kinematic screening, per-set ranking or built-in multi-source validation. But
the ranking as implemented returns 0.730/0.730/0.728 (§5.5). A ranking that
cannot separate a field-confirmed joint set from an exposure surface has not
been demonstrated, whatever its merit in principle.

**Fairness of comparison.** GeoRank_3D is not compared with DSE or qFacet on
equal terms: those tools were run in their own software by the authors with
tuned parameters (Tables 1–2), while GeoRank_3D's operative parameters (22° cone,
10° bandwidth) are not reported at all (§5.6). No common evaluation protocol, no
shared ground truth beyond three manual attitudes, no repeated runs, no
sensitivity analysis of the GRI weights (which was impossible — §5.4).

**Statistical support.** Section 4.4 is the strongest part of the paper and its
statistics are correct. Two caveats the paper does not state:

* The **74.59 % middle-sector share has a bootstrap 95 % CI of roughly 45–92 %**
  (n = 31), and the Side sector's is 3–45 % — the intervals overlap. Fig. 18b
  claims to show bootstrap CIs, but no interval is quoted in the text or the
  abstract, where 74.59 % is stated to four significant figures.
* The Kruskal–Wallis test (p = 0.424) shows no significant difference in median
  block size between sectors, which the paper reports honestly; the volume
  concentration is driven by four blocks.

**Uncertainty and repeated runs.** Monte Carlo machinery exists but no interval
is reported anywhere (§5.8). No repeated runs, no seed sensitivity, no parameter
sensitivity for GeoRank_3D. The DSE (kNN, η_MAX) and RocFall3 (R_n, μ)
sensitivity analyses are reported, but only as prose — no tables, no code, no
outputs.

**Overfitting.** Not applicable in the supervised sense. The analogous risk is
parameter tuning to a known answer: Table 2's qFacet parameters were "tuned
iteratively so that the extracted facets… agreed with the manually mapped sets",
which is disclosed and reasonable, but it means qFacet's agreement with manual
mapping is not independent evidence.

**Physical plausibility.** The 27 m/s and >25 m bounce heights are plausible for
a 36 m near-vertical face; the rigid-body, mass-independent framing is used
consistently and the fragmentation limitation is discussed unusually well
(§5, pp. 53–54). The 15 m barrier result follows.

**External validation.** None. Comparison with Ali et al. (2021) block volumes
is qualitative. The workflow is claimed to be "transferable" on the basis of a
single slope.

**Conclusions vs evidence.** Three conclusions are stronger than the evidence:
(i) "reproduced the dominant persistent joint set independently of the field
survey (6.5° agreement)" — does not survive §5.1; (ii) "ranking … by
geomechanical reliability" — the index does not rank (§5.5); (iii)
"corroborating 26–29 % from the mapped populations" — not an independent
corroboration (§5.13). The volume, trajectory and mitigation conclusions are
proportionate to their evidence.

**Classification of findings.**

* *Code defects:* §5.1, §5.2, §5.3, §5.4, §5.5, §5.7, §5.8, §6 silent failures.
* *Methodological weaknesses:* §5.6, §5.12, §5.17, unreported GRI weight
  sensitivity, no common benchmark protocol.
* *Unsupported paper claims:* C5 (§5.9), C7 (§5.4), C9 (§5.1), C13 (§5.8), the
  "ranking" framing (§5.5), C12 (§5.13), §5.11 arithmetic.
* *Missing reproducibility components:* all data, tests, Section 4.4 and 4.5
  code, GUI-only validation, undocumented decimation (§5.16), no licence.
* *Limitations that do not invalidate the study:* no fragmentation model,
  gravity-only triggering, manual CloudCompare cleaning, single-site validation,
  oblique/flexural toppling out of scope. All are acknowledged in the paper.

---

## 10. Overall verdict

**Insufficient or incomplete for validating the paper's claims**, and for the
GeoRank_3D component specifically, **substantially inconsistent with the paper**.

Evidence: the archive contains code for §3.3 only, and none for §4.1–4.5; the
central validation number (6.5°) rests on a comparison that mixes two coordinate
conventions (§5.1); the reliability index that is the paper's stated contribution
returns a 0.002 spread and cannot fail its own threshold (§5.5); the paper's
configurable-weights claim is contradicted by a one-line grep (§5.4); the DSE
"< 5°" claim is contradicted by the paper's own Table 4 (§5.9); no reported
number is reachable from any command-line script (§5.14, §5.15); and there are
no tests or data at all.

Against that, Section 4.4 reproduces exactly, the numerical core of the
algorithm (Eqs. 2, 3, 6, 7, 10, 11) is correctly implemented, and the
fragmentation and rigid-body limitations are discussed with real care.

### Direct answers

**1. Does the code implement the model described in the paper?**
Partly. Eqs. (1)–(4), (6), (7), (10), (11) yes. Eq. (9) is implemented as
published and as published it is wrong. Eq. (5) is implemented with hard-coded
weights and an input that saturates it. Eq. (12) is implemented but reachable
only from the GUI. The pipeline architecture matches §3.3.2 stage for stage.

**2. Does it support all major methodological claims?**
No. C5, C7, C9, C12 and C13 are not supported; the "reliability ranking" framing
is not supported by the index's behaviour.

**3. Can it reproduce the reported results?**
Section 4.4: yes, and this repository does so exactly — but from code the
authors did not supply. Tables 8–10 and Fig. 22: not without the point cloud,
and Table 8's numbers are defective independently of the data. Sections 4.1–4.3
and 4.5: no code exists.

**4. Are important components missing?**
Yes: all input data; all Section 4.4 statistics code; all Section 4.5 material;
tests; a CLI path to the validation layer; the RocFall3 model; a licence; and
the actual working point count of the GeoRank_3D run.

**5. Are the conclusions justified?**
The block-volume, trajectory and mitigation conclusions are. The three
GeoRank_3D-specific conclusions are not, in their present form.

**6. Would this be acceptable for a rigorous Q1 review?**
No — not in this state. My recommendation would be **major revision**: the
science is salvageable and the site work is valuable, but Table 8 and every
orientation in §4.6 must be recomputed, the reliability index must be shown to
discriminate or the claim withdrawn, and the archive must be completed.

**7. What must be corrected before submission?**
See Section 11.

---

## 11. Required corrections before publication

**Must fix — the paper cannot be assessed without these.**

1. Correct Eq. (9) and recompute every GeoRank_3D orientation, Table 8, Fig. 22,
   Fig. 23 and the associated text in the abstract, §4.6 and the conclusions
   (§5.1). Apply MATLAB patches 1–3 together (§5.7).
2. Confirm the point-cloud axis convention and report the slope face's own
   extracted attitude as a sanity check (§5.1, §5.10).
3. Either demonstrate that the GRI discriminates — by feeding it patch-level
   persistence, de-saturating the abundance term, and publishing a
   weight-sensitivity analysis — or withdraw the "ranking by reliability" framing
   (§5.4, §5.5).
4. Delete or correct the "< 5°" DSE claim (§5.9) and the "1,872 of 475,503"
   arithmetic (§5.11).
5. Report the working point count, applied voxel size and retained fraction of
   the GeoRank_3D run; reconcile with the 65 M-point VOM (§5.16).
6. Report the Monte Carlo confidence intervals in Table 8, or remove claim C13
   and the corresponding row of Table 7 (§5.8).
7. Report the real clustering parameters (22° cone, 10° bandwidth, 12 max sets)
   (§5.6).

**Should fix — required for the reproducibility claim to hold.**

8. Adopt axial statistics for cone membership and set means (§5.3), and remove
   the duplicate orientation convention (§5.2).
9. Merge `GeoRank_3D.m` and `GeoRank_3D_KKH.m`; expose validation and block
   comparison from the CLI (§5.14, §5.15).
10. Deposit the data: the decimated cloud actually analysed, the four reference
    set CSVs, the block inventory, the DSE/qFacet exports, the RocFall3 model and
    its outputs. Archive with a DOI and cite it; replace "available on
    reasonable request" (JRMGE and most Q1 journals now expect a repository).
11. Add tests. At minimum: orientation round-trip, a synthetic cloud with known
    set attitudes, GRI component reconstruction, and the Section 4.4 statistics
    against the published values. This repository supplies all four.
12. Record provenance in every run: version, seed, git commit, resolved options.
13. Add a licence.

**Nice to have.**

14. Report bootstrap CIs on the sector volume shares in the text, not only in a
    figure (§9).
15. Report the KS p-value for Table 10 and the sensitivity of the agreement to
    the exposure-surface ceiling (§5.17).
16. Rephrase the wedge-population statistic to reflect the independence
    assumption (§5.12) and soften "corroborates" (§5.13).
17. Fix the seeder count, the voxel/2.5D terminology, the duplicated README
    item, and remove the dead `_orientation_deviation.m` (§5.18).

---

## 12. Corrected implementation plan

Delivered in this repository:

* **`georank3d/`** — a corrected, executable reference implementation of §3.3,
  dependency-light (numpy, scipy, matplotlib), with every parameter in a
  validated JSON configuration, provenance written on every run, and one
  orientation convention used everywhere. Modules map one-to-one onto §3.3.2's
  nine stages.
* **`georank3d/blockstats.py`** — Section 4.4 and Figs. 16–18, absent from the
  archive, reproducing every published value.
* **`matlab_patches/`** — eight drop-in corrections so the authors can keep
  their MATLAB package and GUI.
* **`tests/`** — 62 tests, including regression guards for each defect above and
  a check of every Section 4.4 statistic against the published number.
* **`data/`** — Tables 4 and 5 transcribed to CSV so the validation and
  block-statistics stages are runnable. Clearly labelled as transcription, not
  as primary data.

**Deliberately not implemented.** RocFall3 (commercial third-party, no
reimplementation is faithful); the CloudCompare manual cleaning and 2.5D volume
computation (interactive, and the paper gives no reproducible parameterisation);
DSE and qFacet (external tools, already published). Nothing has been added that
the paper does not describe.

---

## 13. Complete revised code

See `README.md` for install and run instructions, `georank3d/` for the package,
`scripts/` for the two CLI entry points, and `matlab_patches/PATCHES.md` for the
MATLAB corrections. Every module docstring cites the manuscript section or
equation it implements and, where it differs from the archive, says why.

---

## 14. Information still required from the authors

1. The point cloud actually analysed by GeoRank_3D, and the reconciliation
   between 65 M, 10 M and ~1.7 M points (§5.16).
2. Confirmation of the coordinate frame of the exported cloud, and the
   GeoRank_3D attitude of the slope face itself (§5.1, §5.10).
3. The 4,818 manual orientations, the 1,472 DSE planes, and the qFacet FM/KD
   exports, as CSV.
4. The Dips project files behind Figs. 11–14 and the 26–29 % wedge figures.
5. The CloudCompare 2.5D volume projects behind Table 5, including the buffer
   width and the 90–99 % cell-matching criterion.
6. The RocFall3 model, seeder coordinates, block-shape assumptions, number of
   trajectories per seeder, and the raw outputs behind Figs. 19–21 and the
   "38–40 %" and "56–60 %" figures.
7. The GRI weight-sensitivity analysis that §3.3.1 implies was performed.
8. The Monte Carlo confidence intervals (§5.8).
9. The exact GeoRank_3D options struct used for the published run.
10. Whether Table 8's deviations were produced by
    `GeoRank_3D_validate_multi_source` or computed by hand.
11. The two-analyst manual mapping protocol: how disagreements were resolved and
    what the inter-analyst spread was (§3.2.1 claims this is "a robust ground
    truth"; no agreement statistic is reported).
12. Software versions: MATLAB release, toolbox versions, CloudCompare build,
    DSE 2.5.0.0 parameters as exported, RocFall3 version.

---

## 15. Recommended supplementary materials for a Q1 submission

* **S1 — Data archive with DOI** (Zenodo or similar): decimated point cloud,
  classified cloud, four reference set CSVs, block inventory, DSE/qFacet
  exports, Dips projects, RocFall3 model and outputs.
* **S2 — Code archive with DOI**, versioned and licensed, containing the
  corrected MATLAB package *and* a test suite with expected outputs.
* **S3 — Reproduction guide**: exact commands, the resolved options struct, the
  seed, runtime, and hardware, that regenerate Tables 8–10 and Figs. 16–18 and
  22–23 from S1.
* **S4 — Parameter-sensitivity tables**: GRI weights, `setConeDeg`,
  `densityBandwidthDeg`, `planarityMin`, `normalK`; the DSE kNN/η_MAX study
  currently described only in prose; the RocFall3 R_n/μ study likewise.
* **S5 — Uncertainty report**: Monte Carlo per-set dip and dip-direction
  intervals with the statistic defined, plus bootstrap CIs on the sector volume
  shares.
* **S6 — Manual-mapping protocol**: the two-analyst procedure and an
  inter-analyst agreement statistic.
* **S7 — Convention note**: an explicit statement of the coordinate frame, the
  dip/dip-direction convention and the hemisphere convention, with the
  slope-face attitude as a worked check. Given §5.1, this is the single most
  valuable supplement for a reader.

---

*Prepared as an independent audit. All recomputations are reproducible from this
repository: `python -m pytest` and the two scripts in `scripts/`.*
