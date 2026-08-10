# Drop-in corrections for the submitted MATLAB release

These files are minimal, drop-in replacements for functions in
`03_GeoRank_3D_Code_and_GUI.zip`. They are provided so the authors can keep
their MATLAB package and App Designer GUI while fixing the defects identified
in `../AUDIT_REPORT.md`. Copy each file over the one of the same name in
`GeoRank_3D/functions/` (or `GeoRank_3D/` for the top-level scripts).

Every patch preserves the original function signature, so the GUI and the rest
of the pipeline are unaffected.

**Apply patches 1, 2 and 3 together.** `GeoRank_3D_plot_poles.m` currently
contains a second 180-degree error that cancels the first, so the published
stereonet is geometrically correct even though the tabulated attitudes are
not. Fixing only the orientation convention would silently rotate every
stereonet by 180 degrees.

| # | file | defect | severity |
|---|---|---|---|
| 1 | `GeoRank_3D_orientation.m` | dip direction returned 180 degrees out (manuscript Eq. 9 applied to an upward normal) | critical |
| 2 | `GeoRank_3D_dipdir2normal.m` | inverse of the same wrong convention | critical |
| 3 | `GeoRank_3D_plot_poles.m` | plots the pole at the dip-direction azimuth instead of the pole trend | critical when 1 and 2 are fixed |
| 4 | `GeoRank_3D_kinematic_dispersion.m` | private `local_orient` duplicated the orientation code with the opposite sign | major |
| 5 | `GeoRank_3D_reliability_score.m` | `opts.reliabilityWeights` declared but never read; weights were literals | major |
| 6 | `GeoRank_3D_cluster_normals.m` | Euclidean chord test on hemisphere-flipped normals truncates steeply dipping sets at the `nz = 0` seam | major |
| 7 | `GeoRank_3D_compute_set_metrics.m` | set mean from the arithmetic mean of flipped normals; GRI persistence fed from the set-wide convex hull | major |
| 8 | `GeoRank_3D_uncertainty.m` | `opts.mcAssignThresh_deg` declared but never read; CI is a standard error of the mean across runs | moderate |

## Not patched here, but required before submission

* `GeoRank_3D_KKH.m` does not call `GeoRank_3D_kinematic_dispersion`, so the
  dispersion-aware screening that carries the manuscript's central claim is
  unreachable from `run_GeoRank_3D_KKH.m` and from
  `GeoRank_3D_demo_synthetic.m`. Either delete `GeoRank_3D_KKH.m` and route
  both drivers through `GeoRank_3D.m`, or add the two missing calls.
* The multi-source validation (`GeoRank_3D_validate_multi_source`), the block
  comparison (`GeoRank_3D_compare_block_volumes`) and the reference-CSV reader
  are reachable only from the GUI. Manuscript Table 8 and Table 10 therefore
  cannot be produced from any command-line driver in the archive.
* `GeoRank_3D_orientation_deviation.m` is dead code: nothing calls it. Its
  header claims it produces "JRMGE Table 3", which does not correspond to any
  table in the manuscript.
* There is no analysis code at all for Sections 4.4 and 4.5 (block-volume
  statistics, magnitude-frequency, sector analysis, RocFall3). The Python
  package in this repository supplies the Section 4.4 analysis; RocFall3 is
  commercial third-party software and cannot be reimplemented here.
