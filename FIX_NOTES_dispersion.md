# Upgrade: dispersion-aware (research-gap-driven) kinematic screening

## What was limited
The kinematic tests screened only each set's MEAN orientation. On the Chilas run
every mean-based flag is zero (no planar, toppling or feasible wedge). That is a
correct result for the means — but the paper's argument is that hazard lives in the
DISPERSED TAILS of the orientation populations, and the code did not actually
compute the tails; the 26-29% wedge feasibility came only from the external Dips
analysis.

## What was added (authentic, not tuning)
`GeoRank_3D_kinematic_dispersion.m` screens the FULL orientation population of each
set (its real member normals from the cloud) and the pairwise intersections sampled
across those populations, and reports the FRACTION that is kinematically feasible.
No threshold is changed (friction angle, daylight windows, plunge criterion are all
identical). It measures where feasibility actually resides.

New outputs (printed to the log and written to CSV):
  GeoRank_3D_kinematic_dispersion_sets.csv
     SetID, N, PlanarFeasFrac, TopplingFeasFrac, PlanarFeasMean, TopplingFeasMean
  GeoRank_3D_kinematic_dispersion_wedges.csv
     SetA, SetB, WedgeFeasFrac, WedgeFeasMean, MedianPlunge_deg, P95Plunge_deg

Each row contrasts the population fraction (Frac) with the mean-orientation flag
(Mean), so the paper can state "mean 0% vs population X%" directly from GeoRank_3D
output rather than from a separate tool.

## Honesty note
The population fractions are computed from GeoRank_3D's OWN three sets
(80/017, 73/185, 72/208), so the resulting percentage is GeoRank's internal
measurement and need NOT equal the 26-29% obtained from the manual Dips analysis of
JS-1/2/3 — they are different populations and provide two independent lines of
evidence for the same conclusion. Report whatever the run produces; do not
substitute the manual number for the GeoRank number.

## Option
opts.dispSampleCap = 20000;  % samples per set for the population screening
