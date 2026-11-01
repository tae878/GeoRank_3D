# Fix: GeoRank_3D detected only one discontinuity set

## Symptom
On real slope clouds (e.g. the 10 M-point Chilas cloud, reduced to ~1.69 M by
voxelisation), the run finished quickly and reported only **Set 1** — far worse
than DSE (10 sets) or qFacet (3 sets), so the reliability ranking had nothing to
rank.

## Root cause
`GeoRank_3D_cluster_normals.m` clustered the unit normals with a **single global
epsilon DBSCAN**. On an exposed rock face the normals form a *continuous*
distribution: the large slope-face population angularly bridges the joint-set
poles. With one epsilon (8 deg), DBSCAN's density-reachability chains every mode
into a **single cluster**. This is a well-documented failure of DBSCAN on
directional data — nothing to do with the point cloud or the parameters.

## Fix
The clustering step now finds the **principal poles as peaks of the orientation
density on the sphere** (greedy mean-shift), then assigns every normal to its
nearest pole within a cone — the same family of method DSE (Riquelme et al.,
2014) uses. A prominence test suppresses diffuse noise peaks; `minClusterPts`
and the downstream GRI screening do the rest. The function signature is
unchanged, so the GUI and the rest of the pipeline are unaffected.

Validated on synthetic Chilas-like data: recovers the slope face **plus JS-1
(82/11), JS-2 (67/29) and JS-3 (78/346)** exactly, and generalises to 3-set,
5-set and no-dominant-face cases, deterministically across random seeds, with no
spurious noise sets.

## New tunable options (defaults in GeoRank_3D_default_options.m)
    opts.densityBandwidthDeg = 10;   % mode-seeking kernel half-width [deg]
    opts.setConeDeg          = 22;   % max pole-to-normal angle for set membership
    opts.mergeDeg            = 8;    % merge poles closer than this [deg]
    opts.maxSets             = 12;   % cap on number of sets
    opts.poleProminenceFrac  = 0.03; % min peak density vs largest peak (noise guard)

If a run ever finds fewer than two poles, it automatically falls back to the
legacy DBSCAN/k-means path, so behaviour never regresses below the old code.

## What you will now see in the Command Window / Run Log
    pole 1: dip .. / dipdir ..  (support ..)
    pole 2: dip .. / dipdir ..  (support ..)
    ...
    Density-peak clustering: N discontinuity set(s) retained.
      set 1: .. points
      ...
