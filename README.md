# GeoRank_3D: Reliability-Aware Discontinuity Extraction for Rockfall Source-Zone Analysis

GeoRank_3D is a research-grade MATLAB prototype prepared for discontinuity analysis of UAV/photogrammetric virtual outcrop models. It is designed as an improved, reproducible alternative to conventional DSE-style workflows for complex rock slopes such as the Chilas slope along the Karakoram Highway.

## Why this tool is different

Conventional semi-automated discontinuity extraction tools commonly detect many statistically distinct clusters, but not all clusters are geomechanically meaningful for slope-scale failure. GeoRank_3D adds a failure-analysis layer by linking discontinuity extraction with reliability metrics, persistence, roughness, kinematic screening, and uncertainty assessment.

## Main modules

1. Point-cloud import from TXT/CSV/PLY.
2. Optional voxel/grid downsampling for very large point clouds.
3. Adaptive PCA normal estimation and planarity filtering.
4. Pole conversion to dip/dip-direction.
5. Angular DBSCAN clustering on the unit normal sphere.
6. Patch-level spatial grouping for persistence and roughness estimation.
7. Geomechanical significance index based on persistence, planarity, roughness, and kinematic feasibility.
8. Pairwise wedge-intersection screening using mean discontinuity-set orientations.
9. Monte Carlo uncertainty assessment for cluster stability.
9. Export of classified point cloud, discontinuity-set table, uncertainty table, and diagnostic plots.

## Quick start

```matlab
addpath(genpath('GeoRank_3D'));
results = GeoRank_3D_demo_synthetic();
```

For a real point cloud:

```matlab
opts = GeoRank_3D_default_options();
opts.maxPoints = 10000000;      % Analysis budget (default 10M). Memory-bounded pipeline; see Performance section
opts.voxelSize = 0.05;          % Adjust to your point density
opts.slopeDip = 80;             % Chilas slope approximate value
opts.slopeDipDir = 135;         % Set based on your slope aspect
opts.frictionAngle = 35;
opts.angularEpsDeg = 8;
opts.minClusterPts = 300;

results = GeoRank_3D_KKH('your_point_cloud.xyz', 'GeoRank_3D_output', opts);
```

## Required MATLAB toolboxes

Recommended: Statistics and Machine Learning Toolbox, Computer Vision Toolbox, Image Processing Toolbox. The core functions include fallbacks where possible, but DBSCAN and nearest-neighbour search are best with Statistics and Machine Learning Toolbox.

## Performance & hardware scaling

The pipeline now defaults to an **analysis budget of 10,000,000 points** (`opts.maxPoints = 10000000`) and is built to process 8–10M-point clouds on a mid-range desktop without running out of memory. This replaces the previous 250k default.

What changed to make this possible:

- **Memory-bounded normal estimation.** Normals are computed in auto-sized blocks (controlled by `opts.memBudgetGB`, default 9 GB) using a single multithreaded k-d tree and a vectorized analytic 3×3 symmetric eigensolver instead of a per-point `eig()` loop. Roughness comes directly from the smallest eigenvalue, so neighbourhoods are visited once.
- **Scalable clustering & patching.** Angular DBSCAN and per-set spatial DBSCAN run directly below a size cap (`opts.clusterSampleCap`, `opts.patchSampleCap`); above the cap, a representative subsample defines the sets/patches and labels are propagated to all points by nearest prototype/core. Results stay stable while staying tractable at millions of points.
- **Uniform reduce-to-budget.** When a cloud exceeds `maxPoints`, it is coarsened by **uniform voxel down-sampling** (preserves spatial/orientation statistics) rather than random subsampling. Random subsampling remains only as a fallback (`opts.uniformSubsample = false` to force it).
- **Separate display budget.** Pole/scatter plots decimate to `opts.maxPlotPoles` (default 150,000) so rendering stays smooth on modest GPUs; this never affects the analysis.
- **Advisory memory pre-flight.** `GeoRank_3D_check_memory` estimates the resident footprint and prints a warning if it approaches available RAM (never errors).

### Tuned for the reference machine (Intel i5-6500, 16 GB RAM, GTX 1060 6GB)

The defaults below are chosen for exactly this class of machine:

```matlab
opts.maxPoints      = 10000000;  % 8–10M analysis budget
opts.memBudgetGB    = 9;         % peak working set target (leaves headroom on 16 GB)
opts.normalBlockSize= 0;         % 0 = auto-size blocks from memBudgetGB
opts.useParallel    = false;     % 4 physical cores, no HT; knnsearch is already multithreaded
opts.clusterSampleCap = 200000;  % DBSCAN sample size before propagation
opts.patchSampleCap   = 150000;
opts.maxPlotPoles     = 150000;  % display decimation for the GTX 1060
```

Practical guidance for 16 GB machines:

- Keep other applications closed; the target peak working set is ~9 GB, leaving room for the OS and MATLAB itself.
- Output is streamed to disk; on a mechanical HDD the classified-point-cloud export of a 10M-point run is the slowest step. Writing to the SSD is recommended.
- If you hit memory pressure, lower `opts.memBudgetGB` (e.g. to 6–7) or reduce `opts.maxPoints`; the block sizer adapts automatically.
- For clouds well beyond ~15–20M points, increase `opts.voxelSize` so the reduce-to-budget step coarsens to ~10M before analysis, or tile the cloud spatially and merge set tables.

## Important note

This is a research prototype. Validate outputs against field measurements, manual mapping, scar locations, and VOM interpretation before using the results for engineering design or journal submission.
