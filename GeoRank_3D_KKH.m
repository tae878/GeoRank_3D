function results = GeoRank_3D_KKH(inputFile, outDir, opts)
%GEORANK_3D_KKH Reliability-aware discontinuity extraction from 3D point clouds.
%
%   results = GeoRank_3D_KKH(inputFile, outDir, opts)
%
%   This function provides a full discontinuity extraction workflow:
%   import -> downsample -> normal estimation -> orientation clustering ->
%   patch extraction -> geomechanical significance ranking -> uncertainty.
%
%   The code is intended for research use and must be validated against field
%   structural measurements before engineering design.

if nargin < 2 || isempty(outDir); outDir = fullfile(pwd, 'GeoRank_3D_output'); end
if nargin < 3 || isempty(opts); opts = GeoRank_3D_default_options(); end
if ~exist(outDir, 'dir'); mkdir(outDir); end
rng(opts.randomSeed);

fprintf('GeoRank_3D started.\n');
fprintf('Reading point cloud: %s\n', inputFile);
[P, attr] = GeoRank_3D_read_pointcloud(inputFile);
N0 = size(P,1);
fprintf('Loaded %d points.\n', N0);

[P, sampleInfo] = GeoRank_3D_prepare_points(P, opts);
fprintf('Working cloud: %d points.\n', size(P,1));

fprintf('Estimating PCA normals and local geometric descriptors...\n');
[normals, geom] = GeoRank_3D_estimate_normals(P, opts);
valid = geom.planarity >= opts.planarityMin & geom.surfaceVariation <= opts.surfaceVariationMax & all(isfinite(normals),2);
Pvalid = P(valid,:);
Nvalid = normals(valid,:);
geomValid = GeoRank_3D_filter_struct(geom, valid);
fprintf('Planarity-filtered points: %d / %d.\n', size(Pvalid,1), size(P,1));

fprintf('Clustering normal vectors on unit sphere...\n');
labels = GeoRank_3D_cluster_normals(Nvalid, opts);

fprintf('Computing discontinuity set metrics...\n');
[setTable, patchTable] = GeoRank_3D_compute_set_metrics(Pvalid, Nvalid, labels, geomValid, opts);
wedgeTable = GeoRank_3D_wedge_analysis(setTable, opts);

fprintf('Running uncertainty assessment...\n');
uncertaintyTable = GeoRank_3D_uncertainty(Pvalid, Nvalid, geomValid, opts);

fprintf('Exporting results...\n');
GeoRank_3D_export_results(outDir, P, valid, labels, setTable, patchTable, uncertaintyTable, opts, sampleInfo, attr);
if ~isempty(wedgeTable); writetable(wedgeTable, fullfile(outDir, 'GeoRank_3D_wedge_pairs.csv')); end
GeoRank_3D_plot_poles(outDir, Nvalid, labels, setTable, opts);

results = struct();
results.inputFile = inputFile;
results.outDir = outDir;
results.options = opts;
results.originalPointCount = N0;
results.workingPointCount = size(P,1);
results.validPointCount = size(Pvalid,1);
results.setTable = setTable;
results.patchTable = patchTable;
results.uncertaintyTable = uncertaintyTable;
results.wedgeTable = wedgeTable;
results.sampleInfo = sampleInfo;

fprintf('GeoRank_3D complete.\n');
end
