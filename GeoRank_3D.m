function results = GeoRank_3D(inputFile, outDir, opts)
%GEORANK_3D  Reliability-aware discontinuity set extraction from 3D point clouds.
%   Original GeoRank_3D pipeline — bugs fixed, name updated to GeoRank_3D.
%   All internal GeoRank_3D_ functions preserved exactly.
%
%   USAGE (command line):  edit run_GeoRank_3D.m, set 5 parameters, press F5
%   USAGE (GUI):           type  GeoRank_3D_App
%
%   results = GeoRank_3D(inputFile, outDir, opts)

if nargin < 2 || isempty(outDir); outDir = fullfile(pwd,'GeoRank_3D_output'); end
if nargin < 3 || isempty(opts);   opts   = GeoRank_3D_default_options();       end
if ~exist(outDir,'dir'); mkdir(outDir); end
rng(opts.randomSeed);

fprintf('GeoRank_3D started.\n');

% 1. Read
[P, attr] = GeoRank_3D_read_pointcloud(inputFile);
N0 = size(P,1);
fprintf('Loaded %d points.\n', N0);

% 2. Preprocess
[P, sampleInfo] = GeoRank_3D_prepare_points(P, opts);
fprintf('Working cloud: %d points (voxel %.4g m', size(P,1), sampleInfo.appliedVoxelSize);
if ~strcmp(sampleInfo.budgetMethod,'none')
    fprintf(', reduced-to-budget: %s', sampleInfo.budgetMethod);
end
fprintf(').\n');
GeoRank_3D_check_memory(size(P,1), opts);

% 3. Normals (blocked, analytic eigensolver)
fprintf('Estimating normals...\n');
[normals, geom] = GeoRank_3D_estimate_normals(P, opts);

% 4. Quality filter
valid  = geom.planarity >= opts.planarityMin & ...
         geom.surfaceVariation <= opts.surfaceVariationMax & ...
         all(isfinite(normals),2);
Pv     = P(valid,:);
Nv     = normals(valid,:);
geomV  = GeoRank_3D_filter_struct(geom, valid);
fprintf('Retained %d / %d points.\n', size(Pv,1), size(P,1));

% 5. Cluster
fprintf('Clustering...\n');
labels = GeoRank_3D_cluster_normals(Nv, opts);

% 6. Set + patch metrics
fprintf('Computing set metrics...\n');
[setTable, patchTable] = GeoRank_3D_compute_set_metrics(Pv, Nv, labels, geomV, opts);

% 7. Wedge analysis (mean-orientation) + dispersion-aware (population) screening
wedgeTable = GeoRank_3D_wedge_analysis(setTable, opts);
[setDispTable, wedgeDispTable] = GeoRank_3D_kinematic_dispersion(Nv, labels, setTable, opts);
if ~isempty(setDispTable)
    fprintf('Dispersion-aware kinematics (population vs mean):\n');
    for ii = 1:height(setDispTable)
        fprintf('  Set %d: planar %.1f%% / toppling %.1f%% of population feasible (mean: %d / %d)\n', ...
            setDispTable.SetID(ii), 100*setDispTable.PlanarFeasFrac(ii), ...
            100*setDispTable.TopplingFeasFrac(ii), setDispTable.PlanarFeasMean(ii), ...
            setDispTable.TopplingFeasMean(ii));
    end
    if ~isempty(wedgeDispTable)
        for ii = 1:height(wedgeDispTable)
            fprintf('  Wedge %d-%d: %.1f%% of sampled intersections feasible (mean: %d)\n', ...
                wedgeDispTable.SetA(ii), wedgeDispTable.SetB(ii), ...
                100*wedgeDispTable.WedgeFeasFrac(ii), wedgeDispTable.WedgeFeasMean(ii));
        end
    end
end

% 8. Monte Carlo uncertainty (with cluster tracking)
fprintf('Monte Carlo uncertainty (%d runs)...\n', opts.monteCarloRuns);
refDips=[]; refDipDirs=[];
if ~isempty(setTable)
    refDips    = setTable.MeanDip_deg;
    refDipDirs = setTable.MeanDipDirection_deg;
end
[uncertaintyTable, setCI] = GeoRank_3D_uncertainty(Pv, Nv, geomV, opts, refDips, refDipDirs);

% 9. Export
GeoRank_3D_export_results(outDir, P, valid, labels, setTable, patchTable, ...
    uncertaintyTable, opts, sampleInfo, attr);
if ~isempty(wedgeTable)
    writetable(wedgeTable, fullfile(outDir,'GeoRank_3D_wedge_pairs.csv'));
end
if ~isempty(setDispTable)
    writetable(setDispTable, fullfile(outDir,'GeoRank_3D_kinematic_dispersion_sets.csv'));
end
if ~isempty(wedgeDispTable)
    writetable(wedgeDispTable, fullfile(outDir,'GeoRank_3D_kinematic_dispersion_wedges.csv'));
end

% 10. Figures (Schmidt equal-area)
GeoRank_3D_plot_poles(outDir, Nv, labels, setTable, opts);

% 11. Results struct
results = struct('inputFile',inputFile,'outDir',outDir,'options',opts,...
    'originalN',N0,'workingN',size(P,1),'validN',size(Pv,1),...
    'setTable',setTable,'patchTable',patchTable,...
    'uncertaintyTable',uncertaintyTable,'setCI',setCI,...
    'wedgeTable',wedgeTable,'setDispTable',setDispTable,...
    'wedgeDispTable',wedgeDispTable,'sampleInfo',sampleInfo);

fprintf('GeoRank_3D complete. Output: %s\n', outDir);
end
