%RUN_GEORANK_3D  Command-line driver for GeoRank_3D.
%   Edit the 5 mandatory inputs below, then press F5.
%   For GUI: type  GeoRank_3D_App

clc; clear; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% ---- MANDATORY INPUTS (edit these 5 lines) -------------------------
inputFile     = 'your_pointcloud.xyz';
outDir        = 'GeoRank_3D_output';
slopeDip      = 80;
slopeDipDir   = 135;
frictionAngle = 35;

%% ---- OPTIONAL: adjust parameters -----------------------------------
opts = GeoRank_3D_default_options();
opts.slopeDip      = slopeDip;
opts.slopeDipDir   = slopeDipDir;
opts.frictionAngle = frictionAngle;
% opts.angularEpsDeg   = 8;
% opts.minClusterPts   = 300;
% opts.voxelSize       = 0.05;
% opts.monteCarloRuns  = 50;
% opts.planarKinWin    = 25;
% opts.topplingKinWin  = 30;
% opts.maxPoints       = 10000000;  % analysis budget (default 10M); lower for first tests
% opts.memBudgetGB     = 9;         % peak working-set target in GB (9 suits 16 GB RAM)

%% ---- RUN -----------------------------------------------------------
results = GeoRank_3D(inputFile, outDir, opts);

fprintf('\n========== RESULTS ==========\n');
disp(results.setTable(:,{'SetID','MeanDip_deg','MeanDipDirection_deg',...
    'PlanarFeasible','TopplingFeasible','GeoMechanicalReliability'}));
fprintf('Output: %s\n', outDir);
