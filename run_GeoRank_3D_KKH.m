%% run_GeoRank_3D_KKH.m
% Example driver script for a real Chilas / KKH point cloud.
% Replace inputFile with your exported point cloud path.

clear; clc;
addpath(genpath(fileparts(mfilename('fullpath'))));

inputFile = 'your_point_cloud.xyz';   % XYZ/TXT/CSV/PLY; first three columns must be X Y Z
outDir = fullfile(pwd, 'GeoRank_3D_output');

opts = GeoRank_3D_default_options();
opts.maxPoints = 10000000;     % Analysis budget (default 10M). Memory-bounded; lower for first tests
opts.memBudgetGB = 9;          % Peak working-set target in GB (tune to your RAM; 9 suits 16 GB)
opts.voxelSize = 0.05;         % metres; adjust to point density/GSD
opts.normalK = 30;
opts.planarityMin = 0.20;
opts.surfaceVariationMax = 0.18;
opts.angularEpsDeg = 8;
opts.minClusterPts = 300;
opts.patchRadius = 0.75;       % metres; spatial grouping scale
opts.minPatchPts = 80;
opts.slopeDip = 80;
opts.slopeDipDir = 135;        % update to actual slope dip direction
opts.frictionAngle = 35;
opts.monteCarloRuns = 50;
opts.randomSeed = 42;

results = GeoRank_3D_KKH(inputFile, outDir, opts);

disp(results.setTable);
