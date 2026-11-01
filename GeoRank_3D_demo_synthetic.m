function results = GeoRank_3D_demo_synthetic()
%GEORANK_3D_DEMO_SYNTHETIC Demonstrates GeoRank_3D on a synthetic rock slope.
%   The synthetic cloud contains three discontinuity families plus noise.

rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));
outDir = fullfile(rootDir, 'examples', 'demo_output');
if ~exist(outDir, 'dir'); mkdir(outDir); end
rng(7);

P = [];
colors = [];
planes = {
    [0.15 0.05 1.00], [0 0 0], 4500;  % low-dip bedding-like surface
    [1.00 0.20 0.15], [0 0 0], 3500;  % steep set
    [-0.30 1.00 0.20], [0 0 0], 3000; % steep set
};
for i = 1:size(planes,1)
    n = planes{i,1}; n = n ./ norm(n);
    c = planes{i,2}; N = planes{i,3};
    B = null(n);
    uv = (rand(N,2)-0.5) .* [20 15];
    noise = 0.03*randn(N,1);
    Pi = c + uv(:,1).*B(:,1)' + uv(:,2).*B(:,2)' + noise.*n;
    P = [P; Pi]; %#ok<AGROW>
    colors = [colors; i*ones(N,1)]; %#ok<AGROW>
end
P = [P; 25*(rand(1000,3)-0.5)];
colors = [colors; zeros(1000,1)]; %#ok<NASGU>

inputFile = fullfile(outDir, 'synthetic_rock_slope.xyz');
writematrix(P, inputFile, 'Delimiter', 'space');

opts = GeoRank_3D_default_options();
opts.maxPoints = inf;
opts.voxelSize = 0;
opts.normalK = 25;
opts.angularEpsDeg = 9;
opts.minClusterPts = 150;
opts.patchRadius = 1.2;
opts.minPatchPts = 40;
opts.slopeDip = 80;
opts.slopeDipDir = 135;
opts.frictionAngle = 35;
opts.monteCarloRuns = 20;
opts.randomSeed = 7;

results = GeoRank_3D_KKH(inputFile, outDir, opts);
fprintf('Demo complete. Outputs written to: %s\n', outDir);
end
