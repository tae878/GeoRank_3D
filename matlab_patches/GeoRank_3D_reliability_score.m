function [score, terms] = GeoRank_3D_reliability_score(nPts, planarity, roughness, persistence, kin, opts)
%GEORANK_3D_RELIABILITY_SCORE Geomechanical Reliability Index, Eq. (5).
%
%   CORRECTION (audit patch 5 of 8).
%   1. The five weights are now read from opts.reliabilityWeights. The field
%      existed in GeoRank_3D_default_options but was never referenced: the
%      weights were literals on the score line, so the manuscript's statement
%      that they "are stored as a single configurable vector, allowing
%      sensitivity analysis without modifying the code" did not hold.
%   2. The weights are checked to be non-negative and to sum to unity, as
%      Eq. (5) requires.
%   3. The five terms are returned so a reported GRI can be audited component
%      by component. On slope-scale sets the abundance and persistence terms
%      both saturate at 1, which is why the published Chilas run returned
%      0.730 / 0.730 / 0.728 for three sets the manuscript itself describes as
%      structurally different.
%
%   [SCORE, TERMS] = ... returns TERMS with fields abundance, planarity,
%   roughness, persistence, kinematic and two saturation flags.

w = [0.20 0.25 0.15 0.20 0.20];
if nargin >= 6 && isfield(opts,'reliabilityWeights') && ~isempty(opts.reliabilityWeights)
    w = opts.reliabilityWeights(:).';
end
if numel(w) ~= 5
    error('GeoRank_3D:weights', ...
        ['opts.reliabilityWeights must have five entries ' ...
         '[abundance planarity roughness persistence kinematic]; got %d.'], numel(w));
end
if any(w < 0) || abs(sum(w) - 1) > 1e-9
    error('GeoRank_3D:weights', ...
        ['opts.reliabilityWeights must be non-negative and sum to 1 ' ...
         '(Eq. 5 is a normalised index); got sum = %.6f.'], sum(w));
end

rScale = 0.10;
if nargin >= 6 && isfield(opts,'roughnessScale') && opts.roughnessScale > 0
    rScale = opts.roughnessScale;
end
pRef = 10;
if nargin >= 6 && isfield(opts,'persistenceReference') && opts.persistenceReference > 0
    pRef = opts.persistenceReference;
end
aScale = 5;
if nargin >= 6 && isfield(opts,'abundanceLogScale') && opts.abundanceLogScale > 0
    aScale = opts.abundanceLogScale;
end

np = min(log10(max(nPts,1)) / aScale, 1);
pl = min(max(planarity,0), 1);
rf = exp(-max(roughness,0) / rScale);
pe = min(max(persistence,0) / pRef, 1);
ks = min(max(kin.score,0), 1);

score = w(1)*np + w(2)*pl + w(3)*rf + w(4)*pe + w(5)*ks;
score = min(max(score,0), 1);

terms = struct('abundance',np, 'planarity',pl, 'roughness',rf, ...
               'persistence',pe, 'kinematic',ks, ...
               'abundanceSaturated', np >= 1-1e-12, ...
               'persistenceSaturated', pe >= 1-1e-12);
end
