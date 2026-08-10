function labels = GeoRank_3D_assign_to_poles(N, poles, coneDeg)
%GEORANK_3D_ASSIGN_TO_POLES Axial assignment of normals to orientation poles.
%
%   REPLACEMENT for lines 64-68 of GeoRank_3D_cluster_normals.m
%   (audit patch 6 of 8).
%
%   The submitted code assigned normals with
%       coneChord = 2*sin(deg2rad(coneDeg)/2);
%       [nn, d]   = knnsearch(poles, N, 'K', 1);
%       labels    = nn;  labels(d > coneChord) = 0;
%   after GeoRank_3D_unify_normals had forced every normal to nz >= 0. That is
%   a Euclidean chord test between hemisphere-flipped vectors, which has a
%   discontinuity at the nz = 0 seam. For a set dipping 80 degrees the pole
%   plunges only 10 degrees, so the seam runs through the middle of the set
%   and the members that dip past vertical are flipped to the antipode, land
%   about 160 degrees from their own pole and are discarded as noise.
%
%   Simulation on Fisher-distributed populations about 80/017 gives, for a
%   22-degree cone: 90.2 per cent of the true membership retained at
%   kappa = 60 and 70.2 per cent at kappa = 25, versus 98.8 and 83.6 per cent
%   under the axial test, with the set mean dip biased low by up to 2.7
%   degrees and the 95th-percentile intersection plunge biased low by about
%   7 degrees. The reported N points per set, the GRI abundance term and
%   every dispersion percentile inherit that truncation.
%
%   The correct membership test is the axial angle arccos(|n . p|), which is
%   exactly manuscript Eq. (4) and is what the rest of the paper says is used.
%
%   Usage inside GeoRank_3D_cluster_normals:
%       labels = GeoRank_3D_assign_to_poles(N, poles, coneDeg);
%       labels = drop_small_and_compact(labels, minPts);
%
%   The same substitution is required in find_poles: build the density and
%   support queries on [N; -N] (the antipodally symmetrised sample) so a pole
%   near the equator collects its whole population. The mean-shift update in
%   the submitted code already handles this correctly with its sign(Nb*p.')
%   term, so the pole locations themselves are sound; only the membership and
%   support tests are affected.

if isempty(poles)
    labels = zeros(size(N,1),1);
    return;
end

nn = sqrt(sum(N.^2,2)); nn(nn == 0) = 1;
N  = N ./ nn;
pn = sqrt(sum(poles.^2,2)); pn(pn == 0) = 1;
poles = poles ./ pn;

ang = real(acosd(min(abs(N * poles.'), 1)));      % axial angle, Eq. (4)
[best, idx] = min(ang, [], 2);
labels = idx;
labels(best > coneDeg) = 0;
end
