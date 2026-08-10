% GeoRank_3D_kinematic_dispersion.m -- audit patch 4 of 8 (edit in place)
%
% The private helper local_orient (lines 102-107) duplicates the orientation
% conversion instead of calling GeoRank_3D_orientation, and uses the opposite
% azimuth sign:
%
%   GeoRank_3D_orientation.m:9   ddir = mod(atan2d(-nx, -ny), 360)
%   local_orient (line 106)      ddir = mod(atan2d( nx,  ny), 360)
%
% local_orient is the geologically correct form (see patch 1); the rest of the
% pipeline is the one that is 180 degrees out. Until patch 1 is applied, the
% two coexist, so the dispersion-aware planar-sliding and direct-toppling
% windows are effectively interchanged relative to every other module.
%
% On a synthetic set with mean 70/140 (Fisher kappa = 25, n = 20000) screened
% against a slope face of 80/135 with phi = 35 degrees, the two conventions
% give:
%       flipped azimuth : planar  0.50 %,  toppling 93.27 %
%       correct azimuth : planar 74.79 %,  toppling  4.17 %
%
% At Chilas both conventions happen to return 0.0 per cent for planar and
% toppling on all three sets, so the numbers printed in Table 9 are not
% affected, but the code is wrong in general.
%
% The wedge branch (lines 79-87) is unaffected: trend = atan2d(Lx, Ly) of the
% downward-oriented intersection vector is the correct trend convention for a
% LINE, and it is applied to the true member normals. The 21.1 per cent Set 3
% - Set 1 wedge fraction quoted in the abstract and Table 9 therefore stands
% as computed, subject only to the seam truncation of patch 6.
%
% FIX: delete local_orient and call the shared function, so there is exactly
% one orientation convention in the package.
%
%   line 54:  [dip, ddir] = local_orient(Ns);
%   becomes:  [dip, ddir] = GeoRank_3D_orientation(Ns);
%
% then delete lines 102-107. GeoRank_3D_orientation is already vectorised over
% rows and already unifies the hemisphere, so the behaviour is otherwise
% identical.
%
% SEPARATE ISSUE, same file, lines 76-78:
%       m  = min([size(Na,1), size(Nb,1), nSamp]);
%       Aa = Na(randi(size(Na,1),m,1),:);
%       Bb = Nb(randi(size(Nb,1),m,1),:);
% Random pairing samples the joint distribution of the two orientation
% populations under an assumption of INDEPENDENCE between them. Points from
% the same physical wedge are not independent, and no spatial constraint
% requires the two sampled normals to belong to surfaces that actually
% intersect on the face. The resulting fraction is therefore a property of the
% two orientation distributions, not a frequency of real wedges. That is a
% legitimate quantity, but it must be described as such in the manuscript
% rather than as "21.1 per cent of the critical Set 3 - Set 1 wedge
% population", which reads as a population of physical wedges.
