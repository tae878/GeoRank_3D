function [dip, dipdir] = GeoRank_3D_orientation(normals)
%GEORANK_3D_ORIENTATION Convert plane normals to dip and dip direction.
%   Assumes X = east, Y = north, Z = up.
%
%   CORRECTION (audit patch 1 of 8).
%   The previous release computed
%       dipdir = mod(atan2d(-nx, -ny), 360)
%   from a normal it had already forced into the UPPER hemisphere, following
%   Eq. (9) of the manuscript. That returns the dip direction rotated by 180
%   degrees. The horizontal projection of the UPWARD normal points DOWN dip:
%   for a surface descending eastwards, z = -tan(d)*x, the upward normal is
%   (sin d, 0, cos d), which points east. It is the DOWNWARD (lower
%   hemisphere) pole whose trend is dipdir + 180.
%
%   Verify with the down-dip vector, which must lie in the plane:
%       downdip = [cosd(d)*sind(a), cosd(d)*cosd(a), -sind(d)]
%       n       = [sind(d)*sind(a), sind(d)*cosd(a),  cosd(d)]
%       n * downdip' == 0
%
%   Apply this patch together with GeoRank_3D_dipdir2normal.m and
%   GeoRank_3D_plot_poles.m; see PATCHES.md.

N  = GeoRank_3D_unify_normals(normals);
nx = N(:,1); ny = N(:,2); nz = N(:,3);
dip    = atan2d(sqrt(nx.^2 + ny.^2), abs(nz));    % Eq. (8), unchanged
dipdir = mod(atan2d(nx, ny), 360);                % corrected Eq. (9)
end
