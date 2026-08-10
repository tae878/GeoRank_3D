function n = GeoRank_3D_dipdir2normal(dip_deg, dipdir_deg)
%GEORANK_3D_DIPDIR2NORMAL Dip/dip-direction to upward unit normal vector.
%
%   CORRECTION (audit patch 2 of 8).
%   Inverse of the corrected GeoRank_3D_orientation. The previous release
%   negated the horizontal components, which is the normal of a plane dipping
%   towards dipdir + 180.
%
%   This matters beyond bookkeeping: in GeoRank_3D_validate_multi_source the
%   extracted set enters as its own PCA normal (because dipdir2normal inverts
%   orientation exactly) while each externally supplied reference attitude is
%   converted with this function. With the old signs only the reference was
%   rotated, so the tabulated angular deviation was not the angle between the
%   two planes.

d = deg2rad(dip_deg(:));
a = deg2rad(dipdir_deg(:));
n = [sin(d).*sin(a), sin(d).*cos(a), cos(d)];
rn = sqrt(sum(n.^2,2));
n  = n ./ max(rn, eps);
end
