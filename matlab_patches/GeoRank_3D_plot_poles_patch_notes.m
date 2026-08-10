% GeoRank_3D_plot_poles.m -- audit patch 3 of 8 (edit in place)
%
% APPLY THIS ONLY TOGETHER WITH PATCHES 1 AND 2.
%
% The submitted file contains a second 180-degree error that cancels the one
% in GeoRank_3D_orientation, so the PUBLISHED STEREONET IS CORRECT even though
% the tabulated attitudes are not. Fixing the orientation convention on its own
% would silently rotate Fig. 22 by 180 degrees.
%
% Why the current code happens to be right:
%   * On a lower-hemisphere equal-area net the pole to a plane dipping delta
%     towards alpha plunges 90 - delta towards alpha + 180, so it plots at
%     radius sqrt(2)*sin(delta/2) in the direction alpha + 180.
%   * Lines 24-26 plot at radius sqrt(2)*sind(dip/2) -- correct -- in the
%     direction of the value returned by GeoRank_3D_orientation, which before
%     patch 1 is alpha + 180. Right position, for the wrong reason.
%
% After patches 1 and 2, GeoRank_3D_orientation returns the true alpha, so the
% plotting azimuth must be advanced by 180 degrees explicitly.
%
% Line 21-26, replace:
%       [dip, dipdir] = GeoRank_3D_orientation(Nplot);
%       r = sqrt(2) .* sind(dip ./ 2);
%       x = r .* sind(dipdir);
%       y = r .* cosd(dipdir);
% with:
%       [dip, dipdir] = GeoRank_3D_orientation(Nplot);
%       trend = mod(dipdir + 180, 360);      % lower-hemisphere pole trend
%       r = sqrt(2) .* sind(dip ./ 2);
%       x = r .* sind(trend);
%       y = r .* cosd(trend);
%
% Lines 55-57, the set-mean markers, replace:
%       ri = sqrt(2) * sind(setTable.MeanDip_deg(ii)/2);
%       xi = ri * sind(setTable.MeanDipDirection_deg(ii));
%       yi = ri * cosd(setTable.MeanDipDirection_deg(ii));
% with:
%       ti = mod(setTable.MeanDipDirection_deg(ii) + 180, 360);
%       ri = sqrt(2) * sind(setTable.MeanDip_deg(ii)/2);
%       xi = ri * sind(ti);
%       yi = ri * cosd(ti);
%
% The same two-line change is needed in GeoRank_3D_App.m, method drawPole,
% line 269:
%       r = sqrt(2).*sind(d./2); xc = r.*sind(a); yc = r.*cosd(a);
% becomes
%       t = mod(a + 180, 360);
%       r = sqrt(2).*sind(d./2); xc = r.*sind(t); yc = r.*cosd(t);
