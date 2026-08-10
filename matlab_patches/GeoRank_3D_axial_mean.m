function [axis, eigvals] = GeoRank_3D_axial_mean(N)
%GEORANK_3D_AXIAL_MEAN Principal axis of a set of axial unit vectors.
%
%   NEW HELPER (audit patch 7 of 8, part 1).
%   Discontinuity normals are AXIAL data: n and -n describe the same plane.
%   GeoRank_3D_compute_set_metrics used mean(Ni,1) on normals that had been
%   forced into the upper hemisphere by GeoRank_3D_unify_normals. That is not
%   invariant to the arbitrary sign of each normal, and it is biased for any
%   set whose pole lies near the nz = 0 seam. The pole of an 80-degree plane
%   plunges only 10 degrees, so every steeply dipping set at Chilas sits close
%   to that seam.
%
%   The correct estimator is the eigenvector of the largest eigenvalue of the
%   orientation tensor T = (1/n) * sum(n_i * n_i') (Fisher 1993, Ch. 3), which
%   is sign-invariant by construction.
%
%   [AXIS, EIGVALS] = GEORANK_3D_AXIAL_MEAN(N) returns the unit axis, oriented
%   into the upper hemisphere, and the three eigenvalues in descending order.
%   EIGVALS(1) near 1 indicates a tight cluster; the shape of the eigenvalue
%   spectrum is the standard measure of orientation concentration.

if isempty(N)
    axis = [0 0 1]; eigvals = zeros(1,3); return;
end
nn = sqrt(sum(N.^2,2));
nn(nn == 0) = 1;
N  = N ./ nn;

T = (N.' * N) / size(N,1);
[V, D] = eig((T + T.')/2, 'vector');
[eigvals, ord] = sort(D, 'descend');
axis = V(:, ord(1)).';
if axis(3) < 0
    axis = -axis;
end
eigvals = eigvals(:).';
end
