% GeoRank_3D_compute_set_metrics.m -- audit patch 7 of 8 (edit in place)
%
% Two independent problems.
%
% (a) Line 23:  [dip, dipdir] = GeoRank_3D_orientation(mean(Ni,1));
%
%     mean(Ni,1) is the arithmetic mean of normals already forced into the
%     upper hemisphere. Normals are axial data, so this is not invariant to
%     the arbitrary sign of each vector and it is biased for any set whose
%     pole lies near the nz = 0 seam -- which is every steeply dipping set at
%     Chilas. Use the orientation-tensor principal axis instead:
%
%         axis = GeoRank_3D_axial_mean(Ni);
%         [dip, dipdir] = GeoRank_3D_orientation(axis);
%
%     GeoRank_3D_axial_mean.m is supplied in this directory. On a Fisher
%     population about 88/017 with kappa equivalent to 0.10 rad scatter, the
%     arithmetic mean is displaced by more than 1 degree while the tensor axis
%     stays within 0.5 degree.
%
% (b) Line 31:  reliability = GeoRank_3D_reliability_score(size(Pi,1), ...
%                   meanPlanarity, meanRoughness, persistenceApprox, kin, opts);
%
%     persistenceApprox comes from line 24,
%         [areaApprox, persistenceApprox] = GeoRank_3D_projected_area(Pi, plane.normal);
%     which is the maximum chord of the convex hull of the WHOLE SET projected
%     onto a single best-fit plane. A set spans the entire outcrop, so this is
%     an outcrop dimension of tens of metres, not a discontinuity trace
%     length. Against the 10 m reference in Eq. (5) the persistence term A_s
%     therefore saturates at 1 for every set, on every slope, by construction.
%
%     The manuscript itself says persistence is "aggregated per set as ... the
%     maximum patch persistence" (Section 3.3.3). That quantity is already
%     computed on line 45 as patchPersistenceMax but is never passed to the
%     score. Pass it instead:
%
%         reliability = GeoRank_3D_reliability_score(size(Pi,1), ...
%             meanPlanarity, meanRoughness, patchPersistenceMax, kin, opts);
%
%     This requires moving the reliability call BELOW the patch loop
%     (currently lines 38-48), since patchPersistenceMax is not known until
%     the loop has finished.
%
%     Note the ordering consequence: with A_s saturated and the abundance term
%     O_s = min(log10(N)/5, 1) also saturated for any set above 100,000 points,
%     the published Chilas GRI reduces to 0.40 + 0.25*P_s + 0.15*R_s + 0.20*K_s
%     with K_s = 0 for all three sets (Table 9). That is why Table 8 reports
%     0.730, 0.730 and 0.728: a spread of 0.002 across three sets the paper
%     describes as one field-confirmed joint set, one DSE-corroborated set and
%     one non-structural exposure surface. The index as implemented does no
%     discriminating work at this site, and the GRI >= 0.50 threshold is never
%     binding.
%
% (c) Minor, line 11:
%         patchRows = cell(ceil(size(P,1)/max(opts.minPatchPts,1)) + nSets, 11);
%     For a 1.2-million-point cloud with minPatchPts = 80 this pre-allocates
%     about 15,000 rows x 11 columns of cells whether or not they are used.
%     Harmless, but the comment claims it is a fix for a growth pattern; it is
%     an over-allocation.
