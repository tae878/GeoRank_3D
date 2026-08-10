% GeoRank_3D_uncertainty.m -- audit patch 8 of 8 (edit in place)
%
% Two changes; both are small edits rather than a rewrite, so the patch is
% given as instructions to avoid shipping a near-duplicate of the original.
%
% (a) line 21:  thr = 30;
%     The Monte Carlo cluster-to-set assignment threshold is hard-coded even
%     though opts.mcAssignThresh_deg exists in GeoRank_3D_default_options and
%     the manuscript quotes 30 degrees as a default rather than a constant.
%
%     replace with:
%         thr = 30;
%         if isfield(opts,'mcAssignThresh_deg') && ~isempty(opts.mcAssignThresh_deg)
%             thr = opts.mcAssignThresh_deg;
%         end
%
% (b) lines 93-99: the reported interval is
%         mu +/- 1.96 * std(dv(ok)) / sqrt(sum(ok))
%     i.e. a standard error of the MEAN across Monte Carlo runs. That interval
%     shrinks as 1/sqrt(monteCarloRuns), so it describes the precision of the
%     averaging, not the stability of the orientation, and it can be made
%     arbitrarily narrow by asking for more runs. The dip interval also uses
%     linear statistics while the manuscript cites Fisher (1993).
%
%     Report instead the percentile interval of the run-to-run orientation and
%     the Fisher cone of 95 per cent confidence about the mean axis:
%         setCI.dipCI95(rf,:)    = prctile(dv(ok), [2.5 97.5]);
%         c = dd(ok) - circ_mean_deg(dd(ok));
%         c = mod(c + 180, 360) - 180;
%         setCI.dipDirCI95(rf,:) = mod(circ_mean_deg(dd(ok)) + prctile(c,[2.5 97.5]), 360);
%     with alpha_95 = acosd(1 - ((n-R)/R) * ((1/0.05)^(1/(n-1)) - 1)) computed
%     from the matched unit axes, where R is the resultant length.
%
% Note also that the Monte Carlo perturbs only the normals. It does not re-run
% normal estimation or the planarity filter, so it propagates clustering
% sensitivity but not point-cloud or feature-extraction uncertainty. That
% limitation should be stated in the manuscript.
%
% Finally: no Monte Carlo confidence interval appears anywhere in the
% manuscript. Table 8 has no CI column and Section 4.6 quotes none, so the
% uncertainty claim in Section 3.3.2 and Table 7 is currently unsupported by
% any reported number.
