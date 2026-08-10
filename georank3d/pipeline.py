"""End-to-end GeoRank_3D pipeline (manuscript Section 3.3.2, stages 1-9).

    import and preprocess -> normals and features -> planarity filter ->
    orientation clustering -> spatial patches -> geometric parameters ->
    kinematic screening (mean and dispersion-aware) -> reliability ranking ->
    Monte Carlo uncertainty -> export.

Everything the run depends on is written to the output directory: the resolved
configuration, the package version, the random seed, per-stage point counts and
all result tables, so a run can be repeated exactly from its own output.
"""

from __future__ import annotations

import csv
import json
import platform
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List

import numpy as np

from . import __version__
from .clustering import cluster_normals, spatial_patches
from .config import Config
from .geometry import fit_plane, normal_from_dip, projected_hull
from .kinematics import (dispersion_screening, kinematic_score, planar_toppling_flags,
                         wedge_dispersion, wedge_pairs)
from .pointcloud import estimate_normals, read_point_cloud, reduce_to_budget
from .reliability import reliability_index
from .uncertainty import monte_carlo_uncertainty
from .validation import block_volume_comparison, compare_to_references, read_reference_sets

__all__ = ["run_pipeline", "write_csv"]


def write_csv(path, rows: List[Dict]) -> None:
    """Write a list of flat dicts to CSV; a header-only file if rows is empty."""
    path = Path(path)
    if not rows:
        path.write_text("")
        return
    keys: List[str] = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in keys})


def _log(cfg, message: str) -> None:
    if cfg.verbose:
        print(message, flush=True)


def run_pipeline(cfg: Config) -> Dict:
    """Run the full workflow and return a results dictionary."""
    cfg.validate()
    t0 = time.time()
    out_dir = Path(cfg.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(cfg.random_seed)

    # ---- 1. import and preprocess ------------------------------------
    _log(cfg, f"GeoRank_3D {__version__}: reading {cfg.input_cloud}")
    P_raw = read_point_cloud(cfg.input_cloud)
    P, sample_info = reduce_to_budget(P_raw, cfg)
    _log(cfg, f"  {sample_info['input_count']} points -> {sample_info['working_count']} "
              f"(voxel {sample_info['applied_voxel_size_m']:.4g} m, "
              f"budget: {sample_info['budget_method']})")

    # ---- 2. normals and per-point descriptors ------------------------
    _log(cfg, f"  estimating normals (k = {cfg.normal_k})")
    feats = estimate_normals(P, cfg)

    # ---- 3. planarity / surface-variation filter ---------------------
    keep = (feats.planarity >= cfg.planarity_min) & \
           (feats.surface_variation <= cfg.surface_variation_max) & \
           np.isfinite(feats.normals).all(axis=1)
    Pv, fv = P[keep], feats.subset(keep)
    _log(cfg, f"  retained {keep.sum()} / {len(P)} points after planarity filtering "
              f"({100 * keep.mean():.1f}%)")
    if keep.sum() < cfg.min_cluster_points:
        raise RuntimeError(
            f"only {int(keep.sum())} points survived the planarity filter, fewer than "
            f"min_cluster_points = {cfg.min_cluster_points}. Relax planarity_min "
            f"(currently {cfg.planarity_min}) or surface_variation_max "
            f"(currently {cfg.surface_variation_max}).")

    # ---- 4. orientation clustering ------------------------------------
    _log(cfg, "  density-peak clustering of normals")
    labels, _ = cluster_normals(fv.normals, cfg, rng)
    set_ids = [int(s) for s in np.unique(labels[labels > 0])]
    _log(cfg, f"  {len(set_ids)} orientation set(s) retained")
    if not set_ids:
        raise RuntimeError(
            "no orientation set met min_cluster_points. Lower min_cluster_points, "
            "widen set_cone_deg, or check that the cloud has been cleaned of "
            "vegetation and talus before import.")

    # ---- 5-8. per-set geometry, patches, kinematics, reliability ------
    set_records: List[Dict] = []
    patch_rows: List[Dict] = []
    patch_id = 0
    for sid in set_ids:
        m = labels == sid
        Pi, Ni = Pv[m], fv.normals[m]
        from .clustering import set_orientation
        dip, dipdir, axis = set_orientation(Ni)
        plane = fit_plane(Pi)
        hull_area, hull_persistence = projected_hull(Pi, axis)

        patches = spatial_patches(Pi, cfg, rng)
        patch_area_total, patch_persistence_max, n_patches = 0.0, 0.0, 0
        for p in np.unique(patches[patches > 0]):
            ip = patches == p
            if ip.sum() < cfg.min_patch_points:
                continue
            n_patches += 1
            patch_id += 1
            pl = fit_plane(Pi[ip])
            area, persistence = projected_hull(Pi[ip], pl["normal"])
            patch_area_total += area
            patch_persistence_max = max(patch_persistence_max, persistence)
            patch_rows.append(dict(patch_id=patch_id, set_id=sid, n_points=int(ip.sum()),
                                   dip_deg=pl["dip"], dipdir_deg=pl["dipdir"],
                                   plane_rmse_m=pl["rmse"], projected_area_m2=area,
                                   persistence_m=persistence,
                                   centroid_x=float(pl["centroid"][0]),
                                   centroid_y=float(pl["centroid"][1]),
                                   centroid_z=float(pl["centroid"][2])))

        if n_patches == 0:
            _log(cfg, f"    warning: set {sid} yielded no spatial patch of at least "
                      f"{cfg.min_patch_points} points within {cfg.patch_radius_m} m. Its "
                      f"persistence term is 0 by default, which suppresses its GRI. This is "
                      f"usually a point-density artefact: raise patch_radius_m or lower "
                      f"min_patch_points for sparsely sampled sets.")

        persistence = patch_persistence_max if cfg.persistence_source == "max_patch" else hull_persistence
        planar, toppling, azi, opp = planar_toppling_flags(dip, dipdir, cfg)
        kscore = float(kinematic_score(planar, toppling, cfg)[0])
        gri = reliability_index(int(m.sum()), float(np.mean(fv.planarity[m])),
                                float(np.mean(fv.roughness[m])), persistence, kscore, cfg)

        rec = dict(set_id=sid, n_points=int(m.sum()), dip_deg=float(dip), dipdir_deg=float(dipdir),
                   axis=axis, set_plane_rmse_m=plane["rmse"],
                   mean_planarity=float(np.mean(fv.planarity[m])),
                   mean_roughness_m=float(np.mean(fv.roughness[m])),
                   set_hull_area_m2=hull_area, set_hull_persistence_m=hull_persistence,
                   n_patches=n_patches, patch_area_total_m2=patch_area_total,
                   max_patch_persistence_m=patch_persistence_max,
                   persistence_used_m=float(persistence),
                   persistence_source=cfg.persistence_source,
                   azimuth_difference_deg=float(azi[0]), anti_slope_difference_deg=float(opp[0]),
                   slope_face_deviation_deg=float(np.degrees(np.arccos(min(abs(float(
                       axis @ normal_from_dip(cfg.slope_dip_deg, cfg.slope_dipdir_deg))), 1.0)))),
                   planar_feasible_mean=int(planar[0]), toppling_feasible_mean=int(toppling[0]),
                   kinematic_score=kscore, reliability_weights=list(cfg.reliability_weights))
        rec.update(gri)
        rec["accepted"] = bool(rec["gri"] >= cfg.reliability_threshold)
        set_records.append(rec)

    set_records.sort(key=lambda r: r["gri"], reverse=True)

    wedge_mean = wedge_pairs(set_records, cfg)

    normals_by_set = {}
    for sid in set_ids:
        idx = np.flatnonzero(labels == sid)
        if len(idx) > cfg.dispersion_sample_cap:
            idx = rng.choice(idx, cfg.dispersion_sample_cap, replace=False)
        normals_by_set[sid] = fv.normals[idx]
    set_dispersion = dispersion_screening(normals_by_set, set_records, cfg)
    wedge_disp = wedge_dispersion(normals_by_set, wedge_mean, cfg, rng)

    # ---- 9. Monte Carlo uncertainty -----------------------------------
    _log(cfg, f"  Monte Carlo uncertainty ({cfg.monte_carlo_runs} runs)")
    mc = monte_carlo_uncertainty(fv.normals, set_records, cfg, rng)

    # ---- validation against independent sources ------------------------
    validation = {}
    if cfg.reference_sets:
        sources = {name: read_reference_sets(path) for name, path in cfg.reference_sets.items()}
        validation = compare_to_references(set_records, sources, cfg.match_threshold_deg)
        write_csv(out_dir / "validation_extracted_to_reference.csv", validation["extracted_to_reference"])
        write_csv(out_dir / "validation_reference_to_extracted.csv", validation["reference_to_extracted"])

    block_check = {}
    if cfg.observed_blocks_csv:
        from .blockstats import read_block_inventory
        obs = [r["volume_m3"] for r in read_block_inventory(cfg.observed_blocks_csv)]
        block_check = block_volume_comparison([r["projected_area_m2"] for r in patch_rows], obs, cfg)
        write_csv(out_dir / "block_volume_proxy_sensitivity.csv", block_check["sensitivity"])

    # ---- export --------------------------------------------------------
    tabular_sets = [{k: v for k, v in r.items() if k not in ("axis", "reliability_weights")}
                    for r in set_records]
    write_csv(out_dir / "set_metrics.csv", tabular_sets)
    write_csv(out_dir / "patch_metrics.csv", patch_rows)
    write_csv(out_dir / "wedge_pairs_mean.csv", wedge_mean)
    write_csv(out_dir / "kinematic_dispersion_sets.csv", set_dispersion)
    write_csv(out_dir / "kinematic_dispersion_wedges.csv", wedge_disp)
    write_csv(out_dir / "uncertainty_runs.csv", mc["runs"])
    write_csv(out_dir / "uncertainty_sets.csv", mc["per_set"])

    if cfg.export_classified_cloud:
        full = np.zeros(len(P), dtype=int)
        full[keep] = labels
        np.savetxt(out_dir / "classified_cloud.xyz",
                   np.column_stack([P, full]), fmt="%.6f %.6f %.6f %d")

    if cfg.make_figures:
        from .plotting import pole_plot, reliability_bar
        pole_plot(out_dir / "pole_plot_schmidt.png", fv.normals, labels, set_records,
                  cfg.max_plot_poles, cfg.random_seed)
        reliability_bar(out_dir / "reliability_ranking.png", set_records, cfg.reliability_threshold)

    provenance = dict(georank3d_version=__version__, python=sys.version.split()[0],
                      platform=platform.platform(), numpy=np.__version__,
                      random_seed=cfg.random_seed, runtime_seconds=round(time.time() - t0, 2),
                      stage_counts=dict(sample_info,
                                        retained_after_planarity_filter=int(keep.sum()),
                                        sets_retained=len(set_records),
                                        patches_retained=len(patch_rows)))
    (out_dir / "run_provenance.json").write_text(json.dumps(provenance, indent=2))
    cfg.save(out_dir / "resolved_config.json")

    results = dict(config=asdict(cfg), provenance=provenance, sets=tabular_sets,
                   patches=patch_rows, wedge_mean=wedge_mean,
                   set_dispersion=set_dispersion, wedge_dispersion=wedge_disp,
                   uncertainty=mc, validation=validation, block_check=block_check)
    (out_dir / "results.json").write_text(json.dumps(results, indent=2, default=str))
    _log(cfg, f"GeoRank_3D complete in {provenance['runtime_seconds']} s -> {out_dir}")
    return results
