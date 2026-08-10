"""Clustering, kinematics, reliability and end-to-end pipeline tests."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from georank3d.clustering import cluster_normals, spatial_patches
from georank3d.config import Config, ConfigError
from georank3d.geometry import angular_difference_deg, dip_from_normal, normal_from_dip
from georank3d.kinematics import kinematic_score, planar_toppling_flags
from georank3d.reliability import reliability_index

ROOT = Path(__file__).resolve().parents[1]


def fisher(axis, kappa, n, rng):
    axis = axis / np.linalg.norm(axis)
    w = 1 + np.log(rng.random(n) + (1 - rng.random(n)) * np.exp(-2 * kappa)) / kappa
    th = 2 * np.pi * rng.random(n)
    s = np.sqrt(np.clip(1 - w ** 2, 0, None))
    v = np.stack([s * np.cos(th), s * np.sin(th), w], axis=1)
    z = np.array([0.0, 0.0, 1.0])
    ax = np.cross(z, axis)
    if np.linalg.norm(ax) < 1e-12:
        return v if axis[2] > 0 else -v
    ax /= np.linalg.norm(ax)
    ang = np.arccos(np.clip(z @ axis, -1, 1))
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    R = np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K
    return v @ R.T


# ------------------------------------------------------------------ config
def test_weights_must_sum_to_one():
    with pytest.raises(ConfigError, match="sum to 1"):
        Config(reliability_weights=[0.2, 0.2, 0.2, 0.2, 0.1]).validate()


def test_weights_must_have_five_entries():
    with pytest.raises(ConfigError, match="five entries"):
        Config(reliability_weights=[0.5, 0.5]).validate()


def test_unknown_key_is_rejected():
    with pytest.raises(ConfigError, match="unknown configuration key"):
        Config.from_dict({"slope_dip": 80})


@pytest.mark.parametrize("kwargs,pattern", [
    (dict(normal_k=2), "normal_k"),
    (dict(planarity_min=1.5), "planarity_min"),
    (dict(friction_angle_deg=0), "friction_angle_deg"),
    (dict(persistence_source="whatever"), "persistence_source"),
    (dict(merge_deg=40, set_cone_deg=20), "merge_deg"),
])
def test_invalid_configurations_are_rejected(kwargs, pattern):
    with pytest.raises(ConfigError, match=pattern):
        Config(**kwargs).validate()


def test_config_roundtrip(tmp_path):
    cfg = Config(slope_dipdir_deg=140.0, monte_carlo_runs=3).validate()
    p = tmp_path / "c.json"
    cfg.save(p)
    assert Config.load(p).slope_dipdir_deg == 140.0


# -------------------------------------------------------------- clustering
def test_clustering_recovers_three_known_sets():
    rng = np.random.default_rng(0)
    truth = [(82, 11), (67, 29), (78, 346)]
    N = np.vstack([fisher(normal_from_dip(d, a), 200, 4000, rng) for d, a in truth])
    cfg = Config(min_cluster_points=500, set_cone_deg=18.0, density_bandwidth_deg=8.0).validate()
    labels, poles = cluster_normals(N, cfg, np.random.default_rng(1))
    assert len(poles) == 3
    for d, a in truth:
        devs = [angular_difference_deg(normal_from_dip(d, a), p)[0, 0] for p in poles]
        assert min(devs) < 2.0


def test_clustering_does_not_truncate_a_near_vertical_set():
    """A pole 5 degrees from the n_z = 0 seam must keep both of its flanks.

    The reference release compared hemisphere-flipped normals by Euclidean
    chord, which discards the members that dip past vertical.
    """
    rng = np.random.default_rng(2)
    axis = normal_from_dip(85.0, 17.0)
    N = fisher(axis, 60, 20000, rng)
    cfg = Config(min_cluster_points=500, set_cone_deg=22.0).validate()
    labels, _ = cluster_normals(N, cfg, np.random.default_rng(3))
    within = angular_difference_deg(N, axis)[:, 0] <= cfg.set_cone_deg
    assert (labels > 0).sum() >= 0.97 * within.sum()


def test_clustering_is_deterministic_for_a_fixed_seed():
    rng = np.random.default_rng(4)
    N = np.vstack([fisher(normal_from_dip(d, a), 150, 3000, rng) for d, a in [(80, 17), (60, 200)]])
    cfg = Config(min_cluster_points=400).validate()
    a, _ = cluster_normals(N, cfg, np.random.default_rng(9))
    b, _ = cluster_normals(N, cfg, np.random.default_rng(9))
    assert np.array_equal(a, b)


def test_spatial_patches_separate_disconnected_surfaces():
    rng = np.random.default_rng(6)
    left = rng.uniform([-10, -3, 0], [-4, 3, 5], (1500, 3))
    right = rng.uniform([20, -3, 0], [26, 3, 5], (1500, 3))
    cfg = Config(patch_radius_m=1.0, min_patch_points=10).validate()
    labels = spatial_patches(np.vstack([left, right]), cfg, np.random.default_rng(7))
    assert len(np.unique(labels[labels > 0])) == 2


# -------------------------------------------------------------- kinematics
def test_planar_criterion_requires_daylighting_and_friction():
    cfg = Config(slope_dip_deg=80.0, slope_dipdir_deg=135.0, friction_angle_deg=35.0).validate()
    planar, _, _, _ = planar_toppling_flags([60.0], [135.0], cfg)
    assert bool(planar[0])
    assert not bool(planar_toppling_flags([30.0], [135.0], cfg)[0][0])   # below friction
    assert not bool(planar_toppling_flags([85.0], [135.0], cfg)[0][0])   # steeper than face
    assert not bool(planar_toppling_flags([60.0], [200.0], cfg)[0][0])   # outside window


def test_toppling_criterion_uses_the_anti_slope_direction():
    cfg = Config(slope_dip_deg=80.0, slope_dipdir_deg=135.0, friction_angle_deg=35.0).validate()
    _, toppling, _, _ = planar_toppling_flags([80.0], [315.0], cfg)
    assert bool(toppling[0])
    assert not bool(planar_toppling_flags([80.0], [135.0], cfg)[1][0])


def test_kinematic_score_weights_and_cap():
    cfg = Config().validate()
    assert kinematic_score([True], [False], cfg)[0] == pytest.approx(0.65)
    assert kinematic_score([False], [True], cfg)[0] == pytest.approx(0.45)
    assert kinematic_score([True], [True], cfg)[0] == pytest.approx(1.0)
    assert kinematic_score([False], [False], cfg)[0] == pytest.approx(0.0)


# ------------------------------------------------------------- reliability
def test_reliability_uses_configured_weights():
    """Regression guard: the reference release ignored opts.reliabilityWeights."""
    base = Config().validate()
    other = Config(reliability_weights=[0.0, 1.0, 0.0, 0.0, 0.0]).validate()
    args = dict(n_points=10_000, mean_planarity=0.6, mean_roughness=0.02,
                persistence_m=4.0, kinematic_score=0.65)
    assert reliability_index(**args, cfg=other)["gri"] == pytest.approx(0.6)
    assert reliability_index(**args, cfg=base)["gri"] != pytest.approx(0.6)


def test_reliability_terms_saturate_as_documented():
    cfg = Config().validate()
    r = reliability_index(1_000_000, 0.8, 0.0, 100.0, 0.0, cfg)
    assert r["term_abundance_saturated"] and r["term_persistence_saturated"]
    assert r["gri"] == pytest.approx(0.20 + 0.25 * 0.8 + 0.15 * 1.0 + 0.20 * 1.0)


def test_reliability_components_reconstruct_the_score():
    cfg = Config().validate()
    r = reliability_index(5000, 0.42, 0.05, 3.1, 0.45, cfg)
    rebuilt = sum(w * r[t] for w, t in zip(cfg.reliability_weights,
                                           ["term_abundance", "term_planarity", "term_roughness",
                                            "term_persistence", "term_kinematic"]))
    assert r["gri"] == pytest.approx(rebuilt)


# -------------------------------------------------------------- end to end
@pytest.mark.slow
def test_end_to_end_synthetic(tmp_path):
    cloud = tmp_path / "cloud.xyz"
    subprocess.run([sys.executable, str(ROOT / "scripts" / "make_synthetic_slope.py"),
                    "--output", str(cloud), "--points-per-set", "4500",
                    "--slope-points", "9000", "--noise-points", "600"],
                   check=True, cwd=ROOT, capture_output=True)
    out = tmp_path / "results"
    cfg = dict(input_cloud=str(cloud), output_dir=str(out), voxel_size_m=0.0,
               normal_k=25, min_cluster_points=300, set_cone_deg=18.0,
               patch_radius_m=1.2, min_patch_points=40, monte_carlo_runs=3,
               monte_carlo_subsample=6000, make_figures=False, verbose=False)
    cfg_path = tmp_path / "cfg.json"
    cfg_path.write_text(json.dumps(cfg))
    proc = subprocess.run([sys.executable, str(ROOT / "scripts" / "run_pipeline.py"),
                           "--config", str(cfg_path)], cwd=ROOT, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr

    results = json.loads((out / "results.json").read_text())
    found = [(r["dip_deg"], r["dipdir_deg"]) for r in results["sets"]]
    for d, a in [(82, 11), (67, 29), (78, 346), (80, 135)]:
        ref = normal_from_dip(d, a)
        assert min(angular_difference_deg(ref, normal_from_dip(*f))[0, 0] for f in found) < 3.0

    for name in ("set_metrics.csv", "patch_metrics.csv", "kinematic_dispersion_sets.csv",
                 "kinematic_dispersion_wedges.csv", "uncertainty_sets.csv",
                 "resolved_config.json", "run_provenance.json"):
        assert (out / name).is_file(), name


def test_missing_input_is_reported_clearly(tmp_path):
    cfg_path = tmp_path / "cfg.json"
    cfg_path.write_text(json.dumps(dict(input_cloud=str(tmp_path / "nope.xyz"),
                                        output_dir=str(tmp_path / "o"))))
    proc = subprocess.run([sys.executable, str(ROOT / "scripts" / "run_pipeline.py"),
                           "--config", str(cfg_path)], cwd=ROOT, capture_output=True, text=True)
    assert proc.returncode == 1
    assert "point cloud not found" in proc.stderr


def test_dip_direction_convention_is_shared_by_the_dispersion_module():
    """The population screening must use manuscript Eq. (9), like everything else."""
    from georank3d.kinematics import dispersion_screening
    cfg = Config(slope_dip_deg=80.0, slope_dipdir_deg=135.0, friction_angle_deg=35.0).validate()
    rng = np.random.default_rng(21)
    N = fisher(normal_from_dip(60.0, 135.0), 400, 4000, rng)   # squarely planar-feasible
    rec = dict(set_id=1, planar_feasible_mean=1, toppling_feasible_mean=0)
    rows = dispersion_screening({1: N}, [rec], cfg)
    assert rows[0]["planar_feasible_fraction"] > 0.9
    assert rows[0]["toppling_feasible_fraction"] < 0.05
    assert dip_from_normal(N)[1].mean() == pytest.approx(135.0, abs=2.0)
