"""Run configuration for GeoRank_3D.

Every parameter that affects a result is declared here, is serialisable to
JSON and is written verbatim into the output directory alongside the results.
No parameter is hard-coded anywhere else in the package; in particular the
five Geomechanical Reliability Index weights of manuscript Eq. (5) are read
from :attr:`Config.reliability_weights` (in the reference MATLAB release the
option existed but was never read, the weights being literals inside
``GeoRank_3D_reliability_score.m``).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict, fields
from pathlib import Path
from typing import Any, Dict, List


class ConfigError(ValueError):
    """Raised when a configuration value is missing or out of range."""


@dataclass
class Config:
    # ------------------------------------------------------------ I/O
    input_cloud: str = ""
    output_dir: str = "georank3d_output"
    random_seed: int = 42
    export_classified_cloud: bool = True

    # ------------------------------------------------- preprocessing
    voxel_size_m: float = 0.05          # 0 disables voxel filtering
    max_points: int = 10_000_000        # analysis budget
    uniform_reduce_to_budget: bool = True

    # ------------------------------------- normals and point filter
    normal_k: int = 30                  # k nearest neighbours, Eq. (2)
    planarity_min: float = 0.20         # Eq. (3), P_i
    surface_variation_max: float = 0.18  # Eq. (3), S_i

    # ------------------------------------------- orientation sets
    density_bandwidth_deg: float = 10.0   # mean-shift kernel half-width
    set_cone_deg: float = 22.0            # max axial pole-to-normal angle
    merge_deg: float = 8.0                # merge poles closer than this
    max_sets: int = 12
    min_cluster_points: int = 300
    pole_prominence_frac: float = 0.03
    cluster_sample_cap: int = 200_000

    # -------------------------------------------- spatial patches
    patch_radius_m: float = 0.75
    min_patch_points: int = 80
    patch_sample_cap: int = 150_000

    # ------------------------------------------------ slope and phi
    slope_dip_deg: float = 80.0
    slope_dipdir_deg: float = 135.0
    friction_angle_deg: float = 35.0
    planar_window_deg: float = 25.0       # planar-sliding azimuth window
    toppling_window_deg: float = 30.0     # direct-toppling azimuth window
    wedge_window_deg: float = 35.0        # wedge trend window
    planar_kin_weight: float = 0.65       # contribution to K_s
    toppling_kin_weight: float = 0.45     # contribution to K_s

    # ------------------------------------------------------ GRI, Eq. (5)
    # order: [abundance O, planarity P, roughness R, persistence A, kinematic K]
    reliability_weights: List[float] = field(default_factory=lambda: [0.20, 0.25, 0.15, 0.20, 0.20])
    reliability_threshold: float = 0.50
    abundance_log_scale: float = 5.0      # O_s = min(log10(N)/scale, 1)
    roughness_scale_m: float = 0.10       # R_s = exp(-r/scale)
    persistence_reference_m: float = 10.0  # A_s = min(L/reference, 1)
    persistence_source: str = "max_patch"  # "max_patch" or "set_hull"

    # ---------------------------------------- dispersion screening
    dispersion_sample_cap: int = 20_000

    # ------------------------------------------- Monte Carlo (3.3.2)
    monte_carlo_runs: int = 50
    monte_carlo_subsample: int = 20_000
    mc_rotation_sigma_deg: float = 2.0
    mc_angular_tolerance_sigma_deg: float = 2.0
    mc_min_cluster_jitter: float = 0.25    # +/- fraction
    mc_assignment_threshold_deg: float = 30.0

    # ------------------------------------------------- validation
    reference_sets: Dict[str, str] = field(default_factory=dict)
    match_threshold_deg: float = 10.0      # Buyer and Schubert (2017)
    observed_blocks_csv: str = ""
    block_volume_coefficient: float = 0.20   # V = c * A^exponent, Eq. (12)
    block_volume_exponent: float = 1.5
    block_ceiling_multiple: float = 5.0      # exposure-surface cut-off

    # ---------------------------------------------------- plotting
    max_plot_poles: int = 150_000
    make_figures: bool = True

    # -------------------------------------------------------- misc
    verbose: bool = True

    # ------------------------------------------------------------------
    def validate(self) -> "Config":
        """Check ranges and internal consistency; raise :class:`ConfigError`."""
        def _pos(name, value, allow_zero=False):
            if value < 0 or (value == 0 and not allow_zero):
                raise ConfigError(f"{name} must be > 0 (got {value!r})")

        _pos("normal_k", self.normal_k)
        if self.normal_k < 4:
            raise ConfigError("normal_k must be at least 4 for a stable 3x3 covariance")
        _pos("voxel_size_m", self.voxel_size_m, allow_zero=True)
        _pos("max_points", self.max_points)
        for name in ("planarity_min", "surface_variation_max"):
            v = getattr(self, name)
            if not 0.0 <= v <= 1.0:
                raise ConfigError(f"{name} must lie in [0, 1] (got {v!r})")
        for name in ("density_bandwidth_deg", "set_cone_deg", "merge_deg"):
            v = getattr(self, name)
            if not 0.0 < v < 90.0:
                raise ConfigError(f"{name} must lie in (0, 90) degrees (got {v!r})")
        if self.merge_deg > self.set_cone_deg:
            raise ConfigError("merge_deg must not exceed set_cone_deg")
        if not 0.0 <= self.slope_dip_deg <= 90.0:
            raise ConfigError("slope_dip_deg must lie in [0, 90]")
        if not 0.0 <= self.slope_dipdir_deg < 360.0:
            raise ConfigError("slope_dipdir_deg must lie in [0, 360)")
        if not 0.0 < self.friction_angle_deg < 90.0:
            raise ConfigError("friction_angle_deg must lie in (0, 90)")
        w = list(self.reliability_weights)
        if len(w) != 5:
            raise ConfigError("reliability_weights must have five entries "
                              "[abundance, planarity, roughness, persistence, kinematic]")
        if any(x < 0 for x in w):
            raise ConfigError("reliability_weights must be non-negative")
        if abs(sum(w) - 1.0) > 1e-9:
            raise ConfigError(f"reliability_weights must sum to 1.0 (got {sum(w)!r}); "
                              "manuscript Eq. (5) defines a normalised index")
        if self.persistence_source not in ("max_patch", "set_hull"):
            raise ConfigError("persistence_source must be 'max_patch' or 'set_hull'")
        if not 0.0 <= self.reliability_threshold <= 1.0:
            raise ConfigError("reliability_threshold must lie in [0, 1]")
        _pos("persistence_reference_m", self.persistence_reference_m)
        _pos("roughness_scale_m", self.roughness_scale_m)
        _pos("patch_radius_m", self.patch_radius_m)
        if self.monte_carlo_runs < 0:
            raise ConfigError("monte_carlo_runs must be >= 0")
        return self

    # ------------------------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def save(self, path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True))

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Config":
        known = {f.name for f in fields(cls)}
        unknown = set(data) - known
        if unknown:
            raise ConfigError(f"unknown configuration key(s): {sorted(unknown)}")
        return cls(**data).validate()

    @classmethod
    def load(cls, path) -> "Config":
        p = Path(path)
        if not p.is_file():
            raise ConfigError(f"configuration file not found: {p}")
        try:
            data = json.loads(p.read_text())
        except json.JSONDecodeError as exc:
            raise ConfigError(f"{p} is not valid JSON: {exc}") from exc
        return cls.from_dict(data)
