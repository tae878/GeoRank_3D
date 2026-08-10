"""Geomechanical Reliability Index (manuscript Eq. 5).

    GRI_s = w_O O_s + w_P P_s + w_R R_s + w_A A_s + w_K K_s

with

    O_s = min(log10(N_s) / abundance_log_scale, 1)     point abundance
    P_s = mean per-point planarity of the set          Eq. (3)
    R_s = exp(-mean_roughness / roughness_scale)       roughness quality
    A_s = min(persistence / persistence_reference, 1)  persistence
    K_s = kinematic feasibility score                  Section 3.3.1

The weights are read from the configuration (``reliability_weights``) rather
than hard-coded, so the sensitivity analysis the manuscript describes can be
run without editing source.  Every term is returned alongside the score so
that a reported GRI can be audited component by component; on the Chilas
configuration the abundance and persistence terms both saturate at 1 for every
slope-scale set, which is why the reference release returned near-identical
scores (0.728-0.730) for all three sets.
"""

from __future__ import annotations

from typing import Dict

import numpy as np

__all__ = ["reliability_index"]


def reliability_index(n_points: int, mean_planarity: float, mean_roughness: float,
                      persistence_m: float, kinematic_score: float, cfg) -> Dict[str, float]:
    """Return the GRI and its five components for one discontinuity set."""
    w_o, w_p, w_r, w_a, w_k = cfg.reliability_weights
    o = min(np.log10(max(n_points, 1)) / cfg.abundance_log_scale, 1.0)
    p = float(np.clip(mean_planarity, 0.0, 1.0))
    r = float(np.exp(-max(mean_roughness, 0.0) / cfg.roughness_scale_m))
    a = min(max(persistence_m, 0.0) / cfg.persistence_reference_m, 1.0)
    k = float(np.clip(kinematic_score, 0.0, 1.0))
    gri = w_o * o + w_p * p + w_r * r + w_a * a + w_k * k
    return dict(gri=float(np.clip(gri, 0.0, 1.0)),
                term_abundance=float(o), term_planarity=p, term_roughness=r,
                term_persistence=float(a), term_kinematic=k,
                term_abundance_saturated=bool(o >= 1.0 - 1e-12),
                term_persistence_saturated=bool(a >= 1.0 - 1e-12))
