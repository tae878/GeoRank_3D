"""Validation of extracted sets against independent reference sources.

Two checks are implemented, both described in manuscript Section 3.3.4:

*   angular deviation of each GeoRank_3D set from the nearest set reported by
    manual mapping, qFacet Fast Marching, qFacet KD-tree or DSE, against the
    10-degree acceptance threshold of Buyer and Schubert (2017); and
*   a first-order block-volume consistency check, V = c * A^p (Eq. 12),
    compared with a field block census by a two-sample Kolmogorov-Smirnov test
    on the log volumes.

The matching is reported in both directions.  Reporting only
"GeoRank set -> nearest reference set" hides reference sets that no extracted
set reproduces, which is the quantity a reader needs in order to judge whether
the workflow recovered the structural model.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np

from .geometry import normal_from_dip

__all__ = ["read_reference_sets", "compare_to_references", "block_volume_comparison"]

_DIP_KEYS = ("dip", "meandip", "mean_dip", "dip_deg", "dipangle")
_DIR_KEYS = ("dipdir", "dipdirection", "mean_dipdir", "meandipdir", "dip_dir",
             "dipdir_deg", "dipdirection_deg", "azimuth")
_LABEL_KEYS = ("label", "name", "setid", "set_id", "clusterid", "cluster_id", "id", "set")


def _norm_key(s: str) -> str:
    return "".join(ch for ch in s.strip().lower() if ch.isalnum() or ch == "_")


def read_reference_sets(path) -> List[Dict]:
    """Read a reference-set CSV with dip / dip-direction columns.

    Column names are matched case-insensitively against a list of accepted
    aliases; the dip-direction alias list is tested before the dip list so that
    a column literally named ``DipDirection`` is never consumed as the dip.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"reference-set CSV not found: {p}")
    with p.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"{p} contains no data rows")
    headers = list(rows[0].keys())
    keys = {h: _norm_key(h) for h in headers}

    dir_col = next((h for h in headers if keys[h] in _DIR_KEYS), None)
    dip_col = next((h for h in headers if keys[h] in _DIP_KEYS and h != dir_col), None)
    if dip_col is None or dir_col is None:
        raise ValueError(
            f"{p}: could not identify dip and dip-direction columns.\n"
            f"  columns present : {headers}\n"
            f"  accepted dip    : {_DIP_KEYS}\n"
            f"  accepted dipdir : {_DIR_KEYS}")
    label_col = next((h for h in headers if keys[h] in _LABEL_KEYS), None)

    out = []
    for i, row in enumerate(rows):
        try:
            dip = float(row[dip_col])
            dipdir = float(row[dir_col])
        except (TypeError, ValueError):
            continue
        if not (np.isfinite(dip) and np.isfinite(dipdir)):
            continue
        if not 0.0 <= dip <= 90.0:
            raise ValueError(f"{p} row {i + 2}: dip {dip} outside [0, 90]; "
                             "check that the dip and dip-direction columns are not transposed")
        out.append(dict(label=(row[label_col] if label_col else f"Set{i + 1}"),
                        dip_deg=dip, dipdir_deg=dipdir % 360.0))
    if not out:
        raise ValueError(f"{p}: no valid dip / dip-direction pairs parsed")
    return out


def _axial_deviation(a, b) -> float:
    return float(np.degrees(np.arccos(min(abs(float(np.dot(a, b))), 1.0))))


def compare_to_references(set_records: List[Dict], sources: Dict[str, List[Dict]],
                          threshold_deg: float) -> Dict[str, List[Dict]]:
    """Two-way angular comparison of extracted sets and reference sources."""
    extracted = [(rec["set_id"], normal_from_dip(rec["dip_deg"], rec["dipdir_deg"])) for rec in set_records]

    forward = []
    for rec, (sid, ne) in zip(set_records, extracted):
        row = dict(set_id=sid, dip_deg=rec["dip_deg"], dipdir_deg=rec["dipdir_deg"])
        for name, refs in sources.items():
            devs = [(_axial_deviation(ne, normal_from_dip(r["dip_deg"], r["dipdir_deg"])), r["label"])
                    for r in refs]
            best = min(devs) if devs else (float("nan"), "")
            row[f"dev_{name}_deg"] = best[0]
            row[f"match_{name}"] = best[1]
            row[f"pass_{name}"] = bool(np.isfinite(best[0]) and best[0] <= threshold_deg)
        forward.append(row)

    reverse = []
    for name, refs in sources.items():
        for r in refs:
            nr = normal_from_dip(r["dip_deg"], r["dipdir_deg"])
            devs = [(_axial_deviation(nr, ne), sid) for sid, ne in extracted]
            best = min(devs) if devs else (float("nan"), None)
            reverse.append(dict(source=name, reference_label=r["label"],
                                reference_dip_deg=r["dip_deg"], reference_dipdir_deg=r["dipdir_deg"],
                                nearest_set_id=best[1], deviation_deg=best[0],
                                recovered=bool(np.isfinite(best[0]) and best[0] <= threshold_deg)))
    return dict(extracted_to_reference=forward, reference_to_extracted=reverse)


def block_volume_comparison(patch_areas: Sequence[float], observed_volumes: Sequence[float],
                            cfg) -> Dict:
    """First-order block-volume proxy versus a field block census (Eq. 12).

    Patches larger than ``block_ceiling_multiple`` times the largest observed
    block are reported separately as exposure-scale surfaces rather than
    detachable blocks.  Because that ceiling is derived from the validation
    data, the comparison is repeated over a range of multiples so a reader can
    see how sensitive the agreement is to the cut-off.
    """
    from scipy import stats

    areas = np.asarray([a for a in patch_areas if np.isfinite(a) and a > 0], dtype=float)
    obs = np.asarray([v for v in observed_volumes if np.isfinite(v) and v > 0], dtype=float)
    if areas.size == 0 or obs.size == 0:
        raise ValueError("block-volume comparison needs at least one patch area and one observed volume")

    proxy_all = cfg.block_volume_coefficient * areas ** cfg.block_volume_exponent
    result = dict(n_patches_total=int(proxy_all.size),
                  observed=dict(count=int(obs.size), median=float(np.median(obs)),
                                mean=float(obs.mean()), max=float(obs.max()),
                                std=float(obs.std(ddof=1)) if obs.size > 1 else float("nan")),
                  sensitivity=[])

    for multiple in sorted({1.0, 2.0, cfg.block_ceiling_multiple, 10.0, float("inf")}):
        ceiling = multiple * obs.max()
        kept = proxy_all[proxy_all <= ceiling]
        entry = dict(ceiling_multiple=(None if not np.isfinite(multiple) else multiple),
                     ceiling_m3=(None if not np.isfinite(ceiling) else float(ceiling)),
                     n_kept=int(kept.size), n_excluded=int(proxy_all.size - kept.size))
        if kept.size >= 2:
            ks = stats.ks_2samp(np.log10(kept), np.log10(obs))
            entry.update(median_m3=float(np.median(kept)), mean_m3=float(kept.mean()),
                         max_m3=float(kept.max()), std_m3=float(kept.std(ddof=1)),
                         ks_statistic=float(ks.statistic), ks_pvalue=float(ks.pvalue),
                         median_ratio=float(np.median(kept) / np.median(obs)))
        result["sensitivity"].append(entry)

    primary = next(e for e in result["sensitivity"]
                   if e["ceiling_multiple"] == cfg.block_ceiling_multiple)
    result["primary"] = primary
    return result
