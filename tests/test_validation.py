"""Validation-layer tests, including the manuscript's Table 8 deviations."""

from pathlib import Path

import pytest

from georank3d.config import Config
from georank3d.validation import (block_volume_comparison, compare_to_references,
                                  read_reference_sets)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

# GeoRank_3D sets as reported in manuscript Table 8
GEORANK_SETS = [
    dict(set_id=3, dip_deg=80.0, dipdir_deg=17.0),
    dict(set_id=2, dip_deg=73.0, dipdir_deg=185.0),
    dict(set_id=1, dip_deg=72.0, dipdir_deg=208.0),
]

# Table 8 "Dev." columns, degrees
TABLE_8 = {3: dict(manual=6.5, qfacet_fm=9.8, qfacet_kd=8.6, dse=12.3),
           2: dict(manual=25.7, qfacet_fm=31.0, qfacet_kd=36.4, dse=9.0),
           1: dict(manual=31.3, qfacet_fm=39.8, qfacet_kd=38.4, dse=20.9)}


@pytest.fixture(scope="module")
def sources():
    return {name: read_reference_sets(DATA / f"reference_sets_{name}.csv")
            for name in ("manual", "qfacet_fm", "qfacet_kd", "dse")}


def test_reference_csvs_parse(sources):
    assert len(sources["manual"]) == 3
    assert len(sources["dse"]) == 10
    assert sources["manual"][0]["dip_deg"] == 82.0
    assert sources["manual"][0]["dipdir_deg"] == 11.0


def test_table_8_deviations_are_reproduced(sources):
    """Table 8 is internally consistent with arccos(|n_a . n_b|) on Table 4,
    provided both sides are read as quoted attitudes in the same convention.

    The manuscript rounds Table 4 to whole degrees, so a tolerance of 0.5 deg
    is the most that transcription allows.  See
    :func:`test_reported_attitudes_are_not_pca_normals` for why this is not the
    same comparison the MATLAB pipeline actually performs.
    """
    result = compare_to_references(GEORANK_SETS, sources, 10.0)
    for row in result["extracted_to_reference"]:
        expected = TABLE_8[row["set_id"]]
        for name, value in expected.items():
            assert row[f"dev_{name}_deg"] == pytest.approx(value, abs=0.5), \
                f"set {row['set_id']} vs {name}"


def test_only_set_3_passes_the_10_degree_threshold_against_manual(sources):
    result = compare_to_references(GEORANK_SETS, sources, 10.0)
    passes = {r["set_id"]: r["pass_manual"] for r in result["extracted_to_reference"]}
    assert passes == {3: True, 2: False, 1: False}


def test_reverse_matching_exposes_unrecovered_reference_sets(sources):
    """Manual JS-2 and JS-3 have no GeoRank counterpart within 10 degrees.

    Reporting only the forward direction, as Table 8 does, hides this.
    """
    result = compare_to_references(GEORANK_SETS, sources, 10.0)
    manual = {r["reference_label"]: r for r in result["reference_to_extracted"]
              if r["source"] == "manual"}
    assert manual["JS-1"]["recovered"] is True
    assert manual["JS-2"]["recovered"] is False
    assert manual["JS-3"]["recovered"] is False


def test_dse_major_sets_are_not_within_5_degrees_of_manual_mapping(sources):
    """Section 4.3.1 states the DSE major sets deviate by < 5 deg from manual.

    Computed from Table 4 with the manuscript's own metric they do not.
    """
    manual = {r["label"]: r for r in sources["manual"]}
    dse = {r["label"]: r for r in sources["dse"]}
    result = compare_to_references(
        [dict(set_id=i, dip_deg=dse[f"JS-{i}"]["dip_deg"], dipdir_deg=dse[f"JS-{i}"]["dipdir_deg"])
         for i in (1, 2, 3)],
        {"manual": list(manual.values())}, 5.0)
    devs = {r["set_id"]: r["dev_manual_deg"] for r in result["extracted_to_reference"]}
    assert min(devs.values()) > 5.0
    assert devs[2] > 30.0            # DSE JS-2 (40/002) vs manual JS-2 (67/029)
    assert devs[3] > 15.0            # DSE JS-3 (81/189) vs any manual set


def test_reported_attitudes_are_not_pca_normals(sources):
    """Quantifies the consequence of the Eq. (9) sign error for Table 8.

    ``GeoRank_3D_validate_multi_source.m`` converts the extracted set through
    ``dipdir2normal(orientation(n_pca))``, which returns ``n_pca`` unchanged,
    but converts each reference attitude through the same flipped formula, so
    the reference enters rotated by 180 degrees.  Only one of the two vectors
    is rotated, so the tabulated deviation is not the angle between the planes.
    """
    from georank3d.geometry import angular_difference_deg, normal_from_dip

    manual_js1 = normal_from_dip(82.0, 11.0)
    # what Table 8 reports for Set 3 vs manual: both sides as quoted attitudes
    assert angular_difference_deg(normal_from_dip(80.0, 17.0), manual_js1)[0, 0] \
        == pytest.approx(6.5, abs=0.5)
    # the actual plane behind "Set 3 = 80/017" dips towards 197, not 017
    assert angular_difference_deg(normal_from_dip(80.0, 197.0), manual_js1)[0, 0] \
        == pytest.approx(19.0, abs=0.5)


def test_corrected_convention_reassigns_the_field_match(sources):
    """Under a consistent convention the field match moves from Set 3 to Set 1.

    Set 1 is the population the manuscript describes as a slope-parallel
    exposure surface with "no field counterpart closer than 20.9 deg"; on the
    corrected convention it is the closest match to manual JS-2, at 5.1 deg,
    while Set 3 no longer reproduces any manually mapped set within 10 deg.
    """
    from georank3d.geometry import angular_difference_deg, normal_from_dip

    manual = {r["label"]: normal_from_dip(r["dip_deg"], r["dipdir_deg"]) for r in sources["manual"]}
    true_set1 = normal_from_dip(72.0, 28.0)      # reported 72/208
    true_set3 = normal_from_dip(80.0, 197.0)     # reported 80/017

    best_1 = min((angular_difference_deg(true_set1, n)[0, 0], k) for k, n in manual.items())
    best_3 = min((angular_difference_deg(true_set3, n)[0, 0], k) for k, n in manual.items())
    assert best_1[1] == "JS-2" and best_1[0] == pytest.approx(5.1, abs=0.5)
    assert best_3[0] > 10.0


def test_dip_direction_column_is_not_consumed_as_dip(tmp_path):
    p = tmp_path / "r.csv"
    p.write_text("Label,DipDirection,Dip\nA,185,73\n")
    sets = read_reference_sets(p)
    assert sets[0]["dip_deg"] == 73.0
    assert sets[0]["dipdir_deg"] == 185.0


def test_transposed_columns_raise_rather_than_silently_mislead(tmp_path):
    p = tmp_path / "r.csv"
    p.write_text("Label,Dip,DipDir\nA,185,73\n")
    with pytest.raises(ValueError, match="outside"):
        read_reference_sets(p)


def test_block_volume_proxy_reports_ceiling_sensitivity():
    cfg = Config().validate()
    observed = [0.5, 2.1, 7.2, 40.1]
    areas = [1.0, 4.0, 9.0, 25.0, 400.0]          # the last is an exposure-scale face
    out = block_volume_comparison(areas, observed, cfg)
    assert out["primary"]["ceiling_multiple"] == 5.0
    assert out["primary"]["n_excluded"] == 1       # 0.2 * 400^1.5 = 1600 m3
    assert 0.0 <= out["primary"]["ks_pvalue"] <= 1.0
    multiples = [e["ceiling_multiple"] for e in out["sensitivity"]]
    assert None in multiples and 1.0 in multiples and 10.0 in multiples
