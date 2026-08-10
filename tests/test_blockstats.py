"""Block-inventory statistics tested against the values printed in the manuscript.

Every expected value below is quoted from the manuscript (Section 4.4,
Figs. 16 to 18) and is checked against this implementation running on the
Table 5 inventory in ``data/chilas_block_inventory.csv``.
"""

from pathlib import Path

import pytest

from georank3d.blockstats import analyse_blocks, read_block_inventory

ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "data" / "chilas_block_inventory.csv"


@pytest.fixture(scope="module")
def analysis():
    return analyse_blocks(read_block_inventory(INVENTORY), n_bootstrap=2000, seed=42)


def test_descriptive_statistics_match_section_4_4(analysis):
    d = analysis["descriptive"]
    assert d["count"] == 31
    assert d["sum_m3"] == pytest.approx(222.13, abs=0.01)
    assert d["mean_m3"] == pytest.approx(7.17, abs=0.01)
    assert d["median_m3"] == pytest.approx(2.10, abs=0.01)
    assert d["min_m3"] == pytest.approx(0.06, abs=0.01)
    assert d["max_m3"] == pytest.approx(40.10, abs=0.01)
    assert d["coefficient_of_variation"] == pytest.approx(1.49, abs=0.01)


def test_equivalent_diameter_formula(analysis):
    """D_eq = (6V/pi)^(1/3); Table 5 reports 4.25 m for the 40.10 m3 block."""
    biggest = max(analysis["inventory"], key=lambda r: r["volume_m3"])
    assert biggest["equivalent_diameter_m"] == pytest.approx(4.25, abs=0.01)


def test_lognormal_fit_matches_figure_16(analysis):
    ln = analysis["distributions"]["lognormal"]
    assert ln["ln_mean"] == pytest.approx(0.94, abs=0.01)
    assert ln["ln_std"] == pytest.approx(1.53, abs=0.01)
    assert ln["geometric_mean_m3"] == pytest.approx(2.56, abs=0.01)
    assert ln["ks_statistic"] == pytest.approx(0.120, abs=0.001)
    assert ln["ks_pvalue"] == pytest.approx(0.722, abs=0.002)
    assert ln["shapiro_pvalue"] == pytest.approx(0.637, abs=0.002)
    assert ln["qq_r_squared"] == pytest.approx(0.976, abs=0.001)


def test_information_criteria_select_the_lognormal(analysis):
    models = {m["model"]: m for m in analysis["distributions"]["models"]}
    assert analysis["distributions"]["best_by_aic"] == "lognormal"
    assert analysis["distributions"]["best_by_bic"] == "lognormal"
    assert models["lognormal"]["aic"] == pytest.approx(176.66, abs=0.05)
    assert models["lognormal"]["bic"] == pytest.approx(179.53, abs=0.05)
    for other in ("weibull", "gamma", "exponential", "normal"):
        assert models[other]["aic"] > models["lognormal"]["aic"]


def test_volume_concentration_matches_figure_17(analysis):
    by_threshold = {t["threshold_m3"]: t for t in analysis["magnitude_frequency"]["thresholds"]}
    assert by_threshold[20.0]["n_ge"] == 4
    assert by_threshold[20.0]["volume_share_percent"] == pytest.approx(56.65, abs=0.02)
    assert by_threshold[5.0]["volume_share_percent"] == pytest.approx(87.3, abs=0.05)
    assert by_threshold[1.0]["volume_share_percent"] == pytest.approx(98.2, abs=0.05)
    assert analysis["gini"] == pytest.approx(0.66, abs=0.01)


def test_sector_structure_matches_table_6(analysis):
    sectors = {s["sector"]: s for s in analysis["sectors"]["sectors"]}
    assert sectors["Middle"]["n_blocks"] == 17
    assert sectors["Middle"]["total_volume_m3"] == pytest.approx(165.68, abs=0.01)
    assert sectors["Middle"]["volume_share_percent"] == pytest.approx(74.59, abs=0.01)
    assert sectors["Middle"]["n_blocks_ge_20m3"] == 4
    assert sectors["Side"]["n_blocks"] == 8
    assert sectors["Side"]["volume_share_percent"] == pytest.approx(19.04, abs=0.01)
    assert sectors["Lower"]["n_blocks"] == 6
    assert sectors["Lower"]["volume_share_percent"] == pytest.approx(6.37, abs=0.01)
    assert analysis["sectors"]["kruskal_pvalue"] == pytest.approx(0.424, abs=0.002)


def test_bootstrap_interval_contains_the_point_estimate(analysis):
    for s in analysis["sectors"]["sectors"]:
        lo, hi = s["share_ci95_percent"]
        assert lo <= s["volume_share_percent"] <= hi


def test_rejects_a_non_positive_volume(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("block_id,sector,volume_m3\n1,Middle,0\n")
    with pytest.raises(ValueError, match="finite and positive"):
        read_block_inventory(p)


def test_rejects_a_file_without_a_volume_column(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("block_id,sector\n1,Middle\n")
    with pytest.raises(ValueError, match="no volume column"):
        read_block_inventory(p)
