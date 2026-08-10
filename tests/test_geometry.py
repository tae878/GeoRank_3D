"""Geometry and orientation-convention tests."""

import numpy as np
import pytest

from georank3d.geometry import (angular_difference_deg, dip_from_normal, fit_plane,
                                line_of_intersection, mean_axial_direction, normal_from_dip,
                                projected_hull)


@pytest.mark.parametrize("dip,dipdir", [(0, 0), (10, 45), (45, 90), (72, 208), (80, 17),
                                        (73, 185), (89, 300), (30, 359)])
def test_normal_dip_roundtrip(dip, dipdir):
    d, a = dip_from_normal(normal_from_dip(dip, dipdir))
    assert d == pytest.approx(dip, abs=1e-8)
    if dip > 1e-9:
        assert (a - dipdir + 180) % 360 - 180 == pytest.approx(0.0, abs=1e-8)


def test_upward_normal_points_down_dip():
    """A surface descending eastward has an upward normal pointing east.

    z = -tan(delta) x  ->  tan(delta) x + z = 0  ->  n proportional to
    (tan delta, 0, 1), i.e. (sin delta, 0, cos delta).
    """
    n = np.array([np.tan(np.radians(30.0)), 0.0, 1.0])
    n /= np.linalg.norm(n)
    dip, dipdir = dip_from_normal(n)
    assert dip == pytest.approx(30.0, abs=1e-9)
    assert dipdir == pytest.approx(90.0, abs=1e-9)
    assert normal_from_dip(30.0, 90.0) == pytest.approx(n, abs=1e-12)


def test_normal_is_orthogonal_to_the_down_dip_vector():
    """The down-dip unit vector is (cos d sin a, cos d cos a, -sin d)."""
    for dip, dipdir in [(55, 200), (12, 4), (78, 346)]:
        d, a = np.radians(dip), np.radians(dipdir)
        downdip = np.array([np.cos(d) * np.sin(a), np.cos(d) * np.cos(a), -np.sin(d)])
        assert normal_from_dip(dip, dipdir) @ downdip == pytest.approx(0.0, abs=1e-12)


def test_manuscript_equation_9_is_180_degrees_out():
    """Regression guard documenting the defect in Eq. (9) and its MATLAB code.

    ``GeoRank_3D_orientation.m`` computes mod(atan2(-nx, -ny), 360) from a
    normal it has already forced into the upper hemisphere; the minus signs
    belong with a downward-oriented normal.
    """
    n = normal_from_dip(55.0, 200.0)
    as_published = np.degrees(np.arctan2(-n[0], -n[1])) % 360
    assert abs((as_published - 200.0) % 360 - 180.0) < 1e-9
    assert dip_from_normal(n)[1] == pytest.approx(200.0, abs=1e-9)


def test_dip_equals_arccos_abs_nz():
    """delta = arccos(|c|), manuscript Eq. (8)."""
    n = normal_from_dip(37.0, 123.0)
    assert dip_from_normal(n)[0] == pytest.approx(np.degrees(np.arccos(abs(n[2]))), abs=1e-9)


def test_orientation_is_axial():
    """n and -n describe the same plane and must give the same attitude."""
    n = normal_from_dip(78.0, 346.0)
    assert dip_from_normal(n) == pytest.approx(dip_from_normal(-n), abs=1e-9)


def test_lower_hemisphere_pole_trend_is_dipdir_plus_180():
    """Standard stereonet relation: pole plunges 90-delta towards alpha+180."""
    delta, alpha = 55.0, 200.0
    p, t = np.radians(90.0 - delta), np.radians((alpha + 180.0) % 360.0)
    pole_down = np.array([np.cos(p) * np.sin(t), np.cos(p) * np.cos(t), -np.sin(p)])
    assert -pole_down == pytest.approx(normal_from_dip(delta, alpha), abs=1e-12)


def test_angular_difference_is_antipodally_symmetric():
    a = normal_from_dip(80, 17)
    b = normal_from_dip(72, 208)
    assert angular_difference_deg(a, b)[0, 0] == pytest.approx(
        angular_difference_deg(a, -b)[0, 0], abs=1e-9)


def test_mean_axial_direction_ignores_arbitrary_signs():
    rng = np.random.default_rng(3)
    axis = normal_from_dip(80.0, 17.0)
    N = axis + 0.06 * rng.standard_normal((4000, 3))
    N /= np.linalg.norm(N, axis=1, keepdims=True)
    flipped = np.where(rng.random(len(N))[:, None] < 0.5, -N, N)
    m1, _ = mean_axial_direction(N)
    m2, _ = mean_axial_direction(flipped)
    assert angular_difference_deg(m1, m2)[0, 0] < 1e-6
    assert angular_difference_deg(m1, axis)[0, 0] < 1.0


def test_arithmetic_mean_of_flipped_normals_is_biased_near_the_seam():
    """The reference release averaged hemisphere-flipped normals; this is why."""
    rng = np.random.default_rng(11)
    axis = normal_from_dip(88.0, 17.0)          # pole plunges only 2 degrees
    N = axis + 0.10 * rng.standard_normal((20000, 3))
    N /= np.linalg.norm(N, axis=1, keepdims=True)
    flipped = np.where((N[:, 2] < 0)[:, None], -N, N)
    naive = flipped.mean(axis=0)
    naive /= np.linalg.norm(naive)
    tensor, _ = mean_axial_direction(N)
    assert angular_difference_deg(tensor, axis)[0, 0] < 0.5
    assert angular_difference_deg(naive, axis)[0, 0] > 1.0


def test_fit_plane_recovers_a_known_plane():
    rng = np.random.default_rng(5)
    n = normal_from_dip(55.0, 120.0)
    e1 = np.cross(n, [0, 0, 1.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    uv = rng.uniform(-5, 5, (2000, 2))
    P = uv[:, :1] * e1 + uv[:, 1:] * e2 + 0.001 * rng.standard_normal((2000, 1)) * n
    plane = fit_plane(P)
    assert plane["dip"] == pytest.approx(55.0, abs=0.2)
    assert plane["dipdir"] == pytest.approx(120.0, abs=0.5)
    assert plane["rmse"] < 0.01


def test_projected_hull_area_and_persistence_of_a_square():
    n = np.array([0.0, 0.0, 1.0])
    g = np.linspace(-2, 2, 40)
    X, Y = np.meshgrid(g, g)
    P = np.column_stack([X.ravel(), Y.ravel(), np.zeros(X.size)])
    area, persistence = projected_hull(P, n)
    assert area == pytest.approx(16.0, rel=1e-6)
    assert persistence == pytest.approx(4.0 * np.sqrt(2), rel=1e-6)


def test_line_of_intersection_of_two_vertical_planes_is_vertical():
    n1 = normal_from_dip(90.0, 0.0)
    n2 = normal_from_dip(90.0, 90.0)
    plunge, _ = line_of_intersection(n1, n2)
    assert plunge == pytest.approx(90.0, abs=1e-6)


def test_line_of_intersection_trend_and_plunge():
    n1 = normal_from_dip(60.0, 0.0)
    n2 = normal_from_dip(60.0, 90.0)
    plunge, trend = line_of_intersection(n1, n2)
    assert 0.0 < plunge < 60.0
    assert 0.0 <= trend < 360.0
    # the intersection of two planes dipping N and E must plunge to the NE quadrant
    assert 0.0 < trend < 90.0
