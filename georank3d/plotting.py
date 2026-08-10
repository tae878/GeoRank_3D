"""Figures: Schmidt pole plot, reliability ranking and block statistics.

The pole plot follows the Lambert (Schmidt) equal-area lower-hemisphere
convention.  The pole of a plane dipping ``delta`` towards ``alpha`` plunges
``90 - delta`` towards ``alpha + 180``, so it plots at radius
``sqrt(2) sin(delta/2)`` in the direction ``alpha + 180``.

Note for anyone patching the MATLAB release: ``GeoRank_3D_plot_poles.m`` plots
at azimuth ``alpha`` rather than ``alpha + 180``, but it takes ``alpha`` from
``GeoRank_3D_orientation``, which is itself 180 degrees out.  The two errors
cancel, so the published stereonet is geometrically correct while the
tabulated attitudes are not.  Correcting the orientation convention alone
would silently break the figure; both must be fixed together.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .geometry import dip_from_normal  # noqa: E402

__all__ = ["pole_plot", "reliability_bar", "block_statistics_figures", "schmidt_xy"]


def schmidt_xy(dip_deg, dipdir_deg):
    """Lower-hemisphere equal-area coordinates of the pole to a plane."""
    dip = np.asarray(dip_deg, dtype=float)
    trend = (np.asarray(dipdir_deg, dtype=float) + 180.0) % 360.0
    r = np.sqrt(2.0) * np.sin(np.radians(dip) / 2.0)
    return r * np.sin(np.radians(trend)), r * np.cos(np.radians(trend))


def _net(ax):
    th = np.linspace(0, 2 * np.pi, 361)
    ax.plot(np.sin(th), np.cos(th), "k-", lw=1.3)
    for d in (30, 60):
        rr = np.sqrt(2) * np.sin(np.radians(d / 2))
        ax.plot(rr * np.sin(th), rr * np.cos(th), "--", color="0.72", lw=0.5)
        ax.text(rr * 0.72, rr * 0.72, str(d), fontsize=7, color="0.6", ha="center")
    for label, dx, dy in (("N", 0, 1), ("E", 1, 0), ("S", 0, -1), ("W", -1, 0)):
        ax.plot([0.95 * dx, 1.06 * dx], [0.95 * dy, 1.06 * dy], "k-", lw=1.2)
        ax.text(1.17 * dx, 1.17 * dy, label, ha="center", va="center",
                fontsize=11, fontweight="bold")
    ax.set_aspect("equal")
    ax.axis("off")


def pole_plot(out_path, normals: np.ndarray, labels: np.ndarray, set_records: List[Dict],
              max_poles: int = 150_000, seed: int = 42):
    rng = np.random.default_rng(seed)
    if len(normals) > max_poles:
        sel = rng.choice(len(normals), max_poles, replace=False)
        normals, labels = normals[sel], labels[sel]
    dip, dipdir = dip_from_normal(normals)
    x, y = schmidt_xy(dip, dipdir)

    fig, ax = plt.subplots(figsize=(6.4, 6.6), dpi=150)
    unassigned = labels == 0
    if unassigned.any():
        ax.scatter(x[unassigned], y[unassigned], s=3, c="0.75", alpha=0.4,
                   linewidths=0, label="unassigned")
    for rec in set_records:
        m = labels == rec["set_id"]
        if m.any():
            ax.scatter(x[m], y[m], s=3, alpha=0.55, linewidths=0,
                       label=f"Set {rec['set_id']}  {rec['dip_deg']:.0f}/{rec['dipdir_deg']:03.0f}")
    for rec in set_records:
        mx, my = schmidt_xy(rec["dip_deg"], rec["dipdir_deg"])
        ax.plot(mx, my, "k+", ms=13, mew=2)
        ax.text(mx + 0.05, my + 0.05, f"S{rec['set_id']}", fontsize=9, fontweight="bold")
    _net(ax)
    ax.legend(loc="upper left", bbox_to_anchor=(-0.08, 1.02), fontsize=7, frameon=False,
              markerscale=3)
    ax.set_title("GeoRank_3D: Schmidt equal-area, lower-hemisphere pole plot", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def reliability_bar(out_path, set_records: List[Dict], threshold: float):
    fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=150)
    ids = [f"S{r['set_id']}" for r in set_records]
    terms = ["term_abundance", "term_planarity", "term_roughness",
             "term_persistence", "term_kinematic"]
    names = ["abundance $O_s$", "planarity $P_s$", "roughness $R_s$",
             "persistence $A_s$", "kinematic $K_s$"]
    weights = set_records[0]["reliability_weights"] if set_records else [0.2] * 5
    bottom = np.zeros(len(set_records))
    for t, label, w in zip(terms, names, weights):
        vals = np.array([w * r[t] for r in set_records])
        ax.bar(ids, vals, 0.6, bottom=bottom, label=f"{label} (w={w:g})")
        bottom += vals
    ax.axhline(threshold, ls="--", color="r", lw=1.2)
    ax.text(0.02, threshold + 0.01, f"acceptance threshold = {threshold:g}",
            transform=ax.get_yaxis_transform(), color="r", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Discontinuity set")
    ax.set_ylabel("Geomechanical Reliability Index (stacked contributions)")
    ax.set_title("GeoRank_3D: reliability ranking, decomposed by Eq. (5) term", fontsize=10)
    ax.legend(fontsize=7, ncol=2, frameon=False)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def block_statistics_figures(out_dir, analysis: Dict) -> List[str]:
    """Figures 16, 17 and 18 of the manuscript."""
    from scipy import stats as sps

    out_dir = Path(out_dir)
    records = analysis["inventory"]
    V = np.array([r["volume_m3"] for r in records])
    sectors = np.array([r["sector"] for r in records])
    lnV = np.log(V)
    mu = analysis["distributions"]["lognormal"]["ln_mean"]
    sd = analysis["distributions"]["lognormal"]["ln_std"]
    written = []

    # ---- Figure 16: distribution fitting ----------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), dpi=150)
    ax = axes[0, 0]
    ax.hist(lnV, bins=8, density=True, color="0.7", edgecolor="k", lw=0.5)
    xs = np.linspace(lnV.min() - 1, lnV.max() + 1, 400)
    ax.plot(xs, sps.norm.pdf(xs, mu, sd), "r-", lw=1.8,
            label=f"log-normal fit\n$\\mu$={mu:.2f}, $\\sigma$={sd:.2f}")
    ax.set_xlabel("ln V  (V in m$^3$)")
    ax.set_ylabel("Probability density")
    ax.set_title("(a) Frequency distribution of ln V", fontsize=10)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    xs = np.sort(lnV)
    ax.step(xs, np.arange(1, len(xs) + 1) / len(xs), where="post", label="empirical CDF")
    grid = np.linspace(xs.min() - 1, xs.max() + 1, 400)
    ax.plot(grid, sps.norm.cdf(grid, mu, sd), "r-", lw=1.5, label="fitted log-normal")
    ln = analysis["distributions"]["lognormal"]
    ax.set_title(f"(b) CDF, KS D = {ln['ks_statistic']:.3f}, p = {ln['ks_pvalue']:.3f}", fontsize=10)
    ax.set_xlabel("ln V")
    ax.set_ylabel("Cumulative probability")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    osm, osr = sps.probplot(lnV, dist="norm")[0]
    ax.plot(osm, osr, "o", ms=4)
    ax.plot(osm, mu + sd * osm, "r--", lw=1, label="fitted log-normal")
    ax.legend(fontsize=8)
    ax.set_xlabel("Theoretical normal quantile")
    ax.set_ylabel("Observed ln V quantile")
    ax.set_title(f"(c) Normal Q-Q of ln V, $R^2$ = {ln['qq_r_squared']:.3f}, "
                 f"Shapiro-Wilk p = {ln['shapiro_pvalue']:.3f}", fontsize=9)

    ax = axes[1, 1]
    models = analysis["distributions"]["models"]
    labels = [m["model"] for m in models]
    xpos = np.arange(len(models))
    ax.bar(xpos - 0.2, [m["aic"] for m in models], 0.4, label="AIC")
    ax.bar(xpos + 0.2, [m["bic"] for m in models], 0.4, label="BIC")
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("Information criterion")
    ax.set_title("(d) Model comparison (lower is better)", fontsize=10)
    ax.legend(fontsize=8)
    fig.suptitle("Fig. 16  Statistical characterisation of the detachable block volumes", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "fig16_block_volume_distribution.png"
    fig.savefig(p)
    plt.close(fig)
    written.append(str(p))

    # ---- Figure 17: magnitude-frequency and volume concentration ----------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), dpi=150)
    Vs = np.sort(V)[::-1]
    ax = axes[0]
    ax.loglog(Vs, np.arange(1, len(Vs) + 1), "o", ms=4, label="mapped blocks")
    grid = np.logspace(np.log10(V.min()), np.log10(V.max()), 200)
    ax.loglog(grid, len(V) * (1.0 - sps.norm.cdf(np.log(grid), mu, sd)), "r-", lw=1.5,
              label="fitted log-normal")
    ax.set_xlabel("Block volume V (m$^3$)")
    ax.set_ylabel("N($\\geq$V)")
    ax.set_title("(a) Magnitude-frequency relation", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    share = np.cumsum(Vs) / Vs.sum() * 100.0
    ax.semilogx(Vs, share, "o-", ms=4)
    for t in analysis["magnitude_frequency"]["thresholds"]:
        ax.axvline(t["threshold_m3"], color="0.6", ls=":", lw=1)
        ax.annotate(f"$\\geq${t['threshold_m3']:g} m$^3$: {t['volume_share_percent']:.1f}%",
                    xy=(t["threshold_m3"], t["volume_share_percent"]),
                    xytext=(4, -12), textcoords="offset points", fontsize=8)
    ax.set_xlabel("Block volume V (m$^3$)")
    ax.set_ylabel("Share of total mapped volume (%)")
    ax.set_title(f"(b) Volume-weighted exceedance, Gini = {analysis['gini']:.2f}", fontsize=10)
    ax.grid(alpha=0.3, which="both")
    fig.suptitle("Fig. 17  Magnitude-frequency and volume concentration", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = out_dir / "fig17_magnitude_frequency.png"
    fig.savefig(p)
    plt.close(fig)
    written.append(str(p))

    # ---- Figure 18: sector structure --------------------------------------
    order = [s["sector"] for s in analysis["sectors"]["sectors"]]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), dpi=150)
    ax = axes[0, 0]
    data = [V[sectors == s] for s in order]
    ax.boxplot(data, tick_labels=order, widths=0.5)
    for i, d in enumerate(data, start=1):
        ax.plot(np.full(len(d), i) + np.random.default_rng(1).normal(0, 0.04, len(d)), d,
                "o", ms=4, alpha=0.6)
    ax.set_yscale("log")
    ax.set_ylabel("Block volume (m$^3$)")
    kw = analysis["sectors"]["kruskal_pvalue"]
    ax.set_title(f"(a) Volume by source sector (Kruskal-Wallis p = {kw:.3f})", fontsize=10)

    ax = axes[0, 1]
    xpos = np.arange(len(order))
    shares = [s["volume_share_percent"] for s in analysis["sectors"]["sectors"]]
    los = [s["volume_share_percent"] - s["share_ci95_percent"][0] for s in analysis["sectors"]["sectors"]]
    his = [s["share_ci95_percent"][1] - s["volume_share_percent"] for s in analysis["sectors"]["sectors"]]
    ax.bar(xpos, shares, 0.55, yerr=[los, his], capsize=4, color="#2f6fa8")
    for i, s in enumerate(analysis["sectors"]["sectors"]):
        ax.text(i, shares[i] + 3, f"n = {s['n_blocks']}", ha="center", fontsize=9)
    ax.set_xticks(xpos)
    ax.set_xticklabels(order)
    ax.set_ylabel("Share of total detachable volume (%)")
    ax.set_title("(b) Volume share with bootstrap 95% CI", fontsize=10)

    ax = axes[1, 0]
    classes = ["small", "medium", "large", "extra-large"]
    grid_counts = np.array([[sum(1 for r in records if r["sector"] == s and r["volume_class"] == c)
                             for c in classes] for s in order], dtype=float)
    im = ax.imshow(grid_counts, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(classes, rotation=20, ha="right")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    for i in range(grid_counts.shape[0]):
        for j in range(grid_counts.shape[1]):
            ax.text(j, i, int(grid_counts[i, j]), ha="center", va="center", fontsize=9)
    fig.colorbar(im, ax=ax, label="block count")
    ax.set_title("(c) Sector by volume class", fontsize=10)

    ax = axes[1, 1]
    idx = np.argsort(V)[::-1]
    ax.bar(range(len(V)), V[idx], color="#2f6fa8")
    ax2 = ax.twinx()
    ax2.plot(range(len(V)), np.cumsum(V[idx]) / V.sum() * 100, "r-o", ms=3)
    ax2.set_ylabel("Cumulative share (%)", color="r")
    ax.set_xlabel("Blocks ranked by volume")
    ax.set_ylabel("Volume (m$^3$)")
    ax.set_title("(d) Pareto ranking of block volumes", fontsize=10)
    fig.suptitle("Fig. 18  Sector-wise structure of the detachable block-volume hazard", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / "fig18_sector_structure.png"
    fig.savefig(p)
    plt.close(fig)
    written.append(str(p))
    return written
