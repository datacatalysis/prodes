"""Draws the figures for the Prodes against APBS comparison.

Deliberately separate from the measurement, which lives in
scripts/apbs_comparison.py, so that a figure can be redrawn without a five
minute APBS run::

    source activate prodes
    python scripts/plot_apbs_comparison.py

The per protein figures are drawn from docs/visualisation/apbs_vs_prodes_summary.csv,
which is committed, so anyone can regenerate them without access to the machine
the numbers were measured on. The scatter of individual surface points needs the
per point table for the reference structure, which lives under the git ignored
apbs_results/ and is only present on a machine that has run the comparison; that
figure is skipped with a message when it is absent.

One plot per file, no multi-axis figures, per the DataCatalysis coding guide.
"""

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

# Selected before pyplot is imported, so that the script also runs on a headless server.
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
SUMMARY_FILE = REPO_ROOT / "docs" / "visualisation" / "apbs_vs_prodes_summary.csv"
POINTS_FILE = REPO_ROOT / "apbs_results" / "1BSQ" / "1BSQ_pdie4_points.csv.gz"
FIGURE_DIR = REPO_ROOT / "docs" / "visualisation" / "figures"

PRODES_COLOUR = "#c1272d"
APBS_COLOUR = "#0b6fa4"


def plot_fraction_positive(summary, out_file):
    """Fraction of the surface each method calls positive, one row per protein.

    The most legible statement of the result: Prodes finds far less positive
    surface than APBS on every protein, and none at all on four of them.
    """

    data = summary.sort_values("apbs_fraction_positive")
    positions = np.arange(len(data))

    figure, axes = plt.subplots(figsize=(7, 6))
    axes.hlines(positions, data.prodes_fraction_positive * 100, data.apbs_fraction_positive * 100, color="grey", linewidth=1, zorder=1)
    axes.scatter(data.prodes_fraction_positive * 100, positions, color=PRODES_COLOUR, s=45, zorder=2, label="Prodes")
    axes.scatter(data.apbs_fraction_positive * 100, positions, color=APBS_COLOUR, s=45, zorder=2, label="APBS")

    axes.set_yticks(positions)
    axes.set_yticklabels([f"{row.id}  {row.protein}" for row in data.itertuples()], fontsize=8)
    axes.set_xlabel("surface points with positive electrostatic potential (%)")
    axes.set_title("Prodes finds far less positive surface than APBS", fontsize=11)
    axes.legend(frameon=False, loc="lower right")
    axes.spines[["top", "right"]].set_visible(False)

    figure.tight_layout()
    figure.savefig(out_file, dpi=200)
    plt.close(figure)


def plot_spearman_against_recall(summary, out_file):
    """Rank correlation against recall of the APBS positive surface.

    The point of the figure: every protein sits in the "good agreement" band on
    the horizontal axis while four of them recover none of the positive surface
    at all. A rank correlation cannot see the failure, because it is invariant to
    the additive offset that causes it.
    """

    figure, axes = plt.subplots(figsize=(6.5, 5))
    axes.scatter(summary.spearman_P_A, summary.apbs_positive_recall * 100, color=APBS_COLOUR, s=45)

    for row in summary.itertuples():
        axes.annotate(row.id, (row.spearman_P_A, row.apbs_positive_recall * 100), fontsize=7, xytext=(4, 3), textcoords="offset points")

    axes.axhline(0, color="grey", linewidth=0.8, linestyle=":")
    axes.set_xlabel("Spearman rank correlation, Prodes against APBS")
    axes.set_ylabel("APBS positive surface recovered by Prodes (%)")
    axes.set_title("High rank correlation hides a total failure on positive patches", fontsize=11)
    axes.spines[["top", "right"]].set_visible(False)

    figure.tight_layout()
    figure.savefig(out_file, dpi=200)
    plt.close(figure)


def plot_decomposition(summary, out_file):
    """Where the agreement is lost: the charge model or the solvent physics.

    P is Prodes, C is the same Coulomb kernel on the force field charges, A is
    APBS. P to C isolates the charge model, C to A the solvent physics.
    """

    data = summary.sort_values("spearman_P_A")
    positions = np.arange(len(data))
    width = 0.4

    figure, axes = plt.subplots(figsize=(7, 6))
    axes.barh(positions + width / 2, data.spearman_P_C, height=width, color=PRODES_COLOUR, label="P to C, charge model")
    axes.barh(positions - width / 2, data.spearman_C_A, height=width, color=APBS_COLOUR, label="C to A, solvent physics")

    axes.set_yticks(positions)
    axes.set_yticklabels(data.id, fontsize=8)
    axes.set_xlim(0.5, 1.0)
    axes.set_xlabel("Spearman rank correlation")
    axes.set_title("Both the charge model and the solvent physics cost agreement", fontsize=11)
    axes.legend(frameon=False, loc="lower left")
    axes.spines[["top", "right"]].set_visible(False)

    figure.tight_layout()
    figure.savefig(out_file, dpi=200)
    plt.close(figure)


def plot_point_scatter(points, out_file, label):
    """Every surface point of one protein, Prodes potential against APBS potential.

    The quadrants are the story: the lower right is empty, so Prodes almost never
    claims a positive potential wrongly, while the upper left is full, which is
    the positive surface it misses.
    """

    figure, axes = plt.subplots(figsize=(6, 5.5))
    mesh = axes.hexbin(points.ep_prodes_volts, points.ep_apbs_volts, gridsize=70, bins="log", cmap="viridis", mincnt=1)

    axes.axhline(0, color="white", linewidth=0.8)
    axes.axvline(0, color="white", linewidth=0.8)
    axes.set_xlabel("Prodes electrostatic potential [V]")
    axes.set_ylabel("APBS electrostatic potential [V]")
    axes.set_title(f"{label}: agreement in rank, disagreement in sign", fontsize=11)

    bar = figure.colorbar(mesh, ax=axes)
    bar.set_label("surface points")

    figure.tight_layout()
    figure.savefig(out_file, dpi=200)
    plt.close(figure)


def main():
    """Draws every figure into docs/visualisation/figures/."""

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(SUMMARY_FILE)

    plot_fraction_positive(summary, FIGURE_DIR / "apbs_fraction_positive.png")
    plot_spearman_against_recall(summary, FIGURE_DIR / "apbs_spearman_vs_recall.png")
    plot_decomposition(summary, FIGURE_DIR / "apbs_decomposition.png")
    print(f"wrote three figures from {SUMMARY_FILE.name} into {FIGURE_DIR}")

    if POINTS_FILE.exists():
        plot_point_scatter(pd.read_csv(POINTS_FILE), FIGURE_DIR / "apbs_point_scatter_1BSQ.png", "beta-lactoglobulin (1BSQ)")
        print(f"wrote the point scatter from {POINTS_FILE.name}")
    else:
        print(f"skipping the point scatter: {POINTS_FILE} is not present, run scripts/apbs_comparison.py to produce it")


if __name__ == "__main__":
    main()
