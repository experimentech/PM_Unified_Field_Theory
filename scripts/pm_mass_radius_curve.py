#!/usr/bin/env python3
"""PM Compact-Star Mass–Radius Curve.

Generates M–R curves for PM compact stars using the self-consistent PM EOS
    P = c²(ρ − ρ_nuc)/2
and the PM structure equations (exact within PM; see stellar_structure.py).

Observational constraints overlaid:
  • PSR J0030+0451  (NICER, Riley et al. 2019):  M = 1.34 M☉, R = 12.71 km  (±)
  • PSR J0740+6620  (NICER, Riley et al. 2021):  M = 2.08 M☉, R = 12.35 km  (±)
  • PSR J0952-0607  (heaviest known, Romani et al. 2022): M = 2.35 ± 0.17 M☉
  • GW170817 constraint:  R_1.4 ≈ 11.9 km  (Abbott et al. 2018)

Usage
-----
    python scripts/pm_mass_radius_curve.py               # interactive plot
    python scripts/pm_mass_radius_curve.py --save        # save to results/
    python scripts/pm_mass_radius_curve.py --no-plot     # print table only
"""

import argparse
import sys
from pathlib import Path

import numpy as np

# Allow running from the project root
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.pushing_medium.stellar_structure import (
    RHO_CRIT,
    compute_mr_curve,
    solve_pm_star,
)
from src.pushing_medium.critical_state import RHO_NUC


# ---------------------------------------------------------------------------
# Observational boxes (M ± δM, R ± δR in solar masses and km)
# ---------------------------------------------------------------------------

OBSERVATIONS = [
    {
        "label": "J0030+0451\n(NICER)",
        "M":  1.34,  "dM": 0.15,
        "R":  12.71, "dR": 1.14,
        "color": "steelblue",
        "marker": "D",
    },
    {
        "label": "J0740+6620\n(NICER)",
        "M":  2.08,  "dM": 0.07,
        "R":  12.35, "dR": 0.75,
        "color": "darkorange",
        "marker": "s",
    },
    {
        "label": "J0952-0607\n(heaviest)",
        "M":  2.35,  "dM": 0.17,
        "R":  None,  "dR": None,   # radius unknown
        "color": "crimson",
        "marker": "^",
    },
    {
        "label": "GW170817\nR$_{1.4}$",
        "M":  1.4,   "dM": 0.0,    # fixed horizontal line
        "R":  11.9,  "dR": 1.4,    # estimated spread
        "color": "purple",
        "marker": "x",
    },
]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(n_points: int = 80, save: bool = False, no_plot: bool = False) -> None:
    print("PM Compact-Star Mass–Radius Curve")
    print("=" * 42)
    print(f"  ρ_nuc  = {RHO_NUC:.3e} kg/m³")
    print(f"  ρ_crit = {RHO_CRIT:.3e} kg/m³  (= e × ρ_nuc)")
    print()

    # ------------------------------------------------------------------
    # Compute M–R curve
    # ------------------------------------------------------------------
    print(f"Integrating {n_points} stellar models …", end="", flush=True)
    rho_c_arr, M_arr, R_arr = compute_mr_curve(n_points=n_points)
    valid = ~np.isnan(M_arr)
    print(f"  done. ({valid.sum()}/{n_points} converged)")
    print()

    # ------------------------------------------------------------------
    # Find the maximum-mass star
    # ------------------------------------------------------------------
    if valid.any():
        i_max = np.nanargmax(M_arr)
        M_max = M_arr[i_max]
        R_max_star = R_arr[i_max]
        rho_max = rho_c_arr[i_max]
        print(f"PM maximum-mass star (at ρ_c = ρ_crit = {RHO_CRIT:.3e} kg/m³):")
        print(f"  M_max = {M_max:.3f} M_☉")
        print(f"  R     = {R_max_star:.2f} km")
        print()

    # ------------------------------------------------------------------
    # Print summary table
    # ------------------------------------------------------------------
    print(f"{'ρ_c [kg/m³]':>18}  {'M [M_☉]':>10}  {'R [km]':>10}")
    print("-" * 44)
    step = max(1, n_points // 20)
    for i in range(0, n_points, step):
        if valid[i]:
            print(f"  {rho_c_arr[i]:.4e}    {M_arr[i]:>8.3f}    {R_arr[i]:>8.2f}")
    print()

    if no_plot:
        return

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg" if save else "TkAgg")
        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available — skipping plot.")
        return

    fig, ax = plt.subplots(figsize=(8, 6))

    # PM curve
    ax.plot(
        R_arr[valid], M_arr[valid],
        color="tab:green", linewidth=2.5, label="PM (exact)",
        zorder=5,
    )

    # Mark the PM maximum-mass point
    if valid.any():
        ax.scatter(
            [R_max_star], [M_max],
            color="tab:green", s=100, zorder=6,
            marker="*", linewidths=1.5,
            label=f"PM max: {M_max:.2f} M$_\\odot$, {R_max_star:.1f} km",
        )

    # Mark ρ_crit contour annotation
    ax.annotate(
        r"$\rho_c = \rho_\mathrm{crit} = e\,\rho_\mathrm{nuc}$",
        xy=(R_max_star, M_max),
        xytext=(R_max_star + 1.5, M_max - 0.15),
        fontsize=9,
        arrowprops=dict(arrowstyle="-|>", color="tab:green", lw=1.0),
        color="tab:green",
    )

    # Observational constraints
    for obs in OBSERVATIONS:
        if obs["R"] is not None:
            rect = mpatches.FancyBboxPatch(
                (obs["R"] - obs["dR"], obs["M"] - obs["dM"]),
                2 * obs["dR"], 2 * obs["dM"],
                boxstyle="round,pad=0.05",
                linewidth=1.2,
                edgecolor=obs["color"],
                facecolor=obs["color"],
                alpha=0.18,
                zorder=3,
            )
            ax.add_patch(rect)
            ax.scatter(
                [obs["R"]], [obs["M"]],
                color=obs["color"], marker=obs["marker"],
                s=60, zorder=4, label=obs["label"].replace("\n", " "),
            )
        else:
            # J0952-0607: only mass known — draw horizontal band
            ax.axhspan(
                obs["M"] - obs["dM"], obs["M"] + obs["dM"],
                alpha=0.12, color=obs["color"],
                label=obs["label"].replace("\n", " ") + f" M = {obs['M']:.2f} M$_\\odot$",
            )

    ax.set_xlabel("Radius  [km]", fontsize=13)
    ax.set_ylabel(r"Mass  [$M_\odot$]", fontsize=13)
    ax.set_title("PM Compact-Star Mass–Radius Curve\n"
                 r"EOS: $P = c^2(\rho - \rho_\mathrm{nuc})/2$,   "
                 r"$\rho_c \leq e\,\rho_\mathrm{nuc}$",
                 fontsize=11)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save:
        out_path = Path(__file__).parent.parent / "results" / "pm_mr_curve.pdf"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved to {out_path}")
        out_png = out_path.with_suffix(".png")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        print(f"Saved to {out_png}")
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="PM compact-star mass–radius curve"
    )
    parser.add_argument("--n-points", type=int, default=80,
                        help="Number of stellar models to integrate (default 80)")
    parser.add_argument("--save", action="store_true",
                        help="Save plot to results/ instead of displaying")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip plotting; print table only")
    args = parser.parse_args()

    run(n_points=args.n_points, save=args.save, no_plot=args.no_plot)
