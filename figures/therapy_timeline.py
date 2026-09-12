"""
Therapy timeline and treatment response in refractory PR3-ANCA-associated
vasculitis treated with the BCMA-directed bispecific T-cell engager teclistamab.

Journal-style multi-panel figure (NEJM / Nature / Annals of the Rheumatic
Diseases conventions):

    panel 0  treatment timeline, one lane per therapy
    panel a  oral prednisolone maintenance dose
    panel b  anti-PR3 IgG (log scale)
    panel c  CD19+ B cells

Day 0 = first teclistamab dose.

Usage
-----
    python therapy_timeline.py                    # embedded, de-identified data
    python therapy_timeline.py --xlsx Therapie.xlsx   # re-read source workbook
    python therapy_timeline.py --outdir ../out --formats pdf png tiff
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import numpy as np
from matplotlib.patches import Rectangle

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# ---------------------------------------------------------------------------
# EDITABLE PARAMETERS  --  everything a co-author is likely to change
# ---------------------------------------------------------------------------

# Assay units. Verify against your laboratory report before submission.
UNIT_PR3 = "U/mL"
UNIT_CD19 = "cells/µL"

# Upper limit of normal of the anti-PR3 assay (set to None to hide the band).
PR3_ULN = 2.0

# Lower plotting floor for the log axis; values of 0 / below assay cut-off are
# drawn here as open symbols.
PR3_FLOOR = 0.6

# Guideline target for glucocorticoid maintenance (EULAR/PEXIVAS-style
# low-dose regimen). Set to None to hide.
PRED_TARGET = 5.0

# Imaging timepoints (days relative to first teclistamab dose).
# >>> REPLACE WITH THE EXACT SCAN DATES. <<<
# The boxes in the original figure are schematic rather than to scale, so
# these were placed on the clinical narrative instead: A at first
# presentation, B at re-induction, C at post-teclistamab follow-up.
SCANS = [(-318, "A"), (-97, "B"), (157, "C")]

# Teclistamab: 6 doses, 14 Jan 2026 through 19 Feb 2026 (day 0 to day 36).
# Only the exposure window is known from the source record, so the panel
# shows the window. Fill in TEC_DOSE_DAYS to draw the individual doses.
TEC_WINDOW = (0, 36)
TEC_DOSE_DAYS: list[int] = []

XLIM = (-340, 205)

# ---------------------------------------------------------------------------
# DATA  --  day offsets relative to first teclistamab dose (de-identified)
# ---------------------------------------------------------------------------

# (day, prednisolone mg/day)
PREDNISOLONE = [
    (-318, 500), (-315, 60), (-312, 100), (-303, 60), (-296, 100),
    (-291, 500), (-274, 40), (-254, 25), (-202, 7.5), (-170, 5),
    (-128, 0), (-97, 500), (-93, 250), (-82, 60), (-70, 40),
    (-30, 30), (0, 20), (12, 17.5), (58, 10), (83, 5),
    (111, 5), (156, 4), (188, 2),
]

# Days on which the recorded steroid dose is intravenous pulse / infusion
# premedication rather than the oral maintenance dose. These are lifted out of
# panel (a) and shown in the timeline instead, so that the taper stays legible.
# Each entry coincides with an infusion day in the source record.
PULSE_DAYS = {-318, -312, -296, -291, -97, -93}

# (day, anti-PR3 IgG)
PR3 = [
    (-318, 153), (-315, 28), (-303, 17), (-291, 20), (-254, 4.6),
    (-239, 3.5), (-232, 3.1), (-202, 2.9), (-170, 2.3), (-128, 4.9),
    (-97, 9.4), (-82, 5.2), (-55, 1.8), (-30, 1.5), (-12, 1.1),
    (12, 0), (58, 0), (83, 0), (111, 0), (156, 0), (188, 0),
]

# (day, CD19+ cells)
CD19 = [
    (-303, 1), (-291, 0), (-254, 0), (-202, 0), (-128, 22),
    (-82, 0), (-55, 0), (-12, 0), (12, 0), (58, 0),
    (111, 0), (156, 1), (188, 1),
]

# ---------------------------------------------------------------------------
# THERAPY LANES  --  drawn top to bottom in the timeline panel
# ---------------------------------------------------------------------------

C_PLEX = "#8E6FA8"   # plasma exchange   -- muted violet
C_GC = "#6F7B8A"     # glucocorticoid    -- slate
C_RTX = "#2E7D5B"    # rituximab         -- deep green
C_CYC = "#C08A21"    # cyclophosphamide  -- amber
C_AVA = "#9AA1AA"    # avacopan          -- grey
C_DARA = "#2E6E9E"   # daratumumab       -- blue
C_TEC = "#B5322A"    # teclistamab       -- red (the index therapy)

C_PRED = "#2B4C7E"   # prednisolone series -- blue, as in the original figure
C_PR3 = "#A8412C"
C_CD19 = "#35707E"

INK = "#1A1A1A"
MUTED = "#6B6B6B"
GRID = "#E4E4E4"

LANES = [
    dict(label="Plasma exchange", color=C_PLEX, kind="span",
         spans=[(-318, -315)], note="7 sessions"),
    dict(label="Glucocorticoid pulse", color=C_GC,
         doses=[-318, -312, -296, -291, -97, -93],
         note="0.1–0.5 g i.v., 6 pulses"),
    dict(label="Rituximab", color=C_RTX,
         doses=[-312, -296, -106, 173],
         annots=[(-304, "2 × 1 g"), (-106, "500 mg"), (173, "500 mg")]),
    dict(label="Cyclophosphamide", color=C_CYC,
         doses=[-315, -291, -97, -30, -12],
         annots=[(-303, "0.9 + 1.0 g"), (-97, "1.4 g"), (-21, "2 × 1.4 g")],
         note="6.1 g cumulative"),
    dict(label="Avacopan", color=C_AVA,
         spans=[(-291, 167)], note="30 mg twice daily"),
    dict(label="Daratumumab", color=C_DARA,
         doses=[-93, -82, -76, -70], note="4 × 1800 mg s.c."),
    dict(label="Teclistamab", color=C_TEC,
         spans=[TEC_WINDOW], doses=list(TEC_DOSE_DAYS),
         note="6 doses s.c. over 5 weeks"),
]


# ---------------------------------------------------------------------------
# STYLE
# ---------------------------------------------------------------------------

def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "Helvetica Neue",
                            "Liberation Sans", "DejaVu Sans"],
        "font.size": 7.0,
        "axes.labelsize": 7.5,
        "axes.titlesize": 7.5,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "legend.fontsize": 6.5,
        "axes.linewidth": 0.6,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.6,
        "ytick.major.size": 2.6,
        "xtick.minor.size": 1.5,
        "xtick.minor.width": 0.5,
        "lines.linewidth": 1.1,
        "lines.solid_capstyle": "round",
        "legend.frameon": False,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,   # embed as TrueType -- required by most journals
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

def load_from_xlsx(path: Path):
    """Re-read the source workbook; returns the same triple of series."""
    import openpyxl

    ws = openpyxl.load_workbook(path, data_only=True).worksheets[0]
    rows = list(ws.iter_rows(min_row=2, values_only=True))
    pred, pr3, cd19 = [], [], []
    for r in rows:
        day = r[1]
        if day is None:
            continue
        day = int(day)
        if r[3] is not None:
            pred.append((day, float(r[3])))
        if r[4] is not None:
            pr3.append((day, float(r[4])))
        if r[5] is not None:
            cd19.append((day, float(r[5])))
    return pred, pr3, cd19


def panel_letter(ax, letter: str, dx: float = -0.062, dy: float = 1.0) -> None:
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="left", color=INK)


def tidy(ax, last: bool = False) -> None:
    """NEJM-style axes: no top/right spine, offset left spine, light y grid."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_position(("outward", 3))
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, lw=0.5)
    ax.xaxis.grid(False)
    ax.set_xlim(*XLIM)
    if last:
        ax.spines["bottom"].set_position(("outward", 3))
        ax.set_xticks(np.arange(-300, 201, 100))
        ax.set_xticks(np.arange(-325, 201, 25), minor=True)
        ax.set_xlabel("Days from first teclistamab dose")
    else:
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", which="both", length=0, labelbottom=False)


def mark_index_therapy(ax, band: bool = True) -> None:
    """Red day-0 rule plus the teclistamab exposure band, on every panel."""
    if band:
        ax.axvspan(TEC_WINDOW[0], TEC_WINDOW[1], color=C_TEC, alpha=0.07,
                   lw=0, zorder=0)
    ax.axvline(0, color=C_TEC, lw=0.7, ls=(0, (3.5, 2.5)), zorder=1, alpha=0.85)


# ---------------------------------------------------------------------------
# PANELS
# ---------------------------------------------------------------------------

def draw_timeline(ax) -> None:
    n = len(LANES)
    ax.set_ylim(n - 0.45, -1.25)          # inverted: lane 0 on top
    ax.set_xlim(*XLIM)
    ax.axis("off")
    mark_index_therapy(ax)

    bar_h = 0.30

    for i, lane in enumerate(LANES):
        col = lane["color"]

        # faint baseline rule so each lane reads as its own track
        ax.plot(XLIM, [i, i], color=GRID, lw=0.5, zorder=0,
                solid_capstyle="butt")

        # continuous exposure
        for x0, x1 in lane.get("spans", []):
            width = max(x1 - x0, 9.0)     # keep very short spans visible
            ax.add_patch(Rectangle(
                (x0, i - bar_h / 2), width, bar_h,
                linewidth=0, facecolor=col, alpha=0.9, zorder=3))

        # discrete doses -- lollipops
        for day in sorted(lane.get("doses", [])):
            ax.plot([day, day], [i - 0.26, i + 0.26], color=col, lw=1.3,
                    solid_capstyle="round", zorder=4)
            ax.plot([day], [i - 0.26], marker="o", ms=2.9, mfc=col,
                    mec="white", mew=0.5, zorder=5)

        # dose annotations, placed explicitly so closely spaced infusions are
        # summarised rather than overprinted
        for day, text in lane.get("annots", []):
            ax.annotate(text, (day, i - 0.38), ha="center", va="bottom",
                        fontsize=5.9, color=col, fontweight="bold", zorder=6)

        # lane label, flush right of the axis
        note = lane.get("note")
        label = lane["label"]
        ax.annotate(label, (1.008, i), xycoords=("axes fraction", "data"),
                    ha="left", va="center" if not note else "bottom",
                    fontsize=6.8, color=INK, fontweight="bold",
                    annotation_clip=False)
        if note:
            ax.annotate(note, (1.008, i + 0.06),
                        xycoords=("axes fraction", "data"),
                        ha="left", va="top", fontsize=5.9, color=MUTED,
                        annotation_clip=False)

    # ---- imaging timepoints, on their own rail above the lanes ----
    ax.annotate("Imaging", (XLIM[0] - 6, -1.0), ha="right", va="center",
                fontsize=6.5, color=MUTED, annotation_clip=False)
    for day, letter in SCANS:
        ax.plot([day, day], [-1.0, n - 0.6], color=MUTED, lw=0.5,
                ls=(0, (1.2, 2.0)), zorder=1, alpha=0.55)
        ax.plot([day], [-1.0], marker="o", ms=8.0, mfc="white", mec=MUTED,
                mew=0.7, zorder=6, clip_on=False)
        ax.annotate(letter, (day, -1.0), ha="center", va="center",
                    fontsize=6.2, fontweight="bold", color=INK, zorder=7,
                    annotation_clip=False)


def draw_prednisolone(ax) -> None:
    maint = [(d, v) for d, v in PREDNISOLONE if d not in PULSE_DAYS]
    x = np.array([d for d, _ in maint], float)
    y = np.array([v for _, v in maint], float)

    # carry the last recorded dose to the right edge of the plot
    xs = np.append(x, XLIM[1])
    ys = np.append(y, y[-1])

    mark_index_therapy(ax)
    if PRED_TARGET is not None:
        ax.axhline(PRED_TARGET, color=MUTED, lw=0.6, ls=(0, (2.5, 2.0)),
                   zorder=1)
        ax.annotate(f"≤{PRED_TARGET:g} mg/day target", (XLIM[0] + 6,
                    PRED_TARGET + 1.6), ha="left", va="bottom", fontsize=5.9,
                    color=MUTED)

    ax.fill_between(xs, ys, step="post", color=C_PRED, alpha=0.13, lw=0,
                    zorder=2)
    ax.step(xs, ys, where="post", color=C_PRED, lw=1.2, zorder=3)
    ax.plot(x, y, ls="none", marker="o", ms=3.0, mfc="white", mec=C_PRED,
            mew=1.0, zorder=4, clip_on=False)

    ax.set_ylim(0, 66)
    ax.set_yticks([0, 20, 40, 60])
    ax.set_ylabel("Prednisolone\n(mg/day, oral)")
    tidy(ax)

    # annotate the two clinically meaningful excursions
    ax.annotate("relapse:\ndose re-escalated", (-82, 58),
                xytext=(-168, 52), fontsize=6.0, color=MUTED, va="center",
                ha="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=MUTED,
                                shrinkA=3, shrinkB=3))
    ax.annotate("2 mg/day", (186, 9), fontsize=6.4, color=C_PRED,
                fontweight="bold", ha="center", va="bottom")


def draw_pr3(ax) -> None:
    x = np.array([d for d, _ in PR3], float)
    y = np.array([v for _, v in PR3], float)
    below = y <= 0
    yplot = np.where(below, PR3_FLOOR, y)

    mark_index_therapy(ax)
    if PR3_ULN is not None:
        ax.axhspan(PR3_FLOOR * 0.5, PR3_ULN, color=GRID, alpha=0.75, lw=0,
                   zorder=0)
        ax.annotate("negative range", (XLIM[0] + 6, PR3_ULN * 0.80),
                    ha="left", va="top", fontsize=5.9, color=MUTED, zorder=5)

    ax.plot(x, yplot, color=C_PR3, lw=1.2, zorder=3)
    ax.plot(x[~below], yplot[~below], ls="none", marker="o", ms=3.0,
            mfc=C_PR3, mec="white", mew=0.6, zorder=4)
    ax.plot(x[below], yplot[below], ls="none", marker="o", ms=3.0,
            mfc="white", mec=C_PR3, mew=1.0, zorder=4)

    ax.set_yscale("log")
    ax.set_ylim(PR3_FLOOR * 0.5, 320)
    ax.set_yticks([1, 10, 100])
    ax.set_yticklabels(["1", "10", "100"])
    ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, subs=np.arange(2, 10)))
    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.tick_params(axis="y", which="minor", length=1.3, width=0.5)
    ax.set_ylabel(f"Anti-PR3 IgG\n({UNIT_PR3})")
    tidy(ax)

    ax.annotate("open symbols: below assay cut-off",
                (XLIM[0] + 6, PR3_FLOOR * 0.60), ha="left", va="bottom",
                fontsize=5.9, color=MUTED)
    ax.annotate("153", (-318, 153), xytext=(-310, 190), fontsize=6.4,
                color=C_PR3, fontweight="bold", ha="left", va="bottom")
    ax.annotate("serological remission", (105, 2.4), fontsize=6.4,
                color=C_PR3, ha="center", va="bottom", fontweight="bold")


def draw_cd19(ax) -> None:
    x = np.array([d for d, _ in CD19], float)
    y = np.array([v for _, v in CD19], float)

    mark_index_therapy(ax)
    ax.fill_between(x, y, color=C_CD19, alpha=0.13, lw=0, zorder=2)
    ax.plot(x, y, color=C_CD19, lw=1.2, zorder=3)
    ax.plot(x, y, ls="none", marker="o", ms=3.0, mfc="white", mec=C_CD19,
            mew=1.0, zorder=4, clip_on=False)

    ax.set_ylim(-0.9, 26)
    ax.set_yticks([0, 10, 20])
    ax.set_ylabel(f"CD19+ B cells\n({UNIT_CD19})")
    tidy(ax, last=True)

    ax.annotate("B-cell repopulation", (-128, 22), xytext=(-118, 24),
                fontsize=6.0, color=MUTED, ha="left", va="top")
    ax.annotate("sustained B-cell depletion", (40, 2.2), fontsize=6.2,
                color=C_CD19, ha="center", va="bottom", fontweight="bold")


# ---------------------------------------------------------------------------
# ASSEMBLY
# ---------------------------------------------------------------------------

def build_figure():
    set_style()
    fig = plt.figure(figsize=(7.09, 7.6))   # 180 mm wide, double-column
    gs = fig.add_gridspec(
        4, 1, height_ratios=[2.05, 1.05, 1.0, 0.82],
        hspace=0.18, left=0.105, right=0.775, top=0.955, bottom=0.062)

    ax_tl = fig.add_subplot(gs[0])
    ax_pred = fig.add_subplot(gs[1])
    ax_pr3 = fig.add_subplot(gs[2], sharex=ax_pred)
    ax_cd19 = fig.add_subplot(gs[3], sharex=ax_pred)

    draw_timeline(ax_tl)
    draw_prednisolone(ax_pred)
    draw_pr3(ax_pr3)
    draw_cd19(ax_cd19)

    for ax, letter in zip((ax_tl, ax_pred, ax_pr3, ax_cd19), "abcd"):
        panel_letter(ax, letter)

    # index-therapy callout above the timeline
    ax_tl.annotate(
        "first teclistamab dose", (4, -1.62), ha="left", va="bottom",
        fontsize=6.6, fontweight="bold", color=C_TEC, annotation_clip=False)

    return fig


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--xlsx", type=Path, default=None,
                   help="re-read the source workbook instead of embedded data")
    p.add_argument("--outdir", type=Path, default=Path(__file__).parent / "output")
    p.add_argument("--formats", nargs="+", default=["pdf", "png"],
                   choices=["pdf", "png", "tiff", "svg", "eps"])
    args = p.parse_args()

    if args.xlsx:
        global PREDNISOLONE, PR3, CD19
        PREDNISOLONE, PR3, CD19 = load_from_xlsx(args.xlsx)

    fig = build_figure()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for fmt in args.formats:
        out = args.outdir / f"figure1_therapy_timeline.{fmt}"
        kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if fmt == "tiff" else {}
        fig.savefig(out, **kw)
        print(f"wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
