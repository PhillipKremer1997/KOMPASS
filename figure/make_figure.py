#!/usr/bin/env python3
"""
Figure 3 - BCMA-directed T-cell engager therapy in refractory pulmonary GPA.

Merges the manuscript's Figure 3A (treatment regimen, prednisolone) and
Figure 3B (anti-PR3 IgG, CD19+ B cells) into one panel stack on a shared,
to-scale time axis.

Panel A  treatment timeline (swimmer lanes)
Panel B  anti-PR3 IgG  (broken y-axis)
Panel C  CD19+ B cells
Panel D  glucocorticoid exposure

Style reference: NEJM / Nature Medicine / Ann Rheum Dis reports on CAR-T and
T-cell-engager therapy in autoimmune disease.

Output: Figure3.pdf (vector, for submission), .png (600 dpi), .svg (BioRender)
"""
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

# --------------------------------------------------------------------------
# 0.  PARAMETERS TO BE CONFIRMED BY THE CLINICAL TEAM        <<< CHECK ME
# --------------------------------------------------------------------------
PR3_ULN      = 2.0     # ULN of the anti-PR3 immunosorbent assay (manuscript)
PR3_LLOQ     = 0.4     # limit of quantification; reported "0" is drawn here
PULSE_CUTOFF = 250     # documented doses >= this value = i.v. GC pulse
# chest CT timepoints, Figure 1A-C: 03/2025, 12/2025, 07/2026.  Only the month
# is documented; A is anchored to the baseline date.       <<< exact dates?
SCANS        = {"A": -318, "B": -30, "C": 180}
# teclistamab s.c. step-up per manuscript: d1 0.06, d3 0.3, d5/d12/d19/d26
# 1.5 mg/kg  ->  days 0, 2, 4, 11, 18, 25 relative to the first dose.
# NB the source spreadsheet notes "6x until 19-Feb-2026" (= day +36). <<< check
TEC  = [0, 2, 4, 11, 18, 25]
TEC_WINDOW   = (0, 25)

OUT = Path(__file__).resolve().parent

# --------------------------------------------------------------------------
# 1.  DATA  (read from ./data, which is not under version control)
# --------------------------------------------------------------------------
DATA = OUT / "data"
if not DATA.is_dir():
    raise SystemExit(
        f"missing data directory: {DATA}\n"
        "Patient-level data are kept out of the repository. Restore "
        "observations.csv, treatments.csv and timepoints.csv from the "
        "source spreadsheet before running this script.")


def read_csv(name):
    with open(DATA / name, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def num(v):
    return np.nan if v in (None, "") else float(v)


obs  = sorted(read_csv("observations.csv"), key=lambda r: float(r["day"]))
day  = np.array([float(r["day"]) for r in obs])
pred = np.array([num(r["prednisolone_mg"]) for r in obs])
pr3  = np.array([num(r["anti_pr3_u_ml"]) for r in obs])
cd19 = np.array([num(r["cd19_cells_ul"]) for r in obs])

TX = read_csv("treatments.csv")
lanes = list(dict.fromkeys(r["agent"] for r in TX))          # order of appearance
TIMEPOINTS = read_csv("timepoints.csv")

tec_days = [float(r["start_day"]) for r in TX if r["agent"] == lanes[-1]]
TEC_WINDOW = (min(tec_days), max(tec_days))

XLIM, XTICKS = (-338, 208), np.arange(-300, 201, 60)

# --------------------------------------------------------------------------
# 2.  STYLE
# --------------------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.linewidth": 0.7,
    "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 2.8, "ytick.major.size": 2.8,
    "xtick.direction": "out", "ytick.direction": "out",
    "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})
C = dict(plex="#9E6BA8", cyc="#E1A33C", rtx="#4C9F70", ava="#8C8C8C",
         dara="#3E8FC4", tec="#C0392B", pr3="#B2182B", cd19="#0E7C7B",
         pred="#2B5FA8", grey="#595959", rail="#EDEDED")

# --------------------------------------------------------------------------
# 3.  HELPERS
# --------------------------------------------------------------------------
def block(ax, x, lane, color, w=5.5, h=0.34):
    ax.add_patch(Rectangle((x - w / 2, lane - h / 2), w, h,
                           lw=0, facecolor=color, zorder=3))

def bar(ax, x0, x1, lane, color, h=0.34, alpha=1.0):
    ax.add_patch(Rectangle((x0, lane - h / 2), x1 - x0, h,
                           lw=0, facecolor=color, alpha=alpha, zorder=3))

def lane_label(ax, x, lane, text, ha="left", color=None, dx=6):
    ax.text(x + (dx if ha == "left" else -dx), lane, text, ha=ha, va="center",
            fontsize=5.9, color=color or C["grey"], zorder=4)

def step_series(x, y):
    """LOCF step curve for the oral dose; an i.v. pulse always starts a new
    oral course, so the next documented dose is carried backwards to it."""
    m = ~np.isnan(y)
    xs, ys = x[m], y[m]
    ox, oy = xs[ys < PULSE_CUTOFF], ys[ys < PULSE_CUTOFF]
    px, py = xs[ys >= PULSE_CUTOFF], ys[ys >= PULSE_CUTOFF]
    for p in px:
        nxt = ox[ox > p]
        if len(nxt):
            ox = np.append(ox, p); oy = np.append(oy, oy[ox[:-1] == nxt[0]][0])
    o = np.argsort(ox)
    return ox[o], oy[o], px, py

def letter(ax, s, dy=1.0):
    ax.text(-0.105, dy, s, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="left")

def timemarks(ax):
    ax.axvspan(*TEC_WINDOW, color=C["tec"], alpha=0.055, lw=0, zorder=0)
    ax.axvline(0, color=C["tec"], lw=0.8, ls=(0, (3.5, 2)), zorder=1)

# --------------------------------------------------------------------------
# 4.  FIGURE
# --------------------------------------------------------------------------
fig = plt.figure(figsize=(7.09, 7.0))          # 180 mm journal width
gs = fig.add_gridspec(5, 1, height_ratios=[2.6, 0.45, 1.5, 0.95, 1.15],
                      hspace=0.30, left=0.165, right=0.975,
                      top=0.925, bottom=0.072)
axT = fig.add_subplot(gs[0])
axU = fig.add_subplot(gs[1], sharex=axT)
axL = fig.add_subplot(gs[2], sharex=axT)
axB = fig.add_subplot(gs[3], sharex=axT)
axP = fig.add_subplot(gs[4], sharex=axT)

# ---------- A  treatment timeline -----------------------------------------
y = {n: len(lanes) - 1 - i for i, n in enumerate(lanes)}
LANE_COLOUR = dict(zip(lanes, ["plex", "cyc", "rtx", "ava", "dara", "tec"]))

for n in lanes:
    axT.plot(XLIM, [y[n]] * 2, color=C["rail"], lw=2.2, zorder=0,
             solid_capstyle="butt")

# a dose block is never wider than the gap to the next dose of the same agent
widths = {}
for n in lanes:
    d = sorted(float(r["start_day"]) for r in TX
               if r["agent"] == n and not r["end_day"])
    gaps = [b - a for a, b in zip(d, d[1:])] or [99]
    widths[n] = min(5.5, max(2.6, min(gaps) * 1.2))

for r in TX:
    n, col = r["agent"], C[LANE_COLOUR[r["agent"]]]
    x0 = float(r["start_day"])
    if r["end_day"]:
        x1 = float(r["end_day"])
        pad = 2.5 if x1 - x0 < 30 else 0     # pad short clusters, not long courses
        bar(axT, x0 - pad, x1 + pad, y[n], col, alpha=0.9 if pad == 0 else 1.0)
    else:
        block(axT, x0, y[n], col, w=widths[n])

    if not r["label"]:
        continue
    if r["label_align"] == "inside":
        axT.text((x0 + float(r["end_day"])) / 2, y[n], r["label"], ha="center",
                 va="center", fontsize=5.9, color="white", zorder=4)
    else:
        lane_label(axT, float(r["label_day"]), y[n], r["label"],
                   ha=r["label_align"],
                   color=C["tec"] if n == lanes[-1] else None)

axT.set_yticks([y[n] for n in lanes]); axT.set_yticklabels(lanes)
axT.set_ylim(-0.55, len(lanes) - 0.35)
for s in axT.spines.values():
    s.set_visible(False)
axT.tick_params(axis="y", length=0, pad=3)
axT.tick_params(axis="x", labelbottom=False, length=0)
timemarks(axT)
letter(axT, "A", dy=1.16)

top = len(lanes) - 0.42
for r in TIMEPOINTS:
    x = float(r["day"])
    axT.plot(x, top, marker="v", ms=4.2, color=C["grey"], clip_on=False, zorder=5)
    axT.text(x, top + 0.13, r["label"], ha="center", va="bottom", fontsize=6.5,
             fontweight="bold", color=C["grey"])
axT.text(XLIM[0] - 6, top, TIMEPOINTS[0]["kind"], ha="right", va="center",
         fontsize=6.5, color=C["grey"])
axT.text(4, top + 0.13, "Teclistamab, day 0", ha="left", va="bottom",
         fontsize=6.5, color=C["tec"], fontweight="bold")

# ---------- B  anti-PR3 IgG (broken axis) ---------------------------------
m = ~np.isnan(pr3)
xp, yp = day[m], pr3[m]
zero = yp == 0
yplot = np.where(zero, PR3_LLOQ, yp)
for ax in (axU, axL):
    ax.plot(xp, yplot, color=C["pr3"], lw=1.2, zorder=3)
    ax.plot(xp[~zero], yplot[~zero], "o", ms=3.2, mfc=C["pr3"], mec="white",
            mew=0.5, ls="none", zorder=4)
    ax.plot(xp[zero], yplot[zero], "o", ms=3.2, mfc="white", mec=C["pr3"],
            mew=0.9, ls="none", zorder=4)
    timemarks(ax)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
axL.axhspan(-1.2, PR3_ULN, color=C["pr3"], alpha=0.06, lw=0, zorder=0)
axL.axhline(PR3_ULN, color=C["pr3"], lw=0.6, ls=(0, (2, 2)), alpha=.85, zorder=2)
axL.text(XLIM[0] + 6, PR3_ULN + 0.9, "normal < 2 U/mL", ha="left",
         va="bottom", fontsize=6, color=C["pr3"])

axU.set_ylim(142, 162); axU.set_yticks([150])
axL.set_ylim(-1.2, 31); axL.set_yticks([0, 10, 20, 30])
axU.tick_params(axis="x", labelbottom=False, length=0)
axL.tick_params(axis="x", labelbottom=False)
axU.spines["bottom"].set_visible(False)
for ax, ys in ((axU, (-0.16, 0.16)), (axL, (0.94, 1.14))):
    ax.plot((-0.007, 0.007), ys, transform=ax.transAxes, color="k",
            clip_on=False, lw=0.7)
axL.set_ylabel("Anti-PR3 IgG (U/mL)")
axL.yaxis.set_label_coords(-0.115, 0.62)
letter(axU, "B", dy=1.5)

# ---------- C  CD19+ B cells ----------------------------------------------
mb = ~np.isnan(cd19)
axB.plot(day[mb], cd19[mb], color=C["cd19"], lw=1.2, zorder=3)
axB.plot(day[mb], cd19[mb], "s", ms=3.0, mfc=C["cd19"], mec="white", mew=0.5,
         ls="none", zorder=4)
axB.set_ylim(-1.8, 28); axB.set_yticks([0, 10, 20])
axB.set_ylabel("CD19+ B cells\n(cells/µL)")
axB.yaxis.set_label_coords(-0.115, 0.5)
axB.spines["top"].set_visible(False); axB.spines["right"].set_visible(False)
axB.tick_params(axis="x", labelbottom=False)
timemarks(axB)
axB.text(-140, 22, "B-cell repopulation\nbefore relapse", ha="right",
         va="center", fontsize=5.9, color=C["cd19"])
letter(axB, "C", dy=1.06)

# ---------- D  glucocorticoids --------------------------------------------
ox, oy, px, py = step_series(day, pred)
xs = np.append(ox, day[-1]); ys = np.append(oy, oy[-1])
axP.fill_between(xs, 0, ys, step="post", color=C["pred"], alpha=0.18, lw=0, zorder=2)
axP.step(xs, ys, where="post", color=C["pred"], lw=1.2, zorder=3)
axP.plot(ox, oy, "o", ms=2.6, mfc=C["pred"], mec="white", mew=0.4, ls="none",
         zorder=4)

# i.v. pulses are shown as arrows in a strip above the daily-dose axis;
# pulses given within 10 days of each other are pooled into one arrow
groups, cur = [], [(px[0], py[0])] if len(px) else []
for x, v in zip(px[1:], py[1:]):
    if x - cur[-1][0] <= 10:
        cur.append((x, v))
    else:
        groups.append(cur); cur = [(x, v)]
if cur:
    groups.append(cur)

YMAX, YARR = 60, 66
for g in groups:
    gx = float(np.mean([x for x, _ in g]))
    axP.plot([gx, gx], [YARR - 4.5, YARR], color=C["pred"], lw=0.7,
             clip_on=False, zorder=4)
    axP.plot(gx, YARR - 5.5, marker="v", ms=4.2, color=C["pred"],
             clip_on=False, zorder=4)
    axP.text(gx, YARR + 0.8, ", ".join(str(int(v)) for _, v in g), ha="center",
             va="bottom", fontsize=5.8, color=C["pred"], clip_on=False)
axP.text(XLIM[1] - 4, YARR + 0.8, "i.v. methylprednisolone pulses (mg)",
         ha="right", va="bottom", fontsize=6, color=C["pred"], clip_on=False)

axP.set_ylim(0, YMAX); axP.set_yticks([0, 20, 40, 60])
axP.set_ylabel("Prednisolone\n(mg/day)")
axP.yaxis.set_label_coords(-0.115, 0.40)
axP.spines["top"].set_visible(False); axP.spines["right"].set_visible(False)
timemarks(axP)
letter(axP, "D", dy=1.06)

axP.set_xlim(*XLIM); axP.set_xticks(XTICKS)
axP.set_xticks(np.arange(-330, 201, 30), minor=True)
axP.tick_params(axis="x", which="minor", length=1.6, width=0.6)
for ax in (axT, axU, axL, axB):
    ax.tick_params(axis="x", which="both", length=0, labelbottom=False)
axP.set_xlabel("Days relative to first teclistamab dose")

for ext in ("pdf", "png", "svg"):
    fig.savefig(OUT / f"Figure3.{ext}", dpi=600 if ext == "png" else None,
                bbox_inches="tight", facecolor="white")
print("ok")
