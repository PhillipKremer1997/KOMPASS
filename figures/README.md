# Figure 1 — therapy timeline

`therapy_timeline.py` builds the multi-panel therapy/response figure from the
`Therapie.xlsx` record. Run it with no arguments to use the embedded data:

    python therapy_timeline.py
    python therapy_timeline.py --formats pdf png tiff --outdir output
    python therapy_timeline.py --xlsx /path/to/Therapie.xlsx   # re-read source

Output goes to `output/` as `figure1_therapy_timeline.{pdf,png,tiff}`. The PDF is
vector with TrueType-embedded fonts (`pdf.fonttype = 42`), which is what journals
ask for; the TIFF is LZW-compressed at 600 dpi.

## Layout

    (a) treatment timeline, one lane per therapy
    (b) oral prednisolone maintenance dose
    (c) anti-PR3 IgG, log scale
    (d) CD19+ B cells

All four panels share the x-axis (days from first teclistamab dose), so a
vertical read-off across panels is valid.

## Interpretation decisions worth reviewing

**Steroid pulses are separated from the maintenance dose.** Six of the recorded
prednisolone values fall on infusion days and are induction pulses or infusion
premedication rather than the oral maintenance dose:

| Day  | Dose   | Coincides with            |
|------|--------|---------------------------|
| −318 | 500 mg | induction                 |
| −312 | 100 mg | rituximab 1 g             |
| −296 | 100 mg | rituximab 1 g             |
| −291 | 500 mg | cyclophosphamide 1.0 g    |
| −97  | 500 mg | cyclophosphamide 1.4 g    |
| −93  | 250 mg | daratumumab               |

Plotting them on the same axis as the maintenance dose forces a 0–500 mg scale,
on which the taper from 60 mg to 2 mg/day is compressed into the bottom eighth of
the panel — this is what happens in the original figure. Moving them to panel (a)
lets panel (b) use a 0–60 mg scale and makes the taper the visible signal.

The classification is clinical, not mechanical. It lives in the `PULSE_DAYS` set
at the top of the script; removing a day from that set returns it to panel (b).
Two entries are judgement calls worth a second opinion:

- **Day −30 (30 mg)** is kept as *maintenance* even though 1.4 g
  cyclophosphamide was given the same day, because 30 mg fits the surrounding
  taper (40 → 30 → 20 mg) rather than a pulse.
- **Days −82 (60 mg) and −70 (40 mg)** are kept as *maintenance*: they fall
  between daratumumab doses and read as a genuine re-escalation for the relapse,
  not as premedication.

**The last recorded dose is carried to the right edge** of panel (b) as a step,
since an oral dose persists between clinic visits. Measured points are marked
with open circles, so carried-forward segments are distinguishable from measured
ones.

**Zero anti-PR3 values** are plotted at `PR3_FLOOR` (0.6) with open symbols
because a log axis cannot show zero. They are values reported as 0, i.e. below
the assay cut-off, not true zeros.

## Editable parameters

Everything a co-author is likely to change sits in one block at the top of the
script: units, assay cut-off, glucocorticoid target, scan timepoints, the
teclistamab window and dose days, and the x-limits. Colours are defined just
below as named constants (`C_RTX`, `C_CYC`, …), chosen to stay distinguishable in
greyscale and for the common forms of colour-vision deficiency.

## Data

The day offsets embedded in the script are relative to the first teclistamab
dose; no calendar dates or identifiers are stored here. The source workbook is
deliberately not committed.
