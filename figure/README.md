# Figure 3 — Teclistamab in refractory pulmonary GPA

Day 0 = first teclistamab dose.

    figure/data/observations.csv   day, prednisolone, anti-PR3 IgG, CD19+
    figure/data/treatments.csv     agent, start/end day, dose label
    figure/data/timepoints.csv     imaging timepoints shown in panel A

`figure/data/` holds patient-level data and is **not** under version
control (see `.gitignore`), as are the rendered figures. Restore it from
the source spreadsheet (`Therapie.xlsx`) before running the script.

    python3 make_figure.py     # writes Figure3.pdf / .png / .svg
    Rscript  make_figure.R      # writes Figure3_R.pdf / .png

Outputs
* `Figure3.pdf` — vector, 180 mm wide, fonts embedded (Type 42); submission file
* `Figure3.png` — 600 dpi raster, for drafts and slides
* `Figure3.svg` — editable text, for import into BioRender / Illustrator

All assumptions are parameters at the top of the script (section 0) and are
listed in `Figure3_legend.md`.

## Python and R versions

`make_figure.py` (matplotlib) and `make_figure.R` (base graphics) read the
same CSVs and draw the same figure on the same 7.09 x 7.0 in canvas; use
whichever is easier to maintain. Verified by rendering both at 600 dpi and
comparing pixel by pixel: every data mark, axis, rail, dose block and
shaded band lands within half a pixel. 1.4% of pixels differ by more than
32/255, which falls to 0.34% once a two-pixel alignment tolerance is
allowed - that residue is glyph rasterisation, as matplotlib's FreeType and
R's cairo set the same Liberation Sans with advance widths about 2% apart.

The R script uses base graphics rather than ggplot2 on purpose: the axis
break in panel B and the pulse arrows drawn outside panel D need
device-level coordinate control that a grid-based layout makes awkward.
