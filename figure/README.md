# Figure 3 — Teclistamab in refractory pulmonary GPA

Day 0 = first teclistamab dose.

    figure/data/observations.csv   day, prednisolone, anti-PR3 IgG, CD19+
    figure/data/treatments.csv     agent, start/end day, dose label
    figure/data/timepoints.csv     imaging timepoints shown in panel A

`figure/data/` holds patient-level data and is **not** under version
control (see `.gitignore`), as are the rendered figures. Restore it from
the source spreadsheet (`Therapie.xlsx`) before running the script.

    python3 make_figure.py     # writes Figure3.pdf / .png / .svg

Outputs
* `Figure3.pdf` — vector, 180 mm wide, fonts embedded (Type 42); submission file
* `Figure3.png` — 600 dpi raster, for drafts and slides
* `Figure3.svg` — editable text, for import into BioRender / Illustrator

All assumptions are parameters at the top of the script (section 0) and are
listed in `Figure3_legend.md`.
