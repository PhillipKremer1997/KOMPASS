# TEC_GPA (Figure 3) — Teclistamab in refractory pulmonary GPA

Day 0 = first teclistamab dose.

    figure/data/observations.csv   day, prednisolone, anti-PR3 IgG, CD19+
    figure/data/treatments.csv     agent, start/end day, dose label
    figure/data/timepoints.csv     imaging timepoints shown in panel A

`figure/data/` holds patient-level data and is **not** under version
control (see `.gitignore`), as are the rendered figures. Restore it from
the source spreadsheet (`Therapie.xlsx`) before running the script.

    python3 make_figure.py          # writes TEC_GPA.pdf / .png / .svg
    Rscript  make_figure.R           # writes TEC_GPA_R.pdf / .png
    Rscript  make_figure_ggplot.R    # writes TEC_GPA_ggplot.pdf / .png

Outputs
* `TEC_GPA.pdf` — vector, 180 mm wide, fonts embedded (Type 42); submission file
* `TEC_GPA.png` — 600 dpi raster, for drafts and slides
* `TEC_GPA.svg` — editable text, for import into BioRender / Illustrator

All assumptions are parameters at the top of the script (section 0) and are
listed in `TEC_GPA_legend.md`.

## The three versions

`make_figure.py` (matplotlib) and `make_figure.R` (base graphics) read the
same CSVs and draw the same figure on the same 7.09 x 7.35 in canvas; use
whichever is easier to maintain. Verified by rendering both at 600 dpi and
comparing pixel by pixel: every data mark, axis, rail, dose block and
shaded band lands within half a pixel. 1.4% of pixels differ by more than
32/255, which falls to 0.34% once a two-pixel alignment tolerance is
allowed - that residue is glyph rasterisation, as matplotlib's FreeType and
R's cairo set the same Liberation Sans with advance widths about 2% apart.

`make_figure_ggplot.R` (ggplot2 + patchwork) is a third route for whoever
prefers the grammar-of-graphics style. It draws the same figure but is
visually equivalent rather than pixel-identical, because ggplot's layout
engine places panels and margins its own way. Two consequences are visible:
what the other two draw outside a panel (the imaging row, the pulse arrows)
is inside the panel here, so those y ranges are a little wider; and
ggplot2 < 3.5 cannot draw minor tick marks, so the x axis has major ticks
only. It needs ggplot2 (>= 3.4) and patchwork; the other two scripts need
no packages beyond matplotlib and base R.

Panel D reaches 100 mg/day. An earlier 60 mg limit silently cut off the two
oral 100 mg doses given on rituximab days - worth remembering if the axis
is ever rescaled.

## Running it in RStudio

Open `TEC_GPA.Rproj`, then Source `make_figure.R` or
`make_figure_ggplot.R`. Each script locates its own folder (Rscript, the
Source button and `source()` from the console are all handled) and reads
`data/` from there, so the working directory does not matter; if the folder
is missing the script says so and names the two places it looked.

    install.packages(c("ggplot2", "patchwork"))   # only for the ggplot version

The project is set to UTF-8, which the dose labels need (they contain × and
→). R 4.2 or newer is recommended: earlier versions on Windows did not use
UTF-8 natively and mangle those two characters.

Fonts: the scripts ask for `sans`, which is Arial on Windows and Helvetica
on macOS - what journals want. Only on Linux with cairo is Liberation Sans
requested instead, because there `sans` means DejaVu Sans. Change `FAM` near
the top of either script to use a different face; it has to be a family the
output device knows, or a plain `pdf()` stops with "invalid font type".

Devices: PDF and PNG are opened by `open_pdf()` / `open_png()`, which try
quartz first on macOS, then cairo, then a plain `pdf()` / `png()`. That
order matters - a Mac without XQuartz has no cairo at all, and merely asking
`capabilities("cairo")` there loads the X11 module and warns. The last
resort cannot encode × and →, and says so.
