#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Figure 3 - BCMA-directed T-cell engager therapy in refractory pulmonary GPA.
#
# ggplot2 + patchwork version.  Same data, same design and same canvas as
# make_figure.R (base graphics) and make_figure.py (matplotlib); the layout
# is produced by ggplot's own engine, so it is visually equivalent rather
# than pixel-identical to those two.
#
# Packages: ggplot2 (>= 3.4), patchwork.  Nothing else.
#   install.packages(c("ggplot2", "patchwork"))
#
# Differences from the base-graphics version, all forced by ggplot's layout
# model: annotations that sit outside a panel there (the imaging row, the
# pulse arrows) are inside the panel here, so the y ranges are a little
# wider; and ggplot2 < 3.5 cannot draw minor tick marks, so the x axis has
# major ticks only.
#
# Output: Figure3_ggplot.pdf (vector, for submission), .png (600 dpi)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ---- 0. parameters confirmed by the clinical team --------------------------
PR3_ULN      <- 2.0   # ULN of the anti-PR3 immunosorbent assay (manuscript)
PR3_LLOQ     <- 0.4   # limit of quantification; reported "0" is drawn here
# Confirmed by the authors: doses >= 250 mg are i.v. pulses (drawn as arrows),
# the 100 mg on rituximab days are oral, the daily dose is carried forward
# between documented changes, and after a pulse the curve resumes at the next
# documented oral dose.
PULSE_CUTOFF <- 250

# ---- where is this script, and where is its data? --------------------------
# Resolves for Rscript, for RStudio's Source button, and for source() from the
# console, then falls back to the working directory.
script_dir <- function() {
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  for (i in rev(seq_len(sys.nframe()))) {           # source() keeps it here
    of <- get0("ofile", envir = sys.frame(i), ifnotfound = NULL)
    if (is.character(of) && length(of) == 1L && nzchar(of))
      return(dirname(normalizePath(of)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(try(rstudioapi::isAvailable(), silent = TRUE))) {
    p <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)
    if (!inherits(p, "try-error") && is.character(p) && nzchar(p))
      return(dirname(normalizePath(p)))
  }
  getwd()
}

OUT  <- script_dir()
DATA <- file.path(OUT, "data")
if (!dir.exists(DATA)) DATA <- file.path(getwd(), "data")
if (!dir.exists(DATA))
  stop("cannot find the data folder.\n",
       "  looked next to the script: ", file.path(OUT, "data"), "\n",
       "  and in the working dir:    ", file.path(getwd(), "data"), "\n\n",
       "Put observations.csv, treatments.csv and timepoints.csv into a folder\n",
       "named 'data' next to this script and run it again. In RStudio, open\n",
       "KOMPASS-Figure3.Rproj first, then press Source.\n",
       "Patient-level data are deliberately kept out of the repository.",
       call. = FALSE)

# ---- 1. data ---------------------------------------------------------------
rd <- function(f) read.csv(file.path(DATA, f), colClasses = "character",
                           encoding = "UTF-8", check.names = FALSE)
n  <- function(v) suppressWarnings(as.numeric(ifelse(v == "", NA, v)))

obs <- rd("observations.csv")
obs <- data.frame(day  = n(obs$day),
                  pred = n(obs$prednisolone_mg),
                  pr3  = n(obs$anti_pr3_u_ml),
                  cd19 = n(obs$cd19_cells_ul))
obs <- obs[order(obs$day), ]

TX    <- rd("treatments.csv")
lanes <- unique(TX$agent)
TP    <- rd("timepoints.csv")

TX$start_day <- n(TX$start_day)
TX$end_day   <- n(TX$end_day)
TX$label_day <- n(TX$label_day)
TX$ypos      <- length(lanes) - match(TX$agent, lanes)
TP$day       <- n(TP$day)

TEC_WINDOW <- range(TX$start_day[TX$agent == tail(lanes, 1)])
XLIM   <- c(-338, 208)
XTICKS <- seq(-300, 200, 60)

# ---- 2. style --------------------------------------------------------------
CO <- c(plex = "#9E6BA8", cyc = "#E1A33C", rtx = "#4C9F70", ava = "#8C8C8C",
        dara = "#3E8FC4", tec = "#C0392B", pr3 = "#B2182B", cd19 = "#0E7C7B",
        pred = "#2B5FA8", grey = "#595959", rail = "#EDEDED")
TX$colour <- CO[setNames(c("plex", "cyc", "rtx", "ava", "dara", "tec"),
                         lanes)[TX$agent]]
# Journals ask for Arial or Helvetica. "sans" is Arial on Windows, Helvetica on
# macOS quartz, and Helvetica in a plain pdf() - all fine. Only cairo on Linux
# maps "sans" to DejaVu Sans, so there Liberation Sans (metric-compatible with
# Arial) is asked for instead. The family has to be one the chosen device
# knows: a plain pdf() errors with "invalid font type" on anything else.
FAM <- if (Sys.info()[["sysname"]] == "Linux" &&
           isTRUE(suppressWarnings(capabilities("cairo")))) "Liberation Sans" else "sans"
MU  <- intToUtf8(0xB5)        # micro sign, independent of the source encoding

lw <- function(pt) pt / .pt   # line width in pt -> ggplot linewidth
gs <- function(pt) pt / .pt   # text size in pt   -> geom_text size

theme_fig <- function(...) {
  base <- theme_void() +
    theme(text        = element_text(family = FAM, size = 7, colour = "black"),
          axis.text.y = element_text(size = 7, hjust = 1,
                                     margin = margin(r = 2.2)),
          axis.text.x = element_text(size = 7, margin = margin(t = 1.4)),
          axis.title.y = element_text(size = 7.5, angle = 90,
                                      margin = margin(r = 3.2)),
          axis.title.x = element_text(size = 7.5, margin = margin(t = 3.5)),
          axis.line    = element_line(linewidth = lw(0.7), colour = "black",
                                      lineend = "square"),
          axis.ticks   = element_line(linewidth = lw(0.7), colour = "black"),
          axis.ticks.length = unit(2.8, "pt"),
          # the top and bottom margins double as the gap between panels and
          # as the room patchwork needs for the panel tag
          plot.margin  = margin(11, 4, 11, 2),
          plot.tag     = element_text(size = 9, face = "bold", family = FAM,
                                      hjust = 0, vjust = 1),
          legend.position = "none")
  # overrides are added rather than passed through, so a caller may replace an
  # element the base theme already sets
  if (...length()) base + theme(...) else base
}

# name = NULL everywhere, so a panel never inherits an axis title from the
# name of the variable it happens to be mapped to
scale_x <- function(labels = FALSE, name = NULL)
  scale_x_continuous(name = name, breaks = XTICKS, minor_breaks = NULL,
                     labels = if (isFALSE(labels)) NULL else waiver())

scale_y <- function(name = NULL, ...)
  scale_y_continuous(name = name, expand = c(0, 0), ...)

# limits belong on the coordinate system, not the scale: a scale limit deletes
# the rows outside it (and warns), while a coord limit clips the drawing, so a
# line leaving the panel is still drawn up to the edge
zoom <- function(ylim) coord_cartesian(xlim = XLIM, ylim = ylim, expand = FALSE)

# day 0 and the teclistamab window, drawn under everything else
timemarks <- function()
  list(annotate("rect", xmin = TEC_WINDOW[1], xmax = TEC_WINDOW[2],
                ymin = -Inf, ymax = Inf, fill = CO[["tec"]], alpha = 0.055),
       annotate("segment", x = 0, xend = 0, y = -Inf, yend = Inf,
                colour = CO[["tec"]], linewidth = lw(0.8),
                linetype = "42"))

# ---- 3. panel A : treatment timeline ---------------------------------------
rails  <- data.frame(y = length(lanes) - seq_along(lanes))
points_tx <- subset(TX, is.na(end_day))
spans_tx  <- subset(TX, !is.na(end_day))
# a dose block is never wider than the gap to the next dose of the same agent
w <- sapply(lanes, function(l) {
  d <- sort(points_tx$start_day[points_tx$agent == l])
  min(5.5, max(2.6, min(if (length(d) > 1) diff(d) else 99) * 1.2))
})
points_tx$w <- w[points_tx$agent]
# short clusters get padded, long courses do not
spans_tx$pad <- ifelse(spans_tx$end_day - spans_tx$start_day < 30, 2.5, 0)
lab_tx  <- subset(TX, label != "" & label_align != "inside")
lab_tx$dx <- ifelse(lab_tx$label_align == "left", 6, -6)
lab_tx$h  <- ifelse(lab_tx$label_align == "left", 0, 1)
lab_in  <- subset(TX, label_align == "inside")
TOP <- length(lanes) - 0.42

pA <- ggplot() +
  geom_segment(data = rails, aes(x = XLIM[1], xend = XLIM[2], y = y, yend = y),
               colour = CO[["rail"]], linewidth = lw(2.2), lineend = "butt") +
  timemarks() +
  geom_rect(data = spans_tx,
            aes(xmin = start_day - pad, xmax = end_day + pad,
                ymin = ypos - 0.17, ymax = ypos + 0.17),
            fill = spans_tx$colour, alpha = ifelse(spans_tx$pad == 0, 0.9, 1)) +
  geom_rect(data = points_tx,
            aes(xmin = start_day - w / 2, xmax = start_day + w / 2,
                ymin = ypos - 0.17, ymax = ypos + 0.17),
            fill = points_tx$colour) +
  geom_text(data = lab_in, aes(x = (start_day + end_day) / 2, y = ypos,
                               label = label),
            colour = "white", size = gs(5.9), family = FAM) +
  geom_text(data = lab_tx, aes(x = label_day + dx, y = ypos, label = label,
                               hjust = h),
            colour = ifelse(lab_tx$agent == tail(lanes, 1), CO[["tec"]],
                            CO[["grey"]]),
            size = gs(5.9), family = FAM) +
  geom_point(data = TP, aes(x = day, y = TOP), shape = 25, size = gs(4.2),
             colour = CO[["grey"]], fill = CO[["grey"]]) +
  geom_text(data = TP, aes(x = day, y = TOP + 0.13, label = label),
            vjust = 0, fontface = "bold", size = gs(6.5),
            colour = CO[["grey"]], family = FAM) +
  annotate("text", x = 4, y = TOP + 0.13, label = "Teclistamab, day 0",
           hjust = 0, vjust = 0, fontface = "bold", size = gs(6.5),
           colour = CO[["tec"]], family = FAM) +
  scale_x() +
  scale_y(breaks = c(rev(rails$y), TOP), labels = c(rev(lanes), TP$kind[1])) +
  zoom(c(-0.55, length(lanes) + 0.15)) +
  theme_fig(axis.line = element_blank(), axis.ticks = element_blank()) +
  # one colour per break, so the imaging row reads as grey. ggplot warns that
  # vectorised element_text input is unsupported; it works, and the warning is
  # silenced here rather than printed on every run
  theme(axis.text.y = suppressWarnings(
    element_text(size = 7, hjust = 1, margin = margin(r = 2.2),
                 colour = c(rep("black", length(lanes)), CO[["grey"]]))))

# ---- 4. panel B : anti-PR3 IgG (broken axis) -------------------------------
pr3 <- subset(obs, !is.na(pr3))
pr3$zero <- pr3$pr3 == 0
pr3$y    <- ifelse(pr3$zero, PR3_LLOQ, pr3$pr3)

pr3_layers <- function() list(
  timemarks(),
  geom_line(data = pr3, aes(day, y), colour = CO[["pr3"]], linewidth = lw(1.2)),
  geom_point(data = subset(pr3, !zero), aes(day, y), shape = 21,
             fill = CO[["pr3"]], colour = "white", stroke = lw(0.5),
             size = gs(3.2)),
  geom_point(data = subset(pr3, zero), aes(day, y), shape = 21,
             fill = "white", colour = CO[["pr3"]], stroke = lw(0.9),
             size = gs(3.2)))

pB_up <- ggplot() + pr3_layers() +
  scale_x() +
  scale_y(breaks = 150) +
  zoom(c(142, 162)) +
  theme_fig(axis.line.x = element_blank(), axis.ticks.x = element_blank()) +
  annotate("segment", x = XLIM[1], xend = XLIM[1] + 7, y = 142.8,
           yend = 145.6, linewidth = lw(0.7))

pB_lo <- ggplot() +
  annotate("rect", xmin = XLIM[1], xmax = XLIM[2], ymin = -1.2,
           ymax = PR3_ULN, fill = CO[["pr3"]], alpha = 0.06) +
  pr3_layers() +
  annotate("segment", x = XLIM[1], xend = XLIM[2], y = PR3_ULN, yend = PR3_ULN,
           colour = CO[["pr3"]], linewidth = lw(0.6), linetype = "22") +
  annotate("text", x = XLIM[1] + 6, y = PR3_ULN + 0.9, hjust = 0, vjust = 0,
           label = "normal < 2 U/mL", size = gs(6), colour = CO[["pr3"]],
           family = FAM) +
  annotate("segment", x = XLIM[1], xend = XLIM[1] + 7, y = 29.4, yend = 31.4,
           linewidth = lw(0.7)) +
  scale_x() +
  scale_y(name = "Anti-PR3 IgG (U/mL)", breaks = c(0, 10, 20, 30)) +
  zoom(c(-1.2, 31)) +
  theme_fig(axis.ticks.x = element_blank())

# ---- 5. panel C : CD19+ B cells --------------------------------------------
cd19 <- subset(obs, !is.na(cd19))

pC <- ggplot() + timemarks() +
  geom_line(data = cd19, aes(day, cd19), colour = CO[["cd19"]],
            linewidth = lw(1.2)) +
  geom_point(data = cd19, aes(day, cd19), shape = 22, fill = CO[["cd19"]],
             colour = "white", stroke = lw(0.5), size = gs(3.0)) +
  annotate("text", x = -140, y = 22, hjust = 1, vjust = 0.5, lineheight = 1.2,
           label = "B-cell repopulation\nbefore relapse", size = gs(5.9),
           colour = CO[["cd19"]], family = FAM) +
  scale_x() +
  scale_y(name = paste0("CD19+ B cells\n(cells/", MU, "L)"), breaks = c(0, 10, 20)) +
  zoom(c(-1.8, 28)) +
  theme_fig(axis.ticks.x = element_blank())

# ---- 6. panel D : glucocorticoids ------------------------------------------
# LOCF step curve for the oral dose; an i.v. pulse always starts a new oral
# course, so the next documented dose is carried backwards to it
p  <- subset(obs, !is.na(pred))
ox <- p$day[p$pred <  PULSE_CUTOFF]; oy <- p$pred[p$pred <  PULSE_CUTOFF]
px <- p$day[p$pred >= PULSE_CUTOFF]; py <- p$pred[p$pred >= PULSE_CUTOFF]
for (q in px) {
  nxt <- ox[ox > q]
  if (length(nxt)) { v <- oy[ox == nxt[1]][1]; ox <- c(ox, q); oy <- c(oy, v) }
}
o  <- order(ox); ox <- ox[o]; oy <- oy[o]
xs <- c(ox, max(obs$day)); ys <- c(oy, tail(oy, 1))
k  <- length(xs)
stp <- data.frame(x = rep(xs, each = 2)[-1], y = rep(ys, each = 2)[-(2 * k)])

# pulses given within 10 days of each other are pooled into one arrow
grp   <- cumsum(c(TRUE, diff(px) > 10))
pulse <- do.call(rbind, lapply(split(seq_along(px), grp), function(g)
  data.frame(x = mean(px[g]),
             label = paste(as.integer(py[g]), collapse = ", "))))
YARR <- 115.5

pD <- ggplot() + timemarks() +
  geom_polygon(data = rbind(data.frame(x = stp$x[1], y = 0), stp,
                            data.frame(x = tail(stp$x, 1), y = 0)),
               aes(x, y), fill = CO[["pred"]], alpha = 0.18) +
  geom_line(data = stp, aes(x, y), colour = CO[["pred"]],
            linewidth = lw(1.2)) +
  geom_point(data = data.frame(x = ox, y = oy), aes(x, y), shape = 21,
             fill = CO[["pred"]], colour = "white", stroke = lw(0.4),
             size = gs(2.6)) +
  geom_segment(data = pulse, aes(x = x, xend = x, y = YARR - 7.9, yend = YARR),
               colour = CO[["pred"]], linewidth = lw(0.7)) +
  geom_point(data = pulse, aes(x = x, y = YARR - 9.6), shape = 25,
             colour = CO[["pred"]], fill = CO[["pred"]], size = gs(4.2)) +
  geom_text(data = pulse, aes(x = x, y = YARR + 1.4, label = label), vjust = 0,
            size = gs(5.8), colour = CO[["pred"]], family = FAM) +
  annotate("text", x = XLIM[2] - 4, y = YARR + 1.4, hjust = 1, vjust = 0,
           label = "i.v. methylprednisolone pulses (mg)", size = gs(6),
           colour = CO[["pred"]], family = FAM) +
  scale_x(labels = TRUE, name = "Days relative to first teclistamab dose") +
  scale_y(name = "Prednisolone\n(mg/day)", breaks = seq(0, 100, 25)) +
  zoom(c(0, 122)) +
  theme_fig()

# ---- 7. assemble -----------------------------------------------------------
fig <- pA / pB_up / pB_lo / pC / pD +
  plot_layout(heights = c(2.81, 0.45, 1.5, 0.95, 1.74)) +
  plot_annotation(tag_levels = list(c("A", "B", "", "C", "D")))

# ---- output devices --------------------------------------------------------
# macOS first, on purpose: quartz handles UTF-8 (the labels contain × and →)
# and needs no XQuartz, while asking for cairo on a Mac without XQuartz loads
# the X11 module and warns. Linux and Windows then use cairo, and a plain
# pdf() is the last resort - it cannot encode the arrow glyph.
open_pdf <- function(file, w, h) {
  if (capabilities("aqua")) {
    ok <- tryCatch({ grDevices::quartz(file = file, type = "pdf",
                                       width = w, height = h); TRUE },
                   error = function(e) FALSE)
    if (ok) return(invisible("quartz"))
  }
  if (isTRUE(suppressWarnings(capabilities("cairo")))) {
    grDevices::cairo_pdf(file, width = w, height = h)
    return(invisible("cairo"))
  }
  grDevices::pdf(file, width = w, height = h, useDingbats = FALSE,
                 encoding = "ISOLatin1")
  warning("no quartz or cairo PDF device: × and → will not render in the PDF",
          call. = FALSE)
  invisible("pdf")
}

open_png <- function(file, w, h, res = 600) {
  a <- list(filename = file, width = w, height = h, units = "in",
            res = res, bg = "white")
  if (capabilities("aqua")) {
    ok <- tryCatch({ do.call(grDevices::png, c(a, type = "quartz")); TRUE },
                   error = function(e) FALSE)
    if (ok) return(invisible("quartz"))
  }
  if (isTRUE(suppressWarnings(capabilities("cairo")))) {
    do.call(grDevices::png, c(a, type = "cairo"))
    return(invisible("cairo"))
  }
  do.call(grDevices::png, a)
  invisible("default")
}

open_pdf(file.path(OUT, "Figure3_ggplot.pdf"), 7.09, 7.35)
print(fig); invisible(dev.off())
open_png(file.path(OUT, "Figure3_ggplot.png"), 7.09, 7.35)
print(fig); invisible(dev.off())
cat("ok\n")
