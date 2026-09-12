#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Figure 3 - BCMA-directed T-cell engager therapy in refractory pulmonary GPA.
#
# Merges the manuscript's Figure 3A (treatment regimen, prednisolone) and
# Figure 3B (anti-PR3 IgG, CD19+ B cells) into one panel stack on a shared,
# to-scale time axis.
#
# Base-R port of make_figure.py; both scripts read the same CSVs in ./data
# and produce the same figure.  Base graphics rather than ggplot2, because
# the axis break in panel B and the pulse arrows drawn outside panel D need
# device-level coordinate control.
#
# Output: Figure3.pdf (vector, for submission), .png (600 dpi)
# ---------------------------------------------------------------------------

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
n <- function(v) suppressWarnings(as.numeric(ifelse(v == "", NA, v)))

obs  <- rd("observations.csv")
obs  <- obs[order(n(obs$day)), ]
day  <- n(obs$day)
pred <- n(obs$prednisolone_mg)
pr3  <- n(obs$anti_pr3_u_ml)
cd19 <- n(obs$cd19_cells_ul)

TX    <- rd("treatments.csv")
lanes <- unique(TX$agent)
TP    <- rd("timepoints.csv")

tec        <- n(TX$start_day[TX$agent == tail(lanes, 1)])
TEC_WINDOW <- c(min(tec), max(tec))

XLIM   <- c(-338, 208)
XTICKS <- seq(-300, 200, 60)

# ---- 2. style --------------------------------------------------------------
CO <- c(plex = "#9E6BA8", cyc = "#E1A33C", rtx = "#4C9F70", ava = "#8C8C8C",
        dara = "#3E8FC4", tec = "#C0392B", pr3 = "#B2182B", cd19 = "#0E7C7B",
        pred = "#2B5FA8", grey = "#595959", rail = "#EDEDED")
LANE_COLOUR <- setNames(c("plex", "cyc", "rtx", "ava", "dara", "tec"), lanes)

# Arial, as requested. macOS (quartz) and Windows (GDI) resolve it directly;
# on Linux fontconfig substitutes the metric-compatible Liberation Sans. Only a
# plain pdf() has to fall back to "sans", because it knows just its own
# PostScript families and aborts with "invalid font type" on anything else.
FAM <- if (isTRUE(suppressWarnings(capabilities("aqua"))) ||
           .Platform$OS.type == "windows" ||
           isTRUE(suppressWarnings(capabilities("cairo")))) "Arial" else "sans"

LW <- function(pt) pt * 96 / 72          # matplotlib points -> R lwd units
CX <- function(pt) pt / 7                # point size -> cex (base ps = 7)
MS <- function(pt) pt / 5.25              # marker diameter in pt -> cex
alpha <- function(col, a) {
  v <- col2rgb(col) / 255
  rgb(v[1], v[2], v[3], alpha = a)
}

# panel rectangles in device coordinates, matching the matplotlib gridspec
# (height_ratios 2.6/0.45/1.5/0.95/1.15, hspace 0.30, top 0.925, bottom 0.072)
ratios <- c(2.6, 0.45, 1.5, 0.95, 1.5)
Htot   <- 0.925 - 0.072
axtot  <- Htot / (1 + 4 * 0.30 / length(ratios))
gap    <- 0.30 * axtot / length(ratios)
hts    <- ratios / sum(ratios) * axtot
y0s    <- rev(cumsum(c(0.072, rev(hts + gap)[-length(hts)])))
PANEL  <- lapply(seq_along(hts), function(i) c(0.165, 0.975, y0s[i], y0s[i] + hts[i]))
names(PANEL) <- c("T", "U", "L", "B", "P")

panel <- function(id, ylim, first = FALSE) {
  par(fig = PANEL[[id]], mar = c(0, 0, 0, 0), new = !first)
  plot.new()
  plot.window(xlim = XLIM, ylim = ylim, xaxs = "i", yaxs = "i")
}

tick <- function(pt) -(pt / 72) / par("cin")[2]   # tick length in text lines

# Axis furniture is drawn explicitly rather than through axis()/mtext(), so
# that spines span the full data range (as matplotlib's do, whereas R's
# axis() line stops at the outermost tick) and labels land on the same
# coordinates in both scripts.
TICKLAB_X <- -346.9                      # right edge of the y tick labels

# points -> user units in the panel that is currently open
pt_x <- function(pt) pt / 72 * diff(par("usr")[1:2]) / par("pin")[1]
pt_y <- function(pt) pt / 72 * diff(par("usr")[3:4]) / par("pin")[2]

# xpd = NA: a spine sits exactly on the panel edge, so clipping to the panel
# would cut away half of its width
spine_l <- function(yl) segments(XLIM[1], yl[1], XLIM[1], yl[2],
                                 lwd = LW(0.7), xpd = NA)
spine_b <- function(y)  segments(XLIM[1], y, XLIM[2], y,
                                 lwd = LW(0.7), xpd = NA)

# text(): R anchors a string on its baseline, matplotlib on the bottom of its
# bounding box, hence the descender correction.
lab <- function(x, y, txt, ha = "left", va = "bottom", cex = 1,
                col = "black", srt = 0, font = 1) {
  fs   <- cex * par("ps")
  adjx <- c(left = 0, right = 1, centre = 0.5)[[ha]]
  d    <- if (va == "bottom") 0.212 * fs else 0
  if (srt == 0) {
    text(x, y + pt_y(d), txt, adj = c(adjx, if (va == "centre") 0.5 else 0),
         cex = cex, col = col, font = font, xpd = NA)
  } else {
    text(x + pt_x(d), y, txt, adj = c(0.5, if (va == "centre") 0.5 else 0),
         cex = cex, col = col, font = font, srt = srt, xpd = NA)
  }
}

# a multi-line block: matplotlib stacks lines 1.2 em apart and centres the
# whole block on the anchor; drawn line by line so the leading matches
labn <- function(x, y, lines, ha = "left", cex = 1, col = "black", srt = 0,
                 font = 1) {
  k <- seq_along(lines)
  o <- ((length(lines) - 1) / 2 - (k - 1)) * 1.2 * cex * par("ps")
  for (i in k)
    if (srt == 0)
      lab(x, y + pt_y(o[i]), lines[i], ha = ha, va = "centre", cex = cex,
          col = col, font = font)
    else
      lab(x - pt_x(o[i]), y, lines[i], va = "centre", cex = cex, col = col,
          font = font, srt = srt)
}

yticks <- function(at) {
  axis(2, at = at, labels = FALSE, lwd = 0, lwd.ticks = LW(0.7), tcl = tick(2.8))
  lab(TICKLAB_X, at, format(at), ha = "right", va = "centre")
}
ytitle <- function(txt, x, y)
  labn(x, y, txt, cex = CX(7.5), srt = 90)

# ---- 3. helpers ------------------------------------------------------------
blk <- function(x, lane, col, w = 5.5, h = 0.34)
  rect(x - w / 2, lane - h / 2, x + w / 2, lane + h / 2, col = col, border = NA)

bar <- function(x0, x1, lane, col, h = 0.34, a = 1)
  rect(x0, lane - h / 2, x1, lane + h / 2,
       col = if (a < 1) alpha(col, a) else col, border = NA)

lane_label <- function(x, lane, txt, ha = "left", col = CO["grey"], dx = 6)
  lab(x + if (ha == "left") dx else -dx, lane, txt, ha = ha, va = "centre",
      cex = CX(5.9), col = col)

timemarks <- function() {
  u <- par("usr")
  rect(TEC_WINDOW[1], u[3], TEC_WINDOW[2], u[4],
       col = alpha(CO["tec"], 0.055), border = NA)
  segments(0, u[3], 0, u[4], col = CO["tec"], lwd = LW(0.8), lty = "42")
}

pletter <- function(s, dy = 1) {
  u <- par("usr")
  lab(u[1] - 0.105 * (u[2] - u[1]), u[3] + dy * (u[4] - u[3]), s,
      cex = CX(9), font = 2)
}

step_series <- function(x, y) {
  ok <- !is.na(y); xs <- x[ok]; ys <- y[ok]
  ox <- xs[ys <  PULSE_CUTOFF]; oy <- ys[ys <  PULSE_CUTOFF]
  px <- xs[ys >= PULSE_CUTOFF]; py <- ys[ys >= PULSE_CUTOFF]
  for (p in px) {
    nxt <- ox[ox > p]
    if (length(nxt)) { v <- oy[ox == nxt[1]][1]; ox <- c(ox, p); oy <- c(oy, v) }
  }
  o <- order(ox)
  list(ox = ox[o], oy = oy[o], px = px, py = py)
}

step_xy <- function(x, y) {          # "post" steps, as polygon/line vertices
  k <- length(x)
  list(x = rep(x, each = 2)[-1], y = rep(y, each = 2)[-(2 * k)])
}

# ---- 4. figure -------------------------------------------------------------
draw <- function() {
  par(ps = 7, family = FAM, lend = "butt", xpd = FALSE)

  ## ---- A  treatment timeline ----------------------------------------------
  ymap <- setNames(rev(seq_along(lanes)) - 1, lanes)
  panel("T", c(-0.55, length(lanes) - 0.35), first = TRUE)

  for (l in lanes)
    segments(XLIM[1], ymap[l], XLIM[2], ymap[l], col = CO["rail"], lwd = LW(2.2))
  timemarks()

  # a dose block is never wider than the gap to the next dose of the same agent
  widths <- sapply(lanes, function(l) {
    d <- sort(n(TX$start_day[TX$agent == l & TX$end_day == ""]))
    g <- if (length(d) > 1) diff(d) else 99
    min(5.5, max(2.6, min(g) * 1.2))
  })

  for (i in seq_len(nrow(TX))) {
    r <- TX[i, ]; l <- r$agent; col <- CO[LANE_COLOUR[l]]
    x0 <- n(r$start_day)
    if (r$end_day != "") {
      x1  <- n(r$end_day)
      pad <- if (x1 - x0 < 30) 2.5 else 0   # pad short clusters, not courses
      bar(x0 - pad, x1 + pad, ymap[l], col, a = if (pad == 0) 0.9 else 1)
    } else {
      blk(x0, ymap[l], col, w = widths[l])
    }
    if (r$label == "") next
    if (r$label_align == "inside") {
      lab((x0 + n(r$end_day)) / 2, ymap[l], r$label, ha = "centre",
          va = "centre", cex = CX(5.9), col = "white")
    } else {
      lane_label(n(r$label_day), ymap[l], r$label, ha = r$label_align,
                 col = if (l == tail(lanes, 1)) CO["tec"] else CO["grey"])
    }
  }

  lab(-342.42, ymap, names(ymap), ha = "right", va = "centre")

  top <- length(lanes) - 0.42
  for (i in seq_len(nrow(TP))) {
    x <- n(TP$day[i])
    points(x, top, pch = 25, cex = MS(4.2), col = CO["grey"], bg = CO["grey"],
           xpd = NA)
    lab(x, top + 0.13, TP$label[i], ha = "centre", cex = CX(6.5), font = 2,
        col = CO["grey"])
  }
  lab(-342.42, top, TP$kind[1], ha = "right", va = "centre", cex = CX(6.5),
      col = CO["grey"])
  lab(4, top + 0.13, "Teclistamab, day 0", cex = CX(6.5), font = 2,
      col = CO["tec"])
  pletter("A", dy = 1.16)

  ## ---- B  anti-PR3 IgG (broken axis) --------------------------------------
  ok    <- !is.na(pr3)
  xp    <- day[ok]; yp <- pr3[ok]
  zero  <- yp == 0
  yplot <- ifelse(zero, PR3_LLOQ, yp)

  for (id in c("U", "L")) {
    panel(id, if (id == "U") c(142, 162) else c(-1.2, 31))
    if (id == "L") {
      rect(XLIM[1], -1.2, XLIM[2], PR3_ULN, col = alpha(CO["pr3"], 0.06),
           border = NA)
      segments(XLIM[1], PR3_ULN, XLIM[2], PR3_ULN, col = CO["pr3"],
               lwd = LW(0.6), lty = "22")
    }
    timemarks()
    lines(xp, yplot, col = CO["pr3"], lwd = LW(1.2))
    points(xp[!zero], yplot[!zero], pch = 21, cex = MS(3.2), bg = CO["pr3"],
           col = "white", lwd = LW(0.5))
    points(xp[zero], yplot[zero], pch = 21, cex = MS(3.2), bg = "white",
           col = CO["pr3"], lwd = LW(0.9))
    spine_l(if (id == "U") c(142, 162) else c(-1.2, 31))
    if (id == "L") spine_b(-1.2)
    yticks(if (id == "U") 150 else c(0, 10, 20, 30))
    u <- par("usr")
    if (id == "L") {
      lab(XLIM[1] + 6, PR3_ULN + 0.9, "normal < 2 U/mL", cex = CX(6),
          col = CO["pr3"])
      segments(u[1] - 0.007 * (u[2] - u[1]), u[3] + 0.94 * (u[4] - u[3]),
               u[1] + 0.007 * (u[2] - u[1]), u[3] + 1.14 * (u[4] - u[3]),
               lwd = LW(0.7), xpd = NA)
      ytitle("Anti-PR3 IgG (U/mL)", -406.37, 18.672)
    } else {
      segments(u[1] - 0.007 * (u[2] - u[1]), u[3] - 0.16 * (u[4] - u[3]),
               u[1] + 0.007 * (u[2] - u[1]), u[3] + 0.16 * (u[4] - u[3]),
               lwd = LW(0.7), xpd = NA)
      pletter("B", dy = 1.5)
    }
  }

  ## ---- C  CD19+ B cells ----------------------------------------------------
  panel("B", c(-1.8, 28))
  timemarks()
  ob <- !is.na(cd19)
  lines(day[ob], cd19[ob], col = CO["cd19"], lwd = LW(1.2))
  points(day[ob], cd19[ob], pch = 22, cex = MS(3.0), bg = CO["cd19"],
         col = "white", lwd = LW(0.5))
  spine_l(c(-1.8, 28)); spine_b(-1.8); yticks(c(0, 10, 20))
  ytitle(c("CD19+ B cells", "(cells/\u00B5L)"), -412.54, 13.091)
  pletter("C", dy = 1.06)

  ## ---- D  glucocorticoids --------------------------------------------------
  # the axis has to reach 100: the two oral 100 mg doses on rituximab days
  # were cut off by the earlier 60 mg limit
  YMAX <- 105; YARR <- 115.5

  s <- step_series(day, pred)
  panel("P", c(0, YMAX))
  timemarks()
  xs <- c(s$ox, tail(day, 1)); ys <- c(s$oy, tail(s$oy, 1))
  sp <- step_xy(xs, ys)
  polygon(c(sp$x[1], sp$x, tail(sp$x, 1)), c(0, sp$y, 0),
          col = alpha(CO["pred"], 0.18), border = NA)
  lines(sp$x, sp$y, col = CO["pred"], lwd = LW(1.2))
  points(s$ox, s$oy, pch = 21, cex = MS(2.6), bg = CO["pred"], col = "white",
         lwd = LW(0.4))

  # i.v. pulses as arrows in a strip above the axis; pulses given within
  # 10 days of each other are pooled into one arrow
  grp <- cumsum(c(TRUE, diff(s$px) > 10))
  for (g in split(seq_along(s$px), grp)) {
    gx <- mean(s$px[g])
    segments(gx, YARR - 7.9, gx, YARR, col = CO["pred"], lwd = LW(0.7),
             xpd = NA)
    points(gx, YARR - 9.6, pch = 25, cex = MS(4.2), col = CO["pred"],
           bg = CO["pred"], xpd = NA)
    lab(gx, YARR + 1.4, paste(as.integer(s$py[g]), collapse = ", "),
        ha = "centre", cex = CX(5.8), col = CO["pred"])
  }
  lab(XLIM[2] - 4, YARR + 1.4, "i.v. methylprednisolone pulses (mg)",
      ha = "right", cex = CX(6), col = CO["pred"])

  spine_l(c(0, YMAX)); spine_b(0); yticks(seq(0, 100, 25))
  axis(1, at = XTICKS, lwd = 0, lwd.ticks = LW(0.7), tcl = tick(2.8),
       mgp = c(3, 0.15, 0), cex.axis = 1)
  axis(1, at = seq(-330, 200, 30), labels = FALSE, lwd = 0,
       lwd.ticks = LW(0.6), tcl = tick(1.6))
  ytitle(c("Prednisolone", "(mg/day)"), -412.48, 42.134)
  lab(-66.48, -26.575, "Days relative to first teclistamab dose", ha = "centre",
      va = "centre", cex = CX(7.5))
  pletter("D", dy = 1.06)
}

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

open_pdf(file.path(OUT, "TEC_GPA_R.pdf"), 7.09, 7.35); draw(); invisible(dev.off())
open_png(file.path(OUT, "TEC_GPA_R.png"), 7.09, 7.35); draw(); invisible(dev.off())
cat("ok\n")
