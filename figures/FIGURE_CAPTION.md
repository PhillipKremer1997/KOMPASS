# Figure 1 — draft caption

**Treatment course and response to teclistamab in refractory PR3-ANCA-associated
vasculitis.** Day 0 is the first teclistamab dose; the shaded red band marks the
teclistamab exposure window (6 subcutaneous doses over 5 weeks).
**(a)** Immunosuppressive and targeted therapies. Bars denote continuous
exposure, lollipops individual infusions or injections; circled letters A–C mark
imaging timepoints. **(b)** Oral prednisolone maintenance dose. Intravenous
methylprednisolone pulses given as induction or infusion premedication are shown
in panel (a) rather than here, so that the maintenance taper is not obscured.
The dashed line marks the ≤5 mg/day low-dose glucocorticoid target.
**(c)** Anti-PR3 IgG on a logarithmic scale; the shaded band is the assay
negative range and open symbols denote values below the assay cut-off.
**(d)** CD19+ peripheral B cells.

## Points to verify before submission

- **Assay units** for anti-PR3 IgG (`UNIT_PR3`) and CD19+ cells (`UNIT_CD19`).
- **Anti-PR3 assay cut-off** (`PR3_ULN`, currently 2.0) — assay dependent.
- **Imaging timepoints** (`SCANS`). The boxes in the original figure were
  schematic rather than drawn to scale, so A/B/C were placed on the clinical
  narrative (first presentation, re-induction, post-teclistamab follow-up).
  Replace with the exact scan dates.
- **Teclistamab dose days** (`TEC_DOSE_DAYS`). The source record gives only the
  window (14 Jan – 19 Feb 2026) and the dose count, so only the window is drawn.
  Supplying the dates adds individual dose marks.
- **Classification of steroid doses as pulse vs. maintenance** (`PULSE_DAYS`) —
  see README.
