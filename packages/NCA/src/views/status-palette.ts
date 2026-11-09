/**
 * Status palette — the SINGLE source of truth for the three NCA verdict hues
 * (Accepted / Flagged / Not estimated) plus the neutral.
 *
 * **Radix "teal-success".**
 * All three are Radix Colors step 11 (one perceptual-lightness level, so the UI
 * does not "jump"), each ≥4.5:1 as text on white, and — the reason teal replaces
 * green — they stay distinct under deuteranopia/protanopia because teal keeps a
 * blue channel where green would merge with red/amber (report Key Findings #4/#5).
 *
 * **Why this file exists.** Every place NCA colour is applied — results-grid cell
 * TEXT, header-tooltip legend swatches, facet tile lines/markers/labels, the λz
 * editor header, and the cross-column fragility tint — must read from ONE palette
 * so the grid and the plots can never drift. results-view.ts binds its
 * `HUE_ACCEPTED/FLAGGED/NOT_ESTIMATED` to these constants and maps its band-tier
 * keys onto them via `fillToTextColour`; there is no second literal palette.
 *
 * **Format.** Datagrok colours are 0xAARRGGBB integers. `0xFF` = full opacity.
 * The low 6 nibbles are the Radix hex verbatim (e.g. teal-11 `#008573`).
 *
 * If a dark theme or a print/AAA export is added later, swap these values for the
 * report's dark step-11 / "dark for print" triads — this file is the only edit.
 */

/** Not estimated — Radix red-11 `#CE2C31` (~5.2:1 on white). */
export const STATUS_NOT_ESTIMATED = 0xFFCE2C31;

/** Flagged — Radix amber-11 `#AB6400` (~4.6:1 on white; borderline per Radix issue #42). */
export const STATUS_FLAGGED = 0xFFAB6400;

/** Accepted — Radix teal-11 `#008573` (~4.6:1 on white; CVD-safe replacement for green). */
export const STATUS_ACCEPTED = 0xFF008573;

/** Neutral — Radix slate-11 `#60646C`. For axes/de-emphasised text that must not read as a verdict. */
export const STATUS_NEUTRAL = 0xFF60646C;
