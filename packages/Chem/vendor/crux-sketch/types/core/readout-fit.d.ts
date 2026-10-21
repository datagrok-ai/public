import type { FormulaRun } from './analysis.js';
import type { BoxSize } from './density.js';
import type { AnalysisPanel } from './types.js';
/** The foot's inset from the canvas's sides (`.foot`: 6 px left, right and below). */
export declare const FOOT_INSET_PX = 6;
/**
 * What the readout's line takes besides its values (`.info`): its toggle, 24 px, the 6 px gap after it, 8 px of padding
 * at its end and a 1 px border each side.
 */
export declare const LINE_CHROME_PX = 40;
/**
 * Under this height, in CSS px, a canvas folds the readout whatever its width: its one line, 26 px over a 6 px inset,
 * would take more than a sixth of it. The QA's 220 x 198 box has a 128 px canvas and 140 x 126 a 56 px one; Datagrok's
 * dialog as it opens, 500 x 400, a 330 px one.
 */
export declare const FOLD_BELOW_HEIGHT_PX = 200;
/** The room the foot gives the readout's line in a canvas `width` px wide: the width less its two insets. */
export declare const footWidth: (width: number) => number;
/**
 * Whether a line `lineWidth` px wide fits the foot of a canvas `width` px wide, with a pixel to spare: the line is
 * measured to the whole pixel, and a line a fraction too wide would lose its last character to the ellipsis.
 */
export declare const lineFitsFoot: (lineWidth: number, width: number) => boolean;
/**
 * Whether the readout's one line, `lineWidth` px wide (its values and LINE_CHROME_PX), fits a canvas of `canvas` (R1):
 * the canvas is 200 px tall or more, and the line fits its foot. Where the readouts measured 252 to 359 px, Datagrok's
 * dialog as it opens, 500 x 400, gives them 398 px: they stay open there, as before. A canvas not measured yet, or with
 * no size (a hidden sketcher), fits it: the readout shows as it always did.
 */
export declare function readoutFits(lineWidth: number, canvas: BoxSize | null): boolean;
/** How the readout shows (R1, R2): not at all, folded to its toggle, or open. */
export type ReadoutLook = 'hidden' | 'folded' | 'open';
/**
 * The readout's look: hidden or collapsed as `analysisPanel` says; while it is `'shown'`, open where its line fits or
 * where the user opened it, else folded (R1). With nothing drawn (`empty`), hidden whatever the state, its toggle too
 * (spike ui-polish-oct12, row 4; INFO-009: the ui-critic's "lone chevron button" on an empty canvas).
 */
export declare function readoutLook(state: AnalysisPanel, fits: boolean, opened: boolean, empty?: boolean): ReadoutLook;
/** What a press on the readout's toggle leaves: the setting, and whether the user has the readout open where it does not fit. */
export interface ReadoutPress {
    readonly state: AnalysisPanel;
    readonly opened: boolean;
}
/**
 * A press on the readout's toggle (R2). Open, the readout closes: where its line fits, to `'collapsed'` as ever, a
 * setting its event says; where it does not, back to folded, a look. Closed, it opens: to `'shown'` (a setting only
 * from `'collapsed'`) and, where its line does not fit, open all the same, the user's. Hidden, it has no toggle.
 */
export declare function pressReadout(state: AnalysisPanel, fits: boolean, opened: boolean): ReadoutPress;
/**
 * The user's opening once the canvas has changed size (R2): kept while the box keeps its size class, the line fitting as
 * it did or not fitting as it did; dropped when that changes, so that a readout opened in a small box folds again in
 * the next small box after a box that fitted it.
 */
export declare const openedAfterResize: (opened: boolean, fitted: boolean, fits: boolean) => boolean;
/** A formula's runs as plain text, as a user reads them: superscripts as superscript characters, counts as plain digits. */
export declare const formulaPlain: (runs: readonly FormulaRun[]) => string;
