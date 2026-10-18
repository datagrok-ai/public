import { type ReadoutItem } from './readout-items.js';
import type { AnalysisPanel, AtomColours, AtomNumbers, KeymapProfile, MethylLabels, StereoDisplay, StereoPanel, UnspecifiedDoubleBonds, WarningMarks } from './types.js';
/** The settings, whole: what `settings` reads and writes and what the `settings` event carries. */
export interface SketcherSettings {
    readonly atomColours: AtomColours;
    readonly warningMarks: WarningMarks;
    readonly analysisPanel: AnalysisPanel;
    readonly stereoLabels: StereoDisplay;
    readonly stereoFlag: StereoDisplay;
    readonly unspecifiedDoubleBonds: UnspecifiedDoubleBonds;
    readonly methylLabels: MethylLabels;
    /** Each atom numbered from 1 on the canvas and in the API's SVG (spike stereo-commands): `'hidden'` (the default) or `'shown'`. */
    readonly atomNumbers: AtomNumbers;
    readonly keymap: KeymapProfile;
    /**
     * What the readout shows (spike readout-config; INFO-016): the items checked, by id, in the readout's order
     * (core/readout-items.ts); out of the box the formula and the two masses.
     */
    readonly readoutItems: readonly ReadoutItem[];
    /** The stereochemistry panel (spike stereo-panel): `'closed'` (the default), `'open'` or `'hidden'` (no button, no item). */
    readonly stereoPanel: StereoPanel;
}
/** The defaults: each option's own (Reset puts them back). */
export declare const DEFAULT_SETTINGS: SketcherSettings;
/** The names of the settings, in the dialog's order. */
export declare const SETTING_NAMES: readonly (keyof SketcherSettings)[];
/**
 * Settings as the sketcher takes them (written whole): each value by its option's own rule, what is left out or not
 * understood its default. So `{ atomColours: 'black' }` is the defaults with black atoms.
 */
export declare function toSettings(value: unknown): SketcherSettings;
/** Whether two settings say the same. */
export declare const sameSettings: (a: SketcherSettings, b: SketcherSettings) => boolean;
/** Whether two values of a setting say the same: a list of items by its ids (spike readout-config), any other by its value. */
export declare function sameSetting(name: keyof SketcherSettings, a: unknown, b: unknown): boolean;
