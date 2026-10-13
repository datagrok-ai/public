import type { AnalysisPanel, AtomColours, KeymapProfile, MethylLabels, StereoDisplay, UnspecifiedDoubleBonds, WarningMarks } from './types.js';
/** The settings, whole: what `settings` reads and writes and what the `settings` event carries. */
export interface SketcherSettings {
    readonly atomColours: AtomColours;
    readonly warningMarks: WarningMarks;
    readonly analysisPanel: AnalysisPanel;
    readonly stereoLabels: StereoDisplay;
    readonly stereoFlag: StereoDisplay;
    readonly unspecifiedDoubleBonds: UnspecifiedDoubleBonds;
    readonly methylLabels: MethylLabels;
    readonly keymap: KeymapProfile;
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
