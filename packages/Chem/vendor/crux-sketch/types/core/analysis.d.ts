import type { AnalysisPanel, SketchAnalysis } from './types.js';
/** `analysisPanel` by HTML's rule for enumerated attributes (K12): `collapsed` or `hidden` in any case, `shown` for anything else. */
export declare const toAnalysisPanel: (value: unknown) => AnalysisPanel;
/** One fragment of `massInfo`'s answer: the connected pieces among the atoms it was asked about. */
export interface MassFragment {
    readonly atoms: readonly number[];
    readonly formula: string;
    readonly averageMass: number;
    readonly exactMass: number;
    readonly charge: number;
    readonly radicalElectrons: number;
    /** R labels and attachment points, each an `R` in the formula that weighs nothing (INFO-006). */
    readonly rSites: number;
    readonly dummies: number;
}
/** `massInfo`'s answer, read: a query has no numbers. */
export interface MassInfo {
    readonly query: boolean;
    readonly formula: string;
    readonly averageMass: number;
    readonly exactMass: number;
    readonly fragments: readonly MassFragment[];
}
/** `massInfo`'s JSON, read; null when it is no answer of that shape. */
export declare function parseMassInfo(json: string | null): MassInfo | null;
/** What the API gives of an answer (K13): the formula as the engine writes it, the totals and each fragment's; null for a query. */
export declare function toAnalysis(info: MassInfo | null): SketchAnalysis | null;
/** A run of the readout's formula: plain text, a subscript (a count) or a superscript (an isotope, a charge, a radical dot). */
export interface FormulaRun {
    readonly text: string;
    readonly kind: 'text' | 'sub' | 'sup';
}
/** A minus sign as wide as the plus (U+2212), as the readout writes a negative charge. */
export declare const MINUS = "\u2212";
/** A formula as the readout writes it (INFO-003): each fragment's runs, the fragments joined by dots (INFO-004). */
export declare function formulaRuns(formula: string): FormulaRun[];
/** The molecular weight as the readout shows it: 2 decimals (K12, provisional). */
export declare const weightText: (n: number) => string;
/** The exact mass as the readout shows it: 4 decimals (K12, provisional). */
export declare const exactMassText: (n: number) => string;
