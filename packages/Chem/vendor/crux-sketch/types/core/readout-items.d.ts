import type { FormulaRun, MassInfo } from './analysis.js';
import { formulaPlain } from './readout-fit.js';
import type { Catalog } from './strings.js';
/**
 * An item the readout can show, by its id (`readoutItems`; spike readout-config): the formula and masses, the notations, the
 * counts, the properties and the pointer's items.
 */
export type ReadoutItem = 'formula' | 'molecular-weight' | 'exact-mass' | 'smiles' | 'cxsmiles' | 'smarts' | 'inchikey' | 'atoms' | 'heavy-atoms' | 'bonds' | 'rings' | 'fragments' | 'stereocentres' | 'charge' | 'clogp' | 'tpsa' | 'hbd' | 'hba' | 'rotatable-bonds' | 'lipinski-violations' | 'qed' | 'under-pointer' | 'drawing';
/** A group of the readout's items, by its id: its menu's submenu and its Settings dialog's group. */
export type ReadoutGroup = 'formula-and-mass' | 'notations' | 'counts' | 'properties' | 'pointer';
/** The groups of the readout's items, in the menus' and the Settings dialog's order. */
export declare const READOUT_GROUPS: readonly {
    readonly id: ReadoutGroup;
    readonly items: readonly ReadoutItem[];
}[];
/** Every item, in the readout's order: the line shows the checked ones so. */
export declare const READOUT_ITEMS: readonly ReadoutItem[];
/** What the readout shows out of the box: the formula and the two masses, as before (INFO-001). */
export declare const DEFAULT_READOUT_ITEMS: readonly ReadoutItem[];
/** Whether a text is an item's id. */
export declare const isReadoutItem: (id: unknown) => id is ReadoutItem;
/** The group an item is in. */
export declare const groupOf: (item: ReadoutItem) => ReadoutGroup;
/**
 * `readoutItems` as the sketcher takes it (the option, the property, the `readout-items` attribute, a kept setting): a list
 * of ids, or a text of them apart by spaces or commas; each id it knows, once, in the readout's order; ids it does not
 * know are left out. Anything else (`undefined`, `null`, a number) is the default. An empty list or text is no item.
 */
export declare function toReadoutItems(value: unknown): readonly ReadoutItem[];
/** Whether two lists of items say the same. */
export declare const sameItems: (a: readonly string[], b: readonly string[]) => boolean;
/** An item checked or unchecked in a list (the menu's check boxes): the list, in the readout's order. */
export declare function withItem(items: readonly ReadoutItem[], item: ReadoutItem, checked: boolean): readonly ReadoutItem[];
/**
 * The engine's `descriptors` (crux-core sketch-arom-props), named as RDKit.js's `get_descriptors()`, plus `qed`; null where
 * a name is missing. `query` for a query, which has none.
 */
export interface Descriptors {
    readonly query: boolean;
    readonly CrippenClogP: number | null;
    readonly tpsa: number | null;
    readonly NumHBD: number | null;
    readonly NumHBA: number | null;
    readonly lipinskiHBD: number | null;
    readonly lipinskiHBA: number | null;
    readonly NumRotatableBonds: number | null;
    readonly NumHeavyAtoms: number | null;
    readonly NumAtoms: number | null;
    readonly NumRings: number | null;
    readonly NumAtomStereoCenters: number | null;
    readonly NumUnspecifiedAtomStereoCenters: number | null;
    readonly amw: number | null;
    readonly qed: number | null;
}
/** `descriptors`' JSON, read; null when it is no answer of that shape. */
export declare function parseDescriptors(json: string | null | undefined): Descriptors | null;
/**
 * Lipinski's rule of five broken (Lipinski 1997, as RDKit's `rdMolDescriptors` counts its parts): a molecular weight over
 * 500, a cLogP over 5, more than 5 H-bond donors and more than 10 acceptors (`lipinskiHBD`, `lipinskiHBA`). Null when the
 * engine gives none of them.
 */
export declare function lipinskiViolations(d: Descriptors): number | null;
/** What is under the pointer, or being drawn or moved (INFO-013), as the editor says it. */
export type PointerFacts = {
    readonly kind: 'atom';
    /** Its number, counted from 1 as the atom-number setting and the molfile count it. */
    readonly number: number;
    readonly element: string;
    readonly charge: number;
    readonly h: number;
    /** RDKit's total valence where the engine gives it; null where it does not. */
    readonly valence: number | null;
} | {
    readonly kind: 'bond';
    readonly number: number;
    /** The bond's type in the sketcher's words ("single", "aromatic"). */
    readonly type: string;
    /** Its length, in the drawing's coordinates (Å, as RDKit's depiction writes them). */
    readonly length: number;
} | {
    /** A bond drawn, or atoms moved: its length or the distance, and its direction, counterclockwise from the x axis. */
    readonly kind: 'draw' | 'move';
    readonly length: number;
    readonly angle: number;
};
/**
 * Where a line's values come from (core/sketcher.ts, the element's editor): each read on demand, cached there per
 * molecule version and selection, so that a line asks only for what its checked items need.
 */
export interface ReadoutSource {
    /** `massInfo` of the drawing or the selection; null when the engine cannot say. */
    mass(): MassInfo | null;
    /** A notation of the drawing or the selection: its text, `''` when it cannot be written, null while its module loads. */
    text(format: 'smiles' | 'cxsmiles' | 'smarts' | 'inchikey'): string | null;
    /**
     * `descriptors` of the drawing or the selection; null where the engine has none (before crux-core sketch-arom-props);
     * `'pending'` while they wait for an idle moment (spike readout-config's cost: 6 ms at 405 atoms).
     */
    descriptors(): Descriptors | null | 'pending';
    /** The bonds among the atoms described: the drawing's, or the selection's atoms' among themselves. */
    bonds(): number;
    /** What is under the pointer or being drawn; null for nothing. */
    pointer(): PointerFacts | null;
}
/** An item as the line shows it. */
export interface ShownItem {
    readonly id: ReadoutItem;
    /** Its label before its value in the line ("MW", "SMILES"); `''` for the formula, which has none. */
    readonly label: string;
    /** Its label's full words where the label is short ("Molecular weight" for "MW"): the label's tooltip. */
    readonly labelName?: string;
    /** Its name in the menus and in what the status line says ("Molecular weight", "SMILES"). */
    readonly name: string;
    /** Its value as text, as a user reads it (the formula's superscripts as superscript characters). */
    readonly value: string;
    /** The formula's runs (INFO-003): the readout sets its counts as subscripts and the rest as superscripts. */
    readonly runs?: readonly FormulaRun[];
    /** What a click on it copies: its value. */
    readonly copy: string;
    /** Its tooltip: a query's SMARTS says why there is no formula; a long notation, itself. */
    readonly title?: string;
    /** Its value is on its way (the InChI module loading): `…` shows until it is in. */
    readonly pending?: boolean;
    /** "R not counted" after it: the last mass shown, where R labels or attachment points count in the formula only (INFO-006). */
    readonly note?: string;
    /** A pointer's item, which the line keeps in sight before the others when there is no room for all (ruling 3). */
    readonly pointer?: boolean;
}
/** What the readout shows: nothing, or its items for the drawing or the selection. */
export type ReadoutLine = {
    readonly kind: 'empty';
} | {
    readonly kind: 'items';
    readonly selection: boolean;
    readonly items: readonly ShownItem[];
};
export declare const EMPTY_LINE: ReadoutLine;
/** A formula's runs as plain text, as a user reads them: superscripts as superscript characters, counts as plain digits. */
export declare const runsPlain: typeof formulaPlain;
/** A charge as the readout writes it: `0`, `+1`, `−2` (a minus sign as wide as the plus). */
export declare const chargeText: (n: number) => string;
/** A decimal as the readout writes it: a minus sign for a negative, none for one that rounds to zero (cLogP −0.0014 is 0.00). */
export declare function fixed(n: number, digits: number): string;
/** An angle in whole degrees, counterclockwise from the x axis, in (−180, 180]. */
export declare function angleText(deg: number): string;
/** Each item's name in the menus and the status line, and its label in the line. */
export declare function itemNames(s: Catalog): Readonly<Record<ReadoutItem, {
    readonly name: string;
    readonly label: string;
    readonly labelName?: string;
}>>;
/** A group's name in the menus and the Settings dialog. */
export declare function groupName(group: ReadoutGroup, s: Catalog): string;
/** The longest notation the line shows whole; a longer one is cut with an ellipsis there, its whole text its tooltip. */
export declare const NOTATION_TITLE_FROM = 24;
/**
 * The readout's line (rulings 1 and 2): the checked items that apply to the drawing or the selection (`selection`), in the
 * readout's order, each with its label, its value and what a click copies. Nothing for an empty drawing (INFO-009). A query
 * shows its SMARTS in place of the formula and masses (the query's message its tooltip), and its SMILES, CXSMILES,
 * InChIKey, counts and properties step aside; so does an item the engine gives no value for. Only what the checked items
 * need is asked of `source`.
 */
export declare function readoutLine(checked: readonly ReadoutItem[], source: ReadoutSource, selection: boolean, empty: boolean, s: Catalog): ReadoutLine;
/** What the pointer's item says: an atom's number, element, charge, H count and valence; a bond's type and length; a bond drawn or a move. */
export declare function pointerText(p: PointerFacts, s: Catalog): string;
/**
 * A line whose values on their way (`pending`: the descriptors and the InChIKey waiting for an idle moment) show what the
 * line before showed for them, so that an edit makes them change once, not flash `…` (spike readout-config's cost); still
 * pending, they copy nothing. An item the line before did not show stays `…`, and so does every one where the line
 * before was of the drawing and this one of a selection, or the other way round.
 */
export declare function withStale(line: ReadoutLine, before: ReadoutLine): ReadoutLine;
/** An item as the line reads it, and as Copy all writes it: its label, then its value, then its note. */
export declare const itemText: (i: ShownItem) => string;
/**
 * The line's whole text (Copy all, the tooltip of a line cut short): "Selection" before a selection's, each item's text;
 * an item still on its way with nothing to show (`…`) is left out.
 */
export declare function lineText(line: ReadoutLine, s: Catalog): string;
/**
 * The items a line of `room` px shows (ruling 3), the rest folded into the `…` button (`more` px with its gap): each item's
 * width with its gap; the pointer's items kept first, then the others in order while they fit. The shown items' indexes,
 * ascending; all of them when all fit; none when not even one fits beside the `…` button.
 */
export declare function fitItems(widths: readonly number[], pointer: readonly boolean[], room: number, more: number): number[];
