import type { SketcherSettings } from './settings.js';
import type { SketcherStrings } from './strings.js';
export type { SketcherSettings } from './settings.js';
export type { Message, PluralForms, SketcherStrings } from './strings.js';
/**
 * `molecule` (the default) or `query` (API-021; spike query, Q2): query mode offers the means to draw query features (the
 * query tools, menu items and dialog sections, spike 11b); molecule mode hides them. Either mode reads, draws, keeps and
 * writes a value's query features, and the getters are the same in both.
 */
export type SketcherMode = 'molecule' | 'query';
/** Where a change came from: a user edit, or the host application through the API. */
export type ChangeSource = 'user' | 'api';
/** What a `change` listener receives, once per completed edit. */
export interface ChangeDetail {
    /** `'user'` for a user's edit (its undo and redo by key or button included), `'api'` for the host's write or `undo()`. */
    readonly source: ChangeSource;
}
/**
 * Why a write was refused: a plain object, never an `Error` (API-013). `format` is the format the
 * text was read as (`smiles`, `cxsmiles`, `smarts`, `molV2000`, `molV3000`, `inchi`, `cdxml`,
 * `unknown`), `reaction` for a reaction, which the sketcher never reads (spike formats, F12, IO-031: a reaction
 * SMILES or SMARTS, an RXN file, a V3000 reaction block; "Reactions are not supported"), or `engine` when the engine
 * itself failed: it could not be loaded, or it trapped, and then no sketcher on the page can read or draw anything more.
 */
export interface SketchError {
    readonly format: string;
    readonly message: string;
    /** 0-based character offset of a SMILES or SMARTS error, where the parser knows it. */
    readonly position?: number;
    /** 1-based line of a MOL block error. */
    readonly line?: number;
    /** The atoms a chemistry (sanitization) error is about. */
    readonly atoms?: readonly number[];
    /**
     * Why a format cannot be written or read now, where a host can act on it (milestone-4-fixes, F6): `loading`, its module
     * is still loading (`await loadFormat(format)` and read again); `unavailable`, its module could not be loaded, which
     * `message` says why, until `loadFormat(format)` tries again. Absent from every other error.
     */
    readonly code?: SketchErrorCode;
}
/** A `SketchError`'s code (milestone-4-fixes, F6; spike inchi): the CDXML or InChI module loading, or unavailable after a failed load. */
export type SketchErrorCode = 'loading' | 'unavailable';
/**
 * What a warning is about (spike charges-checks, K11): an atom with a valence RDKit's sanitization rejects (`valence`),
 * aromatic bonds with no Kekulé structure (`aromaticity`), a label kept as an alias, which nothing interprets (`label`),
 * a stereocentre nothing specifies in a structure drawn with stereo bonds (`stereo`, spike stereo, S13), or something a
 * getter's format cannot hold (`format`).
 */
export type SketchWarningKind = 'valence' | 'aromaticity' | 'label' | 'abbreviation' | 'stereo' | 'format';
/**
 * The loss table's rows (spike formats, F2, F11): what a format cannot hold, as a `format` warning's `code` names it.
 * `v3000` is no loss but a notice: under `molV2000`, that the molfile getter (and `value` after an edit) writes V3000, and
 * why (F13, IO-014). What a read drops (`readWarnings`, F10) has codes of its own: `sd-data` (an SD record's data
 * fields), `records` (an SD file's after the first), `rgroup-definitions` (an RG file's or a V3000 block's),
 * `coordinates-3d` (laid out afresh in 2D), `sgroup`, `query`, `cx-section` and `molfile` (what the reader cannot keep).
 * A CDXML read's (spike cdxml, X2, X3): `text`, `graphic`, `arrow`, `picture`, `bracket` (left out), `structure` (a
 * structure RDKit cannot read, left out) and `alias` (a label kept as an alias). CDXML's losses (X6) add `atom-maps` and
 * `stereo` (an allene's or a non-tetrahedral centre's). InChI's (spike inchi, N6, under `inchi`) add `r-label`, `dummy` and
 * `dative`, which with `query`, `hydrogen-bond` and `zero-order` stop the InChI (`""`), and the InChI library's own warnings:
 * `undefined-stereo`, `charges-rearranged`, `unusual-valence`, `ambiguous-stereo` and `library` (any other, in its words).
 */
export type SketchLossCode = 'r-label' | 'dummy' | 'dative' | 'undefined-stereo' | 'charges-rearranged' | 'unusual-valence' | 'ambiguous-stereo' | 'library' | 'alias' | 'r-label-numbers' | 'attachment-points' | 'stereo-groups' | 'chiral-flag' | 'query' | 'sgroup' | 'hydrogen-bond' | 'zero-order' | 'atom-maps' | 'stereo' | 'v3000';
/**
 * A problem of the molecule as it is, or what a getter's format cannot hold of it (decision `rgroups`, RG2, RG3,
 * RGROUP-018; spike charges-checks, K11, CHECK-012): a plain object in `SketchError`'s shape. Nothing is refused or thrown
 * for it: the molecule loads, edits, copies and writes as it is (CHECK-003).
 * - `kind`: `valence`, `aromaticity`, `label`, `abbreviation` (an abbreviation with more bonds than attachments, spike
 *   abbreviations, K1) and `stereo` (an unspecified stereocentre in a structure drawn with stereo bonds, spike stereo, S13)
 *   are the molecule's problems, the atoms the canvas marks; `format` is what a format loses.
 * - `format`: only on a `format` warning, the format that loses something (`smiles`, `cxsmiles`, `smarts`; the molfiles,
 *   `molV2000` and `molV3000`, hold R labels and attachment points whole; `cdxml` and `inchi` once their modules are in). Under `smiles`, what the `smiles` getter's
 *   text loses: when it is CXSMILES (decision `smiles-cx`), only what CXSMILES loses too, so attachment points warn under
 *   `smarts` alone.
 * - `message`: what the status line says of the problem, in ChemDraw's words ("Too many bonds to this unlabeled carbon"), or
 *   what the format loses.
 * - `atoms`: every atom concerned, by index: the over-valent atom, the aromatic atoms with no Kekulé structure, the atom
 *   whose label is kept, the unspecified stereocentre, or a combined R label (R1,R2) written with its lowest number,
 *   attachment points not written.
 * - `code` and `bonds` (spike formats, F11): only on a `format` warning, the loss table's row (`SketchLossCode`, or a read's
 *   code in `readWarnings`) and the bonds concerned (a dative, hydrogen or zero-order bond), `[]` when none.
 */
export interface SketchWarning {
    readonly kind: SketchWarningKind;
    readonly format?: string;
    readonly code?: SketchLossCode | (string & {});
    readonly message: string;
    readonly atoms: readonly number[];
    readonly bonds?: readonly number[];
}
/** One fragment's formula and masses (spike charges-checks, K13; INFO-004): the connected pieces of the atoms analysed. */
export interface SketchFragmentAnalysis {
    /** Its atoms, by index. */
    readonly atoms: readonly number[];
    readonly formula: string;
    readonly molecularWeight: number;
    readonly exactMass: number;
}
/**
 * The formula and masses of the drawing, or of the selected atoms (spike charges-checks, K13; INFO-010). `formula` is
 * RDKit's Hill formula as the engine writes it, each fragment's joined by dots: an isotope `[13C]`, D and T, the charge
 * last (`C2H3O2-.Na+`, `[13C]H4O`); an R label or an attachment point an `R` (`C6H5R`), a radical electron a `•` before
 * the charge (`C2H5•`). `molecularWeight` and `exactMass` are the totals, RDKit's average and exact masses, unrounded; an
 * R weighs nothing (INFO-006).
 */
export interface SketchAnalysis {
    readonly formula: string;
    readonly molecularWeight: number;
    readonly exactMass: number;
    readonly fragments: readonly SketchFragmentAnalysis[];
}
/**
 * Whether the canvas marks the atoms with problems (spike charges-checks, K10, CHECK-011): `'shown'`, the default, or
 * `'hidden'`, which hides the marks only; `warnings` and the status line keep saying them. Any other value is
 * `'shown'`, and reads back so.
 */
export type WarningMarks = 'shown' | 'hidden';
/**
 * The formula and mass readout at the canvas's foot (spike charges-checks, K12; UI-030): `'shown'`, the default,
 * `'collapsed'`, its toggle alone, or `'hidden'`, none at all. Any other value is `'shown'`, and reads back so.
 */
export type AnalysisPanel = 'shown' | 'collapsed' | 'hidden';
/**
 * What a `render` listener receives (API-061; spike datagrok-platform, R5): the canvas's render count once this drawing is
 * shown, as `renderCount` then reads it.
 */
export interface RenderDetail {
    readonly renderCount: number;
}
/**
 * The events a sketcher emits, and what each listener receives. `settings` (spike settings-touch, K6, P9): the settings in
 * force once a user's Apply in the Settings dialog changed them, so that a host may store them; never for a host's write.
 * `render` (API-061; spike datagrok-platform, R5): the canvas shows a new drawing (the first one included, then after an
 * edit, a write, a view change or a resize), once per drawing, as `renderCount` counts them, with the drawing on the page
 * when it comes. What the canvas marks on a drawing it shows (a tool's preview, the hover, the selection) is no render. A
 * host that waits for the sketcher to settle waits on it and on `isPending` (Datagrok's settles do), never on a timeout.
 */
export interface SketcherEvents {
    change: ChangeDetail;
    error: SketchError;
    settings: SketcherSettings;
    render: RenderDetail;
}
/** An event's name, as `on()` and `off()` take it: `change`, `error`, `settings` or `render`. */
export type SketcherEventType = keyof SketcherEvents;
/** A listener of the event `T`, called with that event's detail (`SketcherEvents[T]`). */
export type SketcherListener<T extends SketcherEventType> = (detail: SketcherEvents[T]) => void;
/** A format `setValue` can force; any other text is detected from its shape. */
export type SketcherFormat = 'smiles' | 'cxsmiles' | 'smarts' | 'molV2000' | 'molV3000' | 'inchi' | 'cdxml';
/** How `setValue` reads its text. */
export interface SetValueOptions {
    /** Read the text as this format instead of detecting it (e.g. SMARTS that is also SMILES). */
    format?: SketcherFormat;
}
/**
 * An action a toolbar can hold: it does something once, and is no tool that stays chosen. `element.table`
 * opens the periodic table; `zoom.in`, `zoom.out`, `zoom.reset` (100%) and `zoom.fit` (fit to view)
 * change the view (spike view); `layout` lays the whole drawing out afresh and fits the view to it, and
 * `clean` (Clean Up) lays it out afresh where it is, turned onto 15° (spike layout); `cut`, `copy` and `paste`
 * go through the system clipboard, as Ctrl+X, Ctrl+C and Ctrl+V do (spike clipboard). `structure.enhanced-stereo` opens
 * the Enhanced Stereochemistry dialog for the selected stereocentres, or every one with nothing selected (spike stereo).
 */
export type ToolbarAction = 'clear' | 'undo' | 'redo' | 'element.table' | 'zoom.in' | 'zoom.out' | 'zoom.reset' | 'zoom.fit' | 'layout' | 'clean' | 'cut' | 'copy' | 'paste' | 'structure.enhanced-stereo' | 'keys' | 'settings';
/** A readout a toolbar can hold: it shows a state and takes no input. `zoom.level` shows the zoom as a percentage (spike view). */
export type ToolbarReadout = 'zoom.level';
/**
 * A menu button a toolbar can hold (spike bar-fit, B7, B9): `zoom`, the zoom level as Ketcher shows it, whose menu holds Zoom
 * in, Zoom out, Actual size and Fit to view (it and `zoom.level` show the same level, and a configuration shows the first
 * of them it names); `options`, the gear, whose menu holds Settings…, Query mode (spike reports-oct5, R3) and Keyboard
 * shortcuts….
 */
export type ToolbarMenu = 'zoom' | 'options';
/**
 * An item of a menu that is no control of a bar (spike reports-oct5, R3), which a configuration's `hidden` takes away as it
 * takes a control: `mode`, the gear's Query mode, a check item that switches `mode` between `'molecule'` and `'query'`.
 */
export type ToolbarMenuItem = 'mode';
/** The 118 element symbols, as written. */
export type ElementSymbol = 'H' | 'He' | 'Li' | 'Be' | 'B' | 'C' | 'N' | 'O' | 'F' | 'Ne' | 'Na' | 'Mg' | 'Al' | 'Si' | 'P' | 'S' | 'Cl' | 'Ar' | 'K' | 'Ca' | 'Sc' | 'Ti' | 'V' | 'Cr' | 'Mn' | 'Fe' | 'Co' | 'Ni' | 'Cu' | 'Zn' | 'Ga' | 'Ge' | 'As' | 'Se' | 'Br' | 'Kr' | 'Rb' | 'Sr' | 'Y' | 'Zr' | 'Nb' | 'Mo' | 'Tc' | 'Ru' | 'Rh' | 'Pd' | 'Ag' | 'Cd' | 'In' | 'Sn' | 'Sb' | 'Te' | 'I' | 'Xe' | 'Cs' | 'Ba' | 'La' | 'Ce' | 'Pr' | 'Nd' | 'Pm' | 'Sm' | 'Eu' | 'Gd' | 'Tb' | 'Dy' | 'Ho' | 'Er' | 'Tm' | 'Yb' | 'Lu' | 'Hf' | 'Ta' | 'W' | 'Re' | 'Os' | 'Ir' | 'Pt' | 'Au' | 'Hg' | 'Tl' | 'Pb' | 'Bi' | 'Po' | 'At' | 'Rn' | 'Fr' | 'Ra' | 'Ac' | 'Th' | 'Pa' | 'U' | 'Np' | 'Pu' | 'Am' | 'Cm' | 'Bk' | 'Cf' | 'Es' | 'Fm' | 'Md' | 'No' | 'Lr' | 'Rf' | 'Db' | 'Sg' | 'Bh' | 'Hs' | 'Mt' | 'Ds' | 'Rg' | 'Cn' | 'Nh' | 'Fl' | 'Mc' | 'Lv' | 'Ts' | 'Og';
/** An element tool, named by its symbol in lower case: `element.n`, `element.cl` (spike atoms). */
export type ToolbarElementTool = `element.${Lowercase<ElementSymbol>}`;
/** The rings of the ring tools (spike rings-chains, RING-001), in the ring bar's order. */
export type RingName = 'benzene' | 'cyclopentadiene' | 'cyclohexane' | 'cyclopentane' | 'cyclopropane' | 'cyclobutane' | 'cycloheptane' | 'cyclooctane';
/** A ring tool, named by its ring: `ring.benzene`, `ring.cyclohexane` (spike rings-chains). */
export type ToolbarRingTool = `ring.${RingName}`;
/**
 * A control a toolbar can hold, named by its test id without `toolbar.`: the tools, the element tools
 * (`element.n`), the ring tools (`ring.benzene`), the chain tool (`chain`), the hand tool (`hand`), the
 * Structure Library (`template.library`: it opens the library, and shows pressed while a template from
 * it is the tool), the actions `clear` (clear canvas), `undo`, `redo`, `element.table` (the periodic
 * table) and the zoom actions, and the readout `zoom.level`; the R-group label tool (`rgroup.label`) and the
 * attachment-point tool (`rgroup.attachment`, spike rgroups); the selection tools, the rectangle
 * (`select.rect`), the lasso (`select.lasso`) and the structure tool (`select.structure`, spike select); the S-Group tool
 * (`sgroup`, spike sgroups).
 * `bond.wavy` and `bond.crossed` are the wavy bond and the crossed double bond (spike stereo). The query bond tools,
 * `bond.any`, `bond.single-or-double`, `bond.single-or-aromatic` and `bond.double-or-aromatic`, show in query mode only
 * (spike query-edit, API-021). A misspelt name fails to compile.
 */
export type ToolbarTool = ToolbarAction | ToolbarReadout | ToolbarMenu | 'select.rect' | 'select.lasso' | 'select.structure' | 'hand' | 'erase' | 'bond.single' | 'bond.double' | 'bond.triple' | 'bond.aromatic' | 'bond.wedge' | 'bond.hash' | 'bond.wavy' | 'bond.crossed' | 'bond.any' | 'bond.single-or-double' | 'bond.single-or-aromatic' | 'bond.double-or-aromatic' | 'chain' | 'template.library' | 'rgroup.label' | 'rgroup.attachment' | 'sgroup' | 'charge.plus' | 'charge.minus' | ToolbarRingTool | ToolbarElementTool;
/**
 * A sub-palette of a toolbar: its id (a test-id path of your own, such as `bond.multiple`: its main
 * button is `toolbar.<palette>`), its name, and its tools. Its main button shows the tool used last.
 */
export interface ToolbarPalette {
    readonly palette: string;
    readonly label: string;
    readonly tools: readonly Exclude<ToolbarTool, ToolbarAction | ToolbarReadout>[];
}
/** A control of a group: a tool or an action, or a sub-palette. */
export type ToolbarEntry = ToolbarTool | ToolbarPalette;
/** A group of a toolbar: a labelled `role="toolbar"` (A11Y-003). */
export interface ToolbarGroup {
    readonly label: string;
    readonly tools: readonly ToolbarEntry[];
}
/**
 * What the toolbars around the canvas hold (UI-001, UI-002): the top one, the left one, the right one
 * (the element palette, spike atoms) and the bottom one (the ring bar, spike rings-chains), each a
 * list of groups. A toolbar that is left out, or has no group, is not shown: `{}` shows no toolbars at
 * all, and the canvas takes the whole box. The elements picked in the periodic table join the group
 * that holds `element.table`, before it. Read-only (`as const`) configurations are accepted.
 *
 * The configuration is presentation, not permission: it decides which buttons show. The keys that
 * choose and use tools (j, Shift+J, Shift+X, the number keys and the ring keys) do not read it, so a
 * host that leaves out the ring bar or the chain button still gets rings and chains from the keyboard.
 */
export interface ToolbarConfig {
    readonly top?: readonly ToolbarGroup[];
    readonly left?: readonly ToolbarGroup[];
    readonly right?: readonly ToolbarGroup[];
    readonly bottom?: readonly ToolbarGroup[];
    /**
     * Controls hidden (spike settings-touch, K7, API-027): taken out of every bar that holds them, and their keys turned off
     * with them, the keys that choose a hidden tool or run a hidden action (Ctrl+Z with `undo` hidden, j with
     * `ring.benzene`). A configuration with `hidden` and no bar is the default bars less what it hides. A control merely left
     * out of the bars keeps its keys: the configuration is presentation, `hidden` is the host's word that a control goes.
     */
    readonly hidden?: readonly (ToolbarTool | ToolbarMenuItem)[];
}
/**
 * How atoms are coloured (H36, CANVAS-031, CANVAS-032): `'element'` (the default) draws heteroatoms and
 * the bond halves next to them in their element's colour; `'black'` draws everything black on white.
 * Any other value draws the default, and reads back as `'element'`.
 */
export type AtomColours = 'element' | 'black';
/**
 * The keymap's profile (H32, API-026): `'chemdraw'` (the default) reads keys as ChemDraw does; `'ketcher'` adds Ketcher's
 * keys where they do not clash (so far Ctrl+Shift+M, Copy As MOL). Any other value is ChemDraw's, and reads back so.
 */
export type KeymapProfile = 'chemdraw' | 'ketcher';
/**
 * Whether a part of the drawing shows (spike stereo, S8, S10): `stereoLabels`, the CIP labels ((R), (S), (E), (Z)),
 * `'hidden'` by default; `stereoFlag`, the drawing's flag (Absolute, Racemic, Relative), `'shown'` by default.
 */
export type StereoDisplay = 'shown' | 'hidden';
/**
 * How a value's double bond with no E/Z of its own (a SMILES' `CC=CC`) is drawn (spike stereo, P2, the product owner's of
 * 2026-10-03): `'crossed'` (the default), as ChemDraw draws the molfile RDKit writes for it, or `'plain'`. The getters write
 * it unspecified either way. A double bond the user draws takes its E/Z from its drawing, and a file's "either" bond and a
 * bond the user crossed stay crossed, whatever the option.
 */
export type UnspecifiedDoubleBonds = 'crossed' | 'plain';
/**
 * Whether terminal carbons are labelled CH₃ (spike settings-touch, K6, P10; UI-021's hydrogen labels): `'hidden'` (the
 * default, ChemDraw's and the ACS look's) or `'shown'`, moldraw2d's `explicitMethyl`. Hydrogens on heteroatoms are always
 * written (H23).
 */
export type MethylLabels = 'shown' | 'hidden';
/**
 * How large the toolbars' controls are (spike bar-fit, B1; UI-005): `'auto'` (the default) follows the sketcher's own box,
 * compact (28 px buttons, 20 px icons, a little less padding) under 720 px wide or 600 px tall, such as Datagrok's sketcher
 * dialog, and regular (32 px buttons, 24 px icons) otherwise; `'regular'` and `'compact'` keep that size in any box. No
 * target is ever under 24 px. It changes nothing but sizes: no change event, no undo step, no getter changes. Any other
 * value is `'auto'`, and reads back so.
 */
export type Density = 'auto' | 'regular' | 'compact';
/**
 * What `createSketcher` takes (H29): the initial molecule and mode, the toolbars, and the options the handle reads and
 * changes later by the same names. Every option is optional; left out, it takes its default.
 */
export interface CreateSketcherOptions {
    /** The initial molecule, in any format the sketcher reads. Showing it emits no change event. */
    value?: string;
    /** `'molecule'` (the default) or `'query'` (API-021): the handle's `mode` switches it later. */
    mode?: SketcherMode;
    /**
     * What the toolbars hold (UI-002): `{ top, left, right, bottom }`, each a list of groups `{ label,
     * tools }`, a tool named by its test id without `toolbar.` (`bond.single`, `element.n`,
     * `ring.benzene`), a sub-palette `{ palette, label, tools }`. Left out, the default: Undo and Redo,
     * then clear canvas, then the zoom group (zoom out, the zoom level, zoom in, 100%, fit to view), then
     * Layout and Clean Up, at the top; the selection palette (the rectangle, the lasso and the structure tool),
     * then the hand tool and the eraser, the bond tools, the stereo bonds palette and the chain tool, then the
     * R-group palette (the R-group label and attachment-point tools) on the left; H, C, N, O, S, P, F, Cl, Br, I and the periodic table on the right; the eight ring tools and the
     * Structure Library at the bottom. `{}` shows no toolbars, and the canvas takes the whole box.
     */
    toolbar?: ToolbarConfig;
    /** How atoms are coloured: `'element'` (the default) or `'black'`. The handle's `atomColours` changes it later. */
    atomColours?: AtomColours;
    /** The keymap's profile: `'chemdraw'` (the default) or `'ketcher'`. The handle's `keymap` changes it later. */
    keymap?: KeymapProfile;
    /** Whether the canvas marks the atoms with problems: `'shown'` (the default) or `'hidden'`. The handle's `warningMarks` changes it later. */
    warningMarks?: WarningMarks;
    /** The formula and mass readout: `'shown'` (the default), `'collapsed'` or `'hidden'`. The handle's `analysisPanel` changes it later. */
    analysisPanel?: AnalysisPanel;
    /** The CIP labels: `'hidden'` (the default) or `'shown'` (spike stereo, S10). The handle's `stereoLabels` changes it later. */
    stereoLabels?: StereoDisplay;
    /** The stereo flag: `'shown'` (the default) or `'hidden'` (spike stereo, S8). The handle's `stereoFlag` changes it later. */
    stereoFlag?: StereoDisplay;
    /** A value's unspecified double bonds: `'crossed'` (the default) or `'plain'` (spike stereo, P2). The handle's `unspecifiedDoubleBonds` changes it later. */
    unspecifiedDoubleBonds?: UnspecifiedDoubleBonds;
    /**
     * The locale the sketcher speaks (I18N-002): a BCP 47 tag, `'en'` by default; `'en-XA'` is the built-in pseudo-locale
     * (I18N-003). The handle's `locale` changes it later.
     */
    locale?: string;
    /** Translations of the sketcher's messages, partial: what they leave out is English (I18N-002). The handle's `strings` changes them later. */
    strings?: SketcherStrings;
    /** Terminal carbons labelled CH₃: `'hidden'` (the default) or `'shown'` (spike settings-touch, K6). The handle's `methylLabels` changes it later. */
    methylLabels?: MethylLabels;
    /**
     * The host's keys over the keymap's profile (spike settings-touch, K7, KEY-004, P11): `{ "Ctrl+L": "layout", "t": null }`,
     * a key bound to a built-in action or to nothing. The handle's `keys` changes them later.
     */
    keys?: HostKeys;
    /**
     * All the settings at once (spike settings-touch, K6), as `settings` writes them; applied after the options one by one.
     * Each setting it names, as each option given, is the host's own: it keeps its value over the user's kept one (spike
     * reports-oct5, R5).
     */
    settings?: Partial<SketcherSettings>;
    /**
     * Whether the `cxsmiles` getter writes the coordinates and the wedges as drawn (spike formats, F6): `false` (the default)
     * or `true`. The handle's `cxsmilesCoordinates` changes it later.
     */
    cxsmilesCoordinates?: boolean;
    /**
     * How large the toolbars' controls are (spike bar-fit, B1): `'auto'` (the default, by the sketcher's box), `'regular'` or
     * `'compact'`. The `density` attribute of `<crux-sketch>` sets it too; the handle's `density` changes it later.
     */
    density?: Density;
    /**
     * Whether the user's settings are kept and shared (spike reports-oct5, R5): `true` (the default) or `false`. On, a setting
     * the user changes (the Settings dialog, the readout's toggle, the canvas menu's Show stereochemistry) is kept in the
     * page's localStorage, under `crux-sketch:settings`, and reaches every other sketcher open on the page (and the page's
     * other windows) at once and every sketcher created later before its first drawing, but for the settings its host set
     * itself, by an option, a property or an attribute. Where the page has no storage, or it refuses, the sketcher works as
     * if off, saying nothing. Off, it neither reads nor writes the storage, nor hears nor tells another sketcher. The
     * `persist-settings` attribute of `<crux-sketch>` sets it too (`"false"` turns it off); the handle's changes it later.
     */
    persistSettings?: boolean;
    /**
     * The query elements the host turns off (API-028; spike query-edit): generic symbols (`'A'`, `'AH'`, `'Q'`, `'QH'`,
     * `'X'`, `'XH'`, `'M'`, `'MH'`, `'*'`), `'list'`, `'not-list'` and the query bonds (`'any'`, `'single-or-double'`,
     * `'single-or-aromatic'`, `'double-or-aromatic'`). The handle's `disabledQueryElements` changes them later.
     */
    disabledQueryElements?: readonly string[];
}
/** A point of the canvas as `positions` gives it: CSS px from the top-left corner of the `<crux-sketch>` element. */
export interface DrawnPoint {
    readonly x: number;
    readonly y: number;
}
/**
 * Where the canvas draws the drawing (API-056): each atom's position and each bond's midpoint, by index (the order of the
 * `molfile`'s atoms and bonds, as the test ids `atom.N` and `bond.N` number them), where a press lands on that atom or
 * bond. An atom drawn nowhere (one a contracted abbreviation hides) and a bond drawn nowhere are null. `bondLength` is a
 * bond's drawn length at the current zoom, in px; 0 when there is nothing drawn.
 */
export interface DrawnPositions {
    readonly atoms: readonly (DrawnPoint | null)[];
    readonly bonds: readonly (DrawnPoint | null)[];
    readonly bondLength: number;
}
/** A box of the sketcher as `layout` gives it: CSS px from the top-left corner of the `<crux-sketch>` element. */
export interface DrawnBox {
    readonly x: number;
    readonly y: number;
    readonly width: number;
    readonly height: number;
}
/** A toolbar's place around the canvas: the top bar, the left bar, the element palette on the right, the ring bar below. */
export type ToolbarPlace = 'top' | 'left' | 'right' | 'bottom';
/**
 * Where the sketcher shows its parts now (API-058; spike datagrok-platform, R3, R4), for a host that points at them, as
 * Datagrok's hit areas do: the canvas; each toolbar shown, by its place; each control a toolbar shows, by its name in a
 * toolbar configuration (`bond.single`, `ring.benzene`, `element.n`, `undo`, `zoom`; a sub-palette's button by the
 * palette's name, `select`, `bond.stereo`), where a press lands on it; and the label editor while it is open. What is not
 * shown is left out: a toolbar the configuration leaves out, a control in a closed palette or a folded group, one its bar
 * has scrolled out of view.
 */
export interface SketcherLayout {
    readonly canvas: DrawnBox | null;
    readonly toolbars: Readonly<Partial<Record<ToolbarPlace, DrawnBox>>>;
    readonly controls: Readonly<Record<string, DrawnBox>>;
    readonly labelEditor: DrawnBox | null;
}
/** The atoms and bonds the user has selected (API-059), by index (the `molfile`'s order, as `positions` numbers them), ascending. */
export interface SketchSelection {
    readonly atoms: readonly number[];
    readonly bonds: readonly number[];
}
/**
 * A tool the sketcher can have chosen (API-060), named as a toolbar configuration names it (`select.rect`, `bond.single`,
 * `ring.benzene`, `element.n`, `chain`, `erase`, `template.library` while a template from the Structure Library is the
 * tool), or `element.query`, the query atom picked last in the periodic table (query mode).
 */
export type SketcherTool = Exclude<ToolbarTool, ToolbarAction | ToolbarReadout | ToolbarMenu> | 'element.query';
/**
 * What `createSketcher` resolves to; `<crux-sketch>` has the same members.
 *
 * The format getters (`smiles`, `molfile`, `molV3000`, `cxsmiles`, `smarts`, `cdxml`) and `plainSmiles` (and `value`
 * once the user has edited) are `""` for an empty drawing; `smiles`, `plainSmiles` and `cxsmiles` also for a drawing with
 * a query feature SMILES cannot hold (spike query, Q12; QUERY-019), a `format` warning with `code: 'query'` under `smiles`
 * and `cxsmiles` saying why; `inchi` and `inchiKey` also for a drawing InChI cannot hold (spike inchi, N6), which `warnings`
 * says why. If the engine cannot write a drawing that has atoms in a format, reading that getter throws a `SketchError`
 * (`{ format, message }`), never an `Error`, and throws it again until the drawing changes; `cdxml`, `inchi` and `inchiKey`
 * throw one with a `code` while their module loads or after it failed to (`loading`, `unavailable`).
 */
export interface SketcherApi {
    /**
     * Canonical SMILES, as RDKit writes it, or CXSMILES when plain SMILES would lose something CXSMILES holds
     * (decision `smiles-cx`, API-055): enhanced stereo (ABS, AND and OR groups), attachment points, aliases, atom
     * labels such as generic atoms, or other atom properties read from a CXSMILES. Then it is the `cxsmiles`
     * getter's text, character for character; otherwise it is plain SMILES, as `plainSmiles` is. R labels of one
     * number (`[*:n]`), radicals, isotopes, charges, atom maps and stereo are held by plain SMILES. RDKit, and so
     * Datagrok, reads either. `""` for an empty drawing. Chosen once per drawing state.
     *
     * Query features (spike query, Q12; QUERY-019): generic atoms alone (A, Q, X, M and their H variants) are written as
     * CXSMILES labels; any other query feature (an atom list, a query property, a custom query, a query bond) has no
     * SMILES, and this is `""`, never a simplified molecule: `warnings` names it (`code: 'query'`), and `smarts` and the
     * molfiles hold it.
     */
    readonly smiles: string;
    /**
     * Always plain SMILES, canonical as RDKit writes it, with no CXSMILES extension: what `smiles` returns whenever
     * plain SMILES loses nothing. It loses something exactly when it differs from `smiles`. `""` for an empty drawing, and
     * for a query feature other than generics, as `smiles` (spike query, Q12); generics alone it leaves out (`*C`).
     */
    readonly plainSmiles: string;
    /** MOL V2000 (V3000 when the molecule needs it); `""` for an empty drawing. */
    readonly molfile: string;
    /** MOL V3000; `""` for an empty drawing. */
    readonly molV3000: string;
    /** Canonical CXSMILES; `""` for an empty drawing, and for a query feature other than generics (spike query, Q12), as `smiles`. */
    readonly cxsmiles: string;
    /**
     * SMARTS, as RDKit's MolToSmarts writes the same query (spike query, Q10; QUERY-017): each query feature as RDKit writes
     * it, drawn hydrogens merged into their heavy atoms where the drawing has a query feature (P10, QUERY-022); a drawing
     * without one gives RDKit's explicit-atom SMARTS (QUERY-018).
     */
    readonly smarts: string;
    /**
     * The drawing as SVG (API-022, IO-032): drawn by the canvas's own renderer, in its style and atom colours, at its
     * natural size (the drawing's 100%, whatever the view's zoom), cut to what it paints; without hover, selection,
     * hotspot or warning marks, hit targets or test ids. Each atom and bond group carries an id (`<prefix>atom-<i>`,
     * `<prefix>bond-<i>`), and its ids get a prefix new with each read, so that SVGs put on one page never share an id.
     * Synchronous and current at every moment, as the other getters are; `""` for an empty drawing. It starts at its
     * `<svg>`, with no XML prolog, ready to put in a page or to save as a file. Its labels are glyph outlines: real text
     * waits for the engine.
     */
    readonly svg: string;
    /**
     * The drawing as CDXML (spike cdxml, X4; API-008, IO-020): ChemDraw's ACS Document 1996 style, its colour and font tables,
     * y down at 14.4 pt per bond (IO-018), labels written as ChemDraw writes them, enhanced stereo, R labels as generic
     * nicknames, the wedges as drawn; ChemDraw opens it as the same structures. Synchronous and current at every moment, written
     * once per drawing state; `""` for an empty drawing. CDXML is written by a module of its own, which the page loads lazily
     * the first time a sketcher meets CDXML (a CDXML value, paste or drop, Copy As CDXML, Ctrl+D, a read of this getter, or
     * `loadFormat('cdxml')`): while it loads, reading this for a drawing with atoms throws a `SketchError` (`format:
     * 'cdxml'`) and starts the load, and `await loadFormat('cdxml')` makes it synchronous from then on. What CDXML cannot
     * hold of the drawing is in `warnings` under `cdxml` once the module is in.
     */
    readonly cdxml: string;
    /**
     * The drawing as standard InChI (spike inchi, N1; API-008, IO-024), byte-equal to RDKit.js's `get_inchi()` of the
     * `molfile` getter's text: RDKit's `MolToInchi` with no options (N5). Synchronous and current at every moment, written once
     * per drawing state; `""` for an empty drawing (IO-026), and for a drawing InChI cannot hold rather than a wrong
     * identifier (R labels, dummy atoms and pseudoatoms, query atoms and bonds, dative, hydrogen and zero-order bonds, N6), which
     * `warnings` says why under `inchi`, with what InChI holds only with a loss (enhanced stereo written as absolute, attachment
     * points, aliases, S-groups) and the InChI library's own warnings. InChI is written by a module of its own, which the page
     * loads lazily the first time a sketcher meets InChI (an InChI value, paste or drop, Copy As InChI or InChIKey, a read of
     * this getter or `inchiKey`, or `loadFormat('inchi')`): while it loads, reading this for a drawing with atoms throws a
     * `SketchError` (`format: 'inchi'`, `code: 'loading'`) and starts the load, and `await loadFormat('inchi')` makes it
     * synchronous from then on. A library error (more than 1,024 atoms, say) throws a `SketchError` too.
     */
    readonly inchi: string;
    /**
     * The drawing's InChIKey (spike inchi, N1; API-008): the key of the `inchi` getter's InChI, as RDKit.js's
     * `get_inchikey_for_inchi` writes it; `""` when that InChI is `""`. It loads, throws and warns as `inchi` does.
     */
    readonly inchiKey: string;
    /**
     * Loads the module a format needs, and resolves once that format is read and written synchronously (spikes cdxml, inchi):
     * `'cdxml'` loads the page's CDXML module and `'inchi'` its InChI module, each once per page, and waits for any write held
     * for it; every other format needs none and resolves at once. Rejects with a `SketchError` (`format: 'cdxml'` or
     * `'inchi'`, `code: 'unavailable'`) when the module cannot be loaded. A CDXML or InChI value written before then waits for
     * its module, the writes after it in turn, each applied with its change event once it is in; until then every getter,
     * `value` included, answers for the drawing shown. An initial CDXML or InChI value is read before the sketcher is ready.
     */
    loadFormat(format: SketcherFormat): Promise<void>;
    /**
     * The exact text last accepted, until the first user edit (API-052), and again when every edit since
     * has been undone (while the history still holds the molecule as written). After an edit it is the
     * `molfile` getter's text (API-054): MOL V2000, or V3000 whenever V2000 would lose something V3000 keeps,
     * which `warnings` then says under `molV2000` with the code `v3000`; `""` when the edit emptied the drawing.
     * Writing it replaces the molecule and starts a new history.
     */
    value: string;
    /**
     * Writes `text` as `value` does, read as `options.format` when it is given instead of detected from its shape (API-010:
     * SMARTS that is also SMILES). One change event with `{ source: 'api' }` when it is accepted; a text it cannot read
     * changes nothing, sets `lastError` and emits `error`, never throwing (API-013). Throws on a destroyed sketcher.
     */
    setValue(text: string, options?: SetValueOptions): void;
    /** The last refused write's error, until the next accepted write clears it to `null`. */
    readonly lastError: SketchError | null;
    /**
     * Every problem of the molecule as it is now, then what the getters' formats cannot hold of it, `[]` when there is
     * neither (spike charges-checks, K11, CHECK-012; decision `rgroups`). First the problems, in atom order: an atom
     * whose valence RDKit's sanitization rejects (`kind: 'valence'`), aromatic atoms with no Kekulé structure
     * (`'aromaticity'`), a label kept as an alias (`'label'`), a stereocentre nothing specifies in a structure drawn with
     * a wedge, a hashed wedge or a wavy bond (`'stereo'`, "Unspecified stereocentre"; spike stereo, S13), each with the
     * status line's words and every atom concerned. Then the losses (`'format'`), one per format and loss: a combined R
     * label (R1,R2), which SMILES, CXSMILES and SMARTS write with its lowest number, and attachment points, which SMARTS
     * does not write. Under `smiles`, what the `smiles` getter's text loses (decision `smiles-cx`): it writes CXSMILES
     * for attachment points, so they warn under `smarts` alone; `plainSmiles` has no entries, since it loses something
     * exactly when it differs from `smiles`. Nothing is refused for any of them. Synchronous and current at every moment,
     * as the getters are; it changes only with the molecule, so a `change` listener reads it as it reads the getters.
     * Reading it emits nothing.
     *
     * Spike formats (F2, F11, F13): every format's every loss comes from the engine's one loss table, each `format`
     * warning with its row's `code` and the `bonds` concerned. Under `molV2000`, `code: 'v3000'` says that the molfile
     * getter writes V3000 instead, and why (more than 999 atoms or bonds, stereo groups, a dative bond, …; IO-014).
     */
    readonly warnings: readonly SketchWarning[];
    /**
     * What the last accepted `value` write, paste or dropped file dropped on reading (spike formats, F10, F11; IO-030): an SD
     * record's data fields, the records after the first, an RG file's R-group definitions, and what the engine reports (an
     * unknown S-group type, 3D coordinates laid out afresh, …). Each a `format` warning: `format` the format read, `code`,
     * `message`, `atoms` and `bonds` empty. `[]` when the read kept everything. Replaced by the next accepted write, paste or
     * drop; kept through every other edit and through `destroy()`. No event of its own: the write's `change` is the moment
     * to read it. Reading it emits nothing.
     */
    readonly readWarnings: readonly SketchWarning[];
    /**
     * Whether the `cxsmiles` getter writes the drawing's coordinates and its wedges as drawn (spike formats, F6; IO-009):
     * `false` (the default), the canonical CXSMILES without them, comparable from drawing to drawing; `true` adds `(x,y,)`
     * and `wU`/`wD`. `smiles` and `plainSmiles` never carry them; Copy As CXSMILES always does, so that a paste lands as
     * drawn. The `cxsmiles-coordinates` attribute sets it too. Setting it is no edit: no change event, no undo step, and only
     * the `cxsmiles` getter's text changes.
     */
    cxsmilesCoordinates: boolean;
    /**
     * The formula and masses of the whole drawing (spike charges-checks, K13; INFO-010), as the readout shows them:
     * `{ formula, molecularWeight, exactMass, fragments }`. An empty drawing's formula is `""` and its masses 0; `null` for
     * a query (INFO-005) and when the engine cannot say. Synchronous and current at every moment, read once per molecule;
     * reading it emits nothing.
     */
    readonly analysis: SketchAnalysis | null;
    /**
     * The same for the atoms selected and the two atoms of each bond selected, each with its own hydrogens (INFO-002);
     * `null` with nothing selected. The selection is the user's, not the API's.
     */
    readonly selectionAnalysis: SketchAnalysis | null;
    /**
     * Whether the canvas marks the atoms with problems (CHECK-011): `'shown'` (the default) or `'hidden'`, which hides the
     * marks and nothing else. It reads back the value in force as `atomColours` does; setting it redraws the marks at once,
     * and is no edit: no change event, no undo step, no getter changes.
     */
    warningMarks: WarningMarks;
    /**
     * The formula and mass readout at the canvas's foot (UI-030): `'shown'` (the default), `'collapsed'` (its toggle alone)
     * or `'hidden'`. The toggle changes it too. It reads back the value in force as `atomColours` does; no edit.
     */
    analysisPanel: AnalysisPanel;
    /** Whether the drawing is empty (API-019), as it is after a clear, a write of `""` and the last atom erased. A read. */
    readonly isEmpty: boolean;
    /**
     * Whether the drawing has a query feature (spike query, Q14; API-020): a generic atom, an atom list or NOT list, a query
     * property, a custom query, a query bond or a bond topology; a plain `*` is none. Synchronous and current at every moment;
     * a read, which emits nothing. After `destroy()`, what it last said.
     */
    readonly hasQuery: boolean;
    /** Whether a selected atom or bond has a query feature (API-020); false with nothing selected. The selection is the user's. */
    readonly selectionHasQuery: boolean;
    /**
     * `'molecule'` (the default) or `'query'` (API-021), switched at any time, on the handle as on `<crux-sketch>` (spike
     * query, Q2, P2; spike query-edit, E5), where it is an attribute too: query mode offers the query tools (the query bonds,
     * the periodic table's lists and generics), the menus' query items, the dialogs' query sections and SMARTS pasted;
     * molecule mode hides them, and still reads, draws, keeps and writes a value's query features. Switching it changes no
     * molecule and is no edit (no change event, no undo step, the getters as they were); any other value is `'molecule'`.
     */
    mode: SketcherMode;
    /**
     * The query elements the host turns off (API-028; spike query-edit, E5): generic symbols (`A` … `MH`, `*`), `list`,
     * `not-list` and the query bond names (`any`, `single-or-double`, `single-or-aromatic`, `double-or-aromatic`). Each is
     * shown disabled where it is offered, saying why, and refused, saying why, where it is typed. It reads back the names it
     * took, in that order; anything else is left out. No edit: no change event, no undo step, the getters as they were.
     */
    disabledQueryElements: readonly string[];
    /** The sketcher's version (packages/sketch/package.json). */
    readonly version: string;
    /** The crux-core revision the engine was built from. */
    readonly engineRev: string;
    /** How many times the canvas has rendered (API-042); mirrored as `data-render-count` on the canvas. */
    readonly renderCount: number;
    /**
     * Where the canvas draws each atom and bond now (API-056; spike datagrok-chem, R1): for a host that points at them, as
     * Datagrok's hit areas do. A read, which renders and emits nothing; current at every moment. Empty before the sketcher
     * is ready, for an empty drawing and after `destroy()`.
     */
    readonly positions: DrawnPositions;
    /**
     * Whether the sketcher has something on its way (API-057; spike datagrok-chem, R1): a press or a drag held on the canvas,
     * or a drawing, a preview, a view change, a tooltip or a read of the clipboard waiting for a coming frame, a short wait
     * (a ring's preview waits for the pointer to rest) or a promise. False when the canvas shows all it will until the next
     * input or write: a host or a test that waits for the sketcher to settle waits on this, never on a timeout. A notice's
     * time on the status line is not counted. False after `destroy()`.
     */
    readonly isPending: boolean;
    /**
     * Where the sketcher shows its canvas, its toolbars, each control they show and the label editor now (API-058; spike
     * datagrok-platform, R3, R4): for a host that points at them, as Datagrok's hit areas do. A read, which renders and
     * emits nothing; current at every moment. Empty (no canvas, no toolbars, no controls) before the sketcher is ready, while
     * it is out of the page and after `destroy()`.
     */
    readonly layout: SketcherLayout;
    /**
     * The atoms and bonds selected (API-059; spike datagrok-platform, R4), by index, ascending; none with nothing selected,
     * before the sketcher is ready and after `destroy()`. The selection is the user's: the API reads it and never sets it. A
     * read, which emits nothing.
     */
    readonly selection: SketchSelection;
    /**
     * The tool chosen now (API-060; spike datagrok-platform, R4), however it was chosen (its button, a key, the canvas menu),
     * named as a toolbar configuration names it; `bond.single` at first. After `destroy()`, the last one. A read.
     */
    readonly tool: SketcherTool;
    /**
     * Undoes the last completed edit (API-043) and returns whether it did anything: with nothing to
     * undo it returns `false` and changes nothing. When it undoes, it emits one `change` with
     * `{ source: 'api' }`. Ctrl+Z and the Undo button use the same history, and emit `{ source: 'user' }`.
     * A value the host writes starts a new history, so this never goes back past it. The history holds
     * at most 500,000 atoms plus bonds across its steps; beyond that its oldest steps are dropped. Throws
     * on a destroyed sketcher, as a write does.
     */
    undo(): boolean;
    /** Redoes the last edit that was undone: as `undo()`, the other way. A new edit clears what could be redone. */
    redo(): boolean;
    /** Whether `undo()` would do anything: the history holds a completed edit since the host last wrote a value. Read-only. */
    readonly canUndo: boolean;
    /** Whether `redo()` would do anything: an edit was undone, and none has been made since. Read-only. */
    readonly canRedo: boolean;
    /**
     * How atoms are coloured (CANVAS-032): `'element'` (the default) or `'black'`. It reads back the value
     * in force, by HTML's rule for enumerated attributes: `'black'` only when set to that keyword (in any
     * case), and `'element'` for anything else, nothing set included. Setting it redraws at once, and is
     * no edit: no change event, no undo step, and no getter changes.
     */
    atomColours: AtomColours;
    /**
     * The keymap's profile (H32): `'chemdraw'` (the default) or `'ketcher'`, which adds Ketcher's keys where they do not
     * clash (so far Ctrl+Shift+M, Copy As MOL). It reads back the value in force, as `atomColours` does. Setting it changes
     * the keys at once, and the keys the tooltips and menus name; it is no edit: no change event, no getter changes.
     */
    keymap: KeymapProfile;
    /**
     * The CIP labels (spike stereo, S10; STEREO-007 to STEREO-009): `'hidden'` (the default, P6) or `'shown'`, which draws
     * (R), (S), (r), (s), (E) and (Z) beside the stereocentres and stereo double bonds, the engine's IUPAC 2013 labels, drawn
     * again with every edit, undo and redo. It reads back the value in force, as `atomColours` does: `'shown'` only when set
     * to that keyword. Setting it redraws at once, and is no edit: no change event, no undo step, and no getter changes.
     * The canvas menu's "Show stereochemistry" sets it too. Each sketcher has its own.
     */
    stereoLabels: StereoDisplay;
    /**
     * The stereo flag (spike stereo, S8; STEREO-018): `'shown'` (the default) or `'hidden'`. The flag is one per drawing, in
     * ChemDraw's words: Racemic when every specified stereocentre is in one AND group, Relative when in one OR group, else
     * Absolute when the chiral flag is set. It reads back the value in force: `'hidden'` only when set to that keyword.
     * Setting it redraws at once, and is no edit.
     */
    stereoFlag: StereoDisplay;
    /**
     * How a value's double bond with no E/Z of its own (a SMILES' `CC=CC`) is drawn (spike stereo, P2): `'crossed'` (the
     * default) or `'plain'`; its menu checks Unknown either way, and the getters write it unspecified. A double bond the user
     * draws takes its E/Z from its drawing, and a file's "either" bond and a bond the user crossed stay crossed, whatever it
     * says. It reads back the value in force: `'plain'` only when set to that keyword. Setting it redraws at once, and is no
     * edit: no change event, no undo step, and no getter changes. The `unspecified-double-bonds` attribute sets it too.
     */
    unspecifiedDoubleBonds: UnspecifiedDoubleBonds;
    /**
     * The locale the sketcher speaks (spike keyboard-a11y, K5; I18N-002, API-032): a BCP 47 tag, `'en'` by default, and
     * `'en-XA'` the pseudo-locale, English accented, about 40% longer and bracketed (I18N-003). Its messages are `strings`,
     * English where they leave one out, with plural forms by the locale (I18N-004); test ids, element symbols, SMILES, R
     * labels and key names never change (I18N-005, I18N-006). It reads back the locale in force as `Intl` writes it, and
     * `'en'` for anything `Intl` cannot read. Setting it is no edit: no change event, no getter changes.
     */
    locale: string;
    /**
     * Translations of the sketcher's messages (I18N-002), by their names in the catalog (`packages/sketch/src/core/strings.ts`):
     * text, a template with `{name}` placeholders, or plural forms (`{ one, other }`, `{count}` the number); element names
     * as `element.<Symbol>`. It reads back what the sketcher took (unknown names and malformed messages left out). Setting
     * it is no edit, as `locale`.
     */
    strings: SketcherStrings;
    /**
     * Terminal carbons labelled CH₃ (spike settings-touch, K6, P10): `'hidden'` (the default) or `'shown'`. It reads back the
     * value in force, as `atomColours` does: `'shown'` only when set to that keyword. Setting it redraws at once, and is no
     * edit: no change event, no undo step, and no getter changes.
     */
    methylLabels: MethylLabels;
    /**
     * How large the toolbars' controls are (spike bar-fit, B1; UI-005, HOST-029): `'auto'` (the default) compact in a box under
     * 720 px wide or 600 px tall and regular otherwise, following the box as it changes; `'regular'` or `'compact'` in any
     * box. It reads back the value in force, as `atomColours` does: `'regular'` or `'compact'` only when set to that keyword.
     * Setting it changes the controls' sizes at once, and nothing else: no change event, no undo step, no getter changes. The
     * `density` attribute sets it too. Each sketcher has its own.
     */
    density: Density;
    /**
     * Whether the user's settings are kept in the page's localStorage and shared with the page's other sketchers (spike
     * reports-oct5, R5): `true`, the default, or `false`, read back as the value in force (`false` only when set so). Set
     * before the sketcher is ready, it decides whether the kept settings are taken at creation; turned off later, the
     * sketcher stops keeping and hearing changes; turned on, it starts again. The `persist-settings` attribute sets it too.
     */
    persistSettings: boolean;
    /**
     * The host's keys (spike settings-touch, K7, KEY-004, P11): each key, as written ("Ctrl+L", "Shift+T", "t", "F2";
     * Ctrl means Ctrl or Cmd), bound to a built-in action (`layout`, `undo`, `keys.open`, `tool.erase`,
     * `tool.ring.cyclohexane`, `element.N`, `bond.double`, `atom.properties` …) or to `null`, which turns the key off so that
     * it reaches the page. A key named is the host's alone, whatever the profile gave it. It reads back what the sketcher
     * took (keys it cannot read and actions it does not know left out). The tooltips, menus and shortcuts sheet follow it.
     */
    keys: HostKeys;
    /**
     * The settings (spike settings-touch, K6, P9; UI-021, UI-022): `atomColours`, `warningMarks`, `analysisPanel`,
     * `stereoLabels`, `stereoFlag`, `unspecifiedDoubleBonds`, `methylLabels` and `keymap`, read and written whole: what a
     * written object leaves out takes its default. Each is its own option as well. A user's change emits `settings` with
     * them, for a host to store; and, unless `persistSettings` is off, the sketcher keeps the user's changes in the page's
     * localStorage and shares them with the page's other sketchers (spike reports-oct5, R5), where a host's own value wins.
     * No edit: no change event, no undo step, no getter changes.
     */
    settings: SketcherSettings;
    /** Adds a listener; returns the function that removes it. */
    on<T extends SketcherEventType>(type: T, listener: SketcherListener<T>): () => void;
    /** Removes a listener `on()` added (API-018); a listener it never added is ignored. */
    off<T extends SketcherEventType>(type: T, listener: SketcherListener<T>): void;
    /** Releases the sketcher. The getters keep answering; writes throw. Calling it twice is harmless. */
    destroy(): void;
}
/** A host's keys (KEY-004, P11): each key, as written ("Ctrl+L", "Shift+T", "t", "F2"), bound to a built-in action, or to nothing (`null`). */
export type HostKeys = Readonly<Record<string, string | null>>;
