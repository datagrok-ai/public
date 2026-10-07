import { type SketcherSettings } from '../core/settings.js';
import { type SketcherStrings } from '../core/strings.js';
import type { HostKeys, DrawnPositions, RenderDetail, SketcherLayout, SketcherTool, SketchSelection, AnalysisPanel, AtomColours, ChangeDetail, KeymapProfile, MethylLabels, Density, SetValueOptions, SketcherFormat, SketchAnalysis, SketchError, SketchWarning, WarningMarks, SketcherApi, SketcherEventType, SketcherListener, SketcherMode, StereoDisplay, ToolbarConfig, UnspecifiedDoubleBonds } from '../core/types.js';
export declare const TAG = "crux-sketch";
/**
 * The events of `<crux-sketch>`: an element's, with `change` and `error` as CustomEvents whose
 * `detail` is what the `on()` listeners receive (API-053). `error` is not lib.dom's `ErrorEvent`.
 */
export interface CruxSketchEventMap extends Omit<HTMLElementEventMap, 'change' | 'error'> {
    change: CustomEvent<ChangeDetail>;
    error: CustomEvent<SketchError>;
    /** A user's Apply in the Settings dialog changed the settings (spike settings-touch, K6): the settings in force. */
    settings: CustomEvent<SketcherSettings>;
    /** The canvas shows a new drawing (API-061; spike datagrok-platform, R5): its render count. */
    render: CustomEvent<RenderDetail>;
}
/** `<crux-sketch>`'s listeners, typed by its events (`CruxSketchEventMap`): `change`, `error`, `settings`, `render`. */
export interface CruxSketchElement {
    addEventListener<K extends keyof CruxSketchEventMap>(type: K, listener: (this: CruxSketchElement, event: CruxSketchEventMap[K]) => unknown, options?: boolean | AddEventListenerOptions): void;
    addEventListener(type: string, listener: EventListenerOrEventListenerObject, options?: boolean | AddEventListenerOptions): void;
    removeEventListener<K extends keyof CruxSketchEventMap>(type: K, listener: (this: CruxSketchElement, event: CruxSketchEventMap[K]) => unknown, options?: boolean | EventListenerOptions): void;
    removeEventListener(type: string, listener: EventListenerOrEventListenerObject, options?: boolean | EventListenerOptions): void;
}
/**
 * `<crux-sketch>` (H29; API-002), defined when the module is imported: the sketcher as a custom element, with every member
 * of the handle `createSketcher` gives (`SketcherApi`), and `ready`, which resolves once it is ready, and `toolbar`, its
 * toolbars' configuration. Its attributes are `value`, `mode`, `atom-colours`, `keymap`, `warning-marks`, `analysis-panel`,
 * `stereo-labels`, `stereo-flag`, `unspecified-double-bonds`, `locale`, `methyl-labels`, `cxsmiles-coordinates`, `density`
 * and `persist-settings`, each setting the property of its name; changes made before it is ready are applied in order
 * once it is (API-003). Its events are DOM events too (`change`, `error`, `settings`, `render`; API-053, API-061). Its
 * typings declare it an `HTMLElement` with this API: Lit, which it is built with and which is bundled, is no part of its
 * type, so that a host type-checks it with nothing installed (spike api-v1, R1).
 */
export declare class CruxSketchElement extends HTMLElement implements SketcherApi {
    #private;
    /**
     * Resolves once the sketcher is ready: the engine loaded, the value applied and drawn. Rejects if
     * the engine cannot be loaded (which is also reported as an `error` event).
     */
    readonly ready: Promise<void>;
    /**
     * How atoms are coloured: `'element'` (the default) or `'black'` (H36). It reads back the value in
     * force, by HTML's rule for enumerated attributes (fix round 3): `'black'` only when set to that
     * keyword (in any case), `'element'` for anything else. Setting it redraws at once: no change event,
     * no undo step, and the getters as they were. The `atom-colours` attribute sets it too, and keeps
     * what was written. The element palette's symbols, the periodic table's (open or not) and the label
     * preview follow it as the canvas does (ATOM-033), each sketcher by its own.
     */
    get atomColours(): AtomColours;
    set atomColours(value: AtomColours);
    /**
     * The CIP labels (S10, STEREO-009): `'hidden'` (the default, P6) or `'shown'`, read back as the value in force by HTML's
     * rule for enumerated attributes, as `atomColours` is. Setting it redraws at once: no change event, no undo step, and
     * the getters as they were. The `stereo-labels` attribute sets it too, and so does the canvas menu's Show stereochemistry.
     */
    get stereoLabels(): StereoDisplay;
    set stereoLabels(value: StereoDisplay);
    /**
     * The stereo flag (S8, STEREO-018): `'shown'` (the default) or `'hidden'`, read back as the value in force. Setting it
     * redraws at once, and is no edit. The `stereo-flag` attribute sets it too.
     */
    get stereoFlag(): StereoDisplay;
    set stereoFlag(value: StereoDisplay);
    /**
     * How a value's double bond with no E/Z of its own is drawn (spike stereo, P2): `'crossed'` (the default) or `'plain'`,
     * read back as the value in force. Setting it redraws at once, and is no edit. The `unspecified-double-bonds` attribute
     * sets it too.
     */
    get unspecifiedDoubleBonds(): UnspecifiedDoubleBonds;
    set unspecifiedDoubleBonds(value: UnspecifiedDoubleBonds);
    /**
     * The toolbar configuration: `{ top, left, right, bottom }`, each a list of groups `{ label, tools }`,
     * the right one the element palette (UI-002, spike atoms) and the bottom one the ring bar (spike
     * rings-chains). `{}` shows no toolbars; `null` restores the default. Set it before the sketcher is
     * ready.
     *
     * It is presentation, not permission: it decides which buttons show, and the keys that choose and use
     * tools do not read it. A host that leaves out the ring bar or the chain button still gets rings and
     * chains from j, Shift+J, Shift+X, the number keys and the ring keys (reading 31).
     */
    get toolbar(): ToolbarConfig;
    set toolbar(config: ToolbarConfig | null | undefined);
    /**
     * The query elements the host turns off (API-028; spike query-edit, E5): generic symbols (`A` … `MH`, `*`), `list`,
     * `not-list` and the query bond names (`any`, `single-or-double`, `single-or-aromatic`, `double-or-aromatic`). Each is shown
     * disabled where it is offered (the toolbar's tools, the periodic table, the menus, Bond Properties), its tooltip saying
     * why, and refused, saying why, where it is typed (the label editor, x). It reads back the names it took, in that order;
     * anything else is left out. No edit: no change event, the getters as they were.
     */
    get disabledQueryElements(): readonly string[];
    set disabledQueryElements(value: readonly string[]);
    /** Plain SMILES, or CXSMILES where plain SMILES would lose something (spike smiles-cx, X1): the core's. */
    get smiles(): string;
    /** Always plain SMILES (spike smiles-cx, X3): the core's. */
    get plainSmiles(): string;
    get molfile(): string;
    get molV3000(): string;
    get cxsmiles(): string;
    get smarts(): string;
    /** The drawing as SVG (spike copy-as, P5; API-022, IO-032): the core's, its ids new with each read. */
    get svg(): string;
    /** The drawing as CDXML (spike cdxml, X4; API-008): the core's. It throws while the page's CDXML module loads: `loadFormat('cdxml')`. */
    get cdxml(): string;
    /** The drawing as standard InChI (spike inchi, N1; API-008): the core's. It throws while the page's InChI module loads: `loadFormat('inchi')`. */
    get inchi(): string;
    /** The drawing's InChIKey (spike inchi, N1; API-008): the core's, the key of the `inchi` getter's InChI; it throws as that getter does. */
    get inchiKey(): string;
    /** Loads the module a format needs (spikes cdxml, inchi: CDXML's, InChI's): the core's; resolves once the format is read and written synchronously. */
    loadFormat(format: SketcherFormat): Promise<void>;
    /**
     * The keymap's profile (H32): `'chemdraw'`, the default, ChemDraw's keys, or `'ketcher'`, which adds Ketcher's (so far
     * Ctrl+Shift+M, Copy As MOL; spike 14 the rest). It reads back the value in force, by HTML's rule for enumerated
     * attributes, as `atomColours` does; the `keymap` attribute sets it too. Setting it changes the keys at once, and the
     * keys the tooltips and the menus name: no change event, no undo step, and the getters as they were.
     */
    get keymap(): KeymapProfile;
    set keymap(value: KeymapProfile);
    /**
     * The host's keys over the profile (spike settings-touch, K7, KEY-004, P11): `{ "Ctrl+L": "layout", "t": null }`, each
     * key it names its own, bound to a built-in action or turned off. It reads back what the sketcher took. Setting it
     * changes the keys at once, and the keys the tooltips, the menus and the shortcuts sheet name; it is no edit.
     */
    get keys(): HostKeys;
    set keys(value: HostKeys);
    /**
     * Terminal carbons labelled CH₃ (K6, P10): `'hidden'` (the default) or `'shown'`, read back as the value in force by
     * HTML's rule for enumerated attributes. Setting it redraws at once, and is no edit. The `methyl-labels` attribute sets it too.
     */
    get methylLabels(): MethylLabels;
    set methylLabels(value: MethylLabels);
    /**
     * How large the toolbars' controls are (spike bar-fit, B1): `'auto'` (the default), `'regular'` or `'compact'`, read back as
     * the value in force by HTML's rule for enumerated attributes. Setting it redraws the bars at once, and is no edit: no
     * change event, no undo step, no getter changes. The `density` attribute sets it too.
     */
    get density(): Density;
    set density(value: Density);
    /**
     * The settings (K6, P9; UI-021, UI-022): every view option at once, read and written whole, what a written object leaves
     * out taking its default. Writing them is no edit and says nothing: the `settings` event is the user's Apply's alone.
     */
    get settings(): SketcherSettings;
    set settings(value: SketcherSettings);
    /**
     * Whether the user's settings are kept and shared (spike reports-oct5, R5): `true`, the default, or `false`, which neither
     * reads nor writes the page's storage and neither hears nor tells the page's other sketchers. The `persist-settings`
     * attribute sets it too (`"false"` turns it off). Set before the sketcher is ready, it decides whether the kept settings
     * are taken; turned off later, the sketcher stops keeping and hearing; turned on later, it starts again.
     */
    get persistSettings(): boolean;
    set persistSettings(value: boolean);
    /**
     * The locale the sketcher speaks (K5; I18N-002, API-032): a BCP 47 tag, `'en'` by default; `'en-XA'` is the
     * pseudo-locale, English accented, longer and bracketed (I18N-003). Its messages are the host's `strings`, English where
     * they leave one out; plurals follow the locale (I18N-004). It reads back the locale in force as `Intl` writes it, and
     * `'en'` for anything `Intl` cannot read. The `locale` attribute sets it too. Setting it shows the chrome in it at once:
     * no change event, no undo step, and the getters as they were. Each sketcher has its own.
     */
    get locale(): string;
    set locale(value: string);
    /**
     * The host's translations (K5, I18N-002): messages of the catalog by name, text, `{name}` templates or plural forms, and
     * element names as `element.<Symbol>`; what they leave out is English. It reads back what the sketcher took of them
     * (unknown names and malformed messages left out), frozen. Setting it shows the chrome in them at once, as `locale` does.
     */
    get strings(): SketcherStrings;
    set strings(value: SketcherStrings);
    get value(): string;
    set value(text: string);
    setValue(text: string, options?: SetValueOptions): void;
    get lastError(): SketchError | null;
    get warnings(): readonly SketchWarning[];
    /** What the last accepted write, paste or drop dropped on reading (spike formats, F10, F11; IO-030): the core's. */
    get readWarnings(): readonly SketchWarning[];
    /**
     * Whether the `cxsmiles` getter writes the coordinates and the wedges as drawn (spike formats, F6; IO-009): `false`, the
     * default, or `true`; the `cxsmiles-coordinates` attribute sets it too (present: on, unless it says `false`). Setting it
     * is no edit: no change event, no undo step, and only the `cxsmiles` getter's text changes.
     */
    get cxsmilesCoordinates(): boolean;
    set cxsmilesCoordinates(on: boolean);
    /** The formula and masses of the drawing (spike charges-checks, K13; INFO-010): the core's. */
    get analysis(): SketchAnalysis | null;
    /** The formula and masses of the selected atoms and of the selected bonds' atoms (INFO-002, INFO-010); null with nothing selected. */
    get selectionAnalysis(): SketchAnalysis | null;
    /** Whether the drawing has a query feature (spike query, Q14, API-020): the core's; after `destroy()`, what it last said. */
    get hasQuery(): boolean;
    /**
     * Whether a selected atom or bond has a query feature (spike query, Q14, API-020): the engine's `hasQueryIn` of the atoms
     * and bonds selected; false with nothing selected, and after `destroy()`.
     */
    get selectionHasQuery(): boolean;
    /**
     * Whether the canvas marks the atoms with problems (spike charges-checks, K10, CHECK-011): `'shown'`, the default, or
     * `'hidden'`, which hides the marks only (`warnings`, `data-problem` and the status line stay). It reads back the value
     * in force, by HTML's rule for enumerated attributes, as `atomColours` does; the `warning-marks` attribute sets it too.
     * Setting it marks the canvas again at once: no change event, no undo step, and the getters as they were.
     */
    get warningMarks(): WarningMarks;
    set warningMarks(value: WarningMarks);
    /**
     * The formula and mass readout at the canvas's foot (spike charges-checks, K12; UI-030): `'shown'`, the default,
     * `'collapsed'`, its toggle alone, or `'hidden'`. The toggle sets it too, and so does the `analysis-panel` attribute; it
     * reads back the value in force, by HTML's rule for enumerated attributes. No edit: no change event, no undo step.
     */
    get analysisPanel(): AnalysisPanel;
    set analysisPanel(value: AnalysisPanel);
    get isEmpty(): boolean;
    get mode(): SketcherMode;
    /**
     * Query mode or molecule mode (spike query, Q2, P2): switched at any time; no edit, no change event, the getters as they
     * were. Query mode shows the query tools, the menus' query items and the dialogs' query sections, and reads SMARTS pasted;
     * molecule mode hides them (spike query-edit), a query tool in hand giving way to the single bond tool.
     */
    set mode(mode: SketcherMode);
    get version(): string;
    get engineRev(): string;
    get renderCount(): number;
    /**
     * Where the canvas draws each atom and bond now (API-056; spike datagrok-chem, R1): the editor's drawing, in CSS px from
     * this element's top-left corner (core/positions.ts). Empty before the sketcher is ready, for an empty drawing and after
     * `destroy()`. A read: it renders and emits nothing.
     */
    get positions(): DrawnPositions;
    /**
     * Whether something is on its way (API-057; spike datagrok-chem, R1): a press or a drag held on the canvas (the editor's
     * or the view's), or a frame, a task, a timer or a wait the sketcher's lifetime holds (a preview, a view change, a
     * tooltip's rest, a clipboard read); a notice's time aside. False after `destroy()`.
     */
    get isPending(): boolean;
    /**
     * Where the sketcher shows its parts now (API-058; spike datagrok-platform, R3, R4): the canvas, each toolbar, each
     * control a toolbar shows (by its name in a configuration, its test id without `toolbar.`) and the label editor while it
     * is open, in CSS px from this element's top-left corner (core/parts.ts). Empty before the sketcher is ready, out of the
     * page and after `destroy()`. A read: it renders and emits nothing.
     */
    get layout(): SketcherLayout;
    /**
     * The atoms and bonds the user has selected (API-059; spike datagrok-platform, R4), by index, ascending: the editor's
     * selection in the drawing shown. None before the sketcher is ready and after `destroy()`. A read.
     */
    get selection(): SketchSelection;
    /** The tool chosen now (API-060; spike datagrok-platform, R4), as a configuration names it; after `destroy()`, the last one. */
    get tool(): SketcherTool;
    undo(): boolean;
    redo(): boolean;
    get canUndo(): boolean;
    get canRedo(): boolean;
    on<T extends SketcherEventType>(type: T, listener: SketcherListener<T>): () => void;
    off<T extends SketcherEventType>(type: T, listener: SketcherListener<T>): void;
    /**
     * Releases the sketcher; the element stays where it is, empty and inert. Nothing it started outlives it (spike
     * host-isolation, I1): its lifetime ends first, cancelling every frame, task and timer still pending, disconnecting
     * its resize observer and letting go of every promise it awaits; then each part drops what it holds; and last the
     * page lays out (`releaseLayout`), so that its layout lets go of the chrome taken out of the shadow root at once.
     */
    destroy(): void;
}
declare global {
    interface HTMLElementTagNameMap {
        'crux-sketch': CruxSketchElement;
    }
}
