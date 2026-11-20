import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type * as CruxSketch from '../../vendor/crux-sketch/crux-sketch';
import type {SketcherFormat, SketcherHandle} from '../../vendor/crux-sketch/crux-sketch';
import '../../css/crux-sketcher.css';

type Notation = 'smiles' | 'molblock' | 'molblockV3000' | 'smarts';

const CRUX_FORMAT: {[n in Notation]: SketcherFormat} = {
  smiles: 'smiles',
  molblock: 'molV2000',
  molblockV3000: 'molV3000',
  smarts: 'smarts',
};

const DG_NOTATION: {[n in Notation]: DG.chem.Notation} = {
  smiles: DG.chem.Notation.Smiles,
  molblock: DG.chem.Notation.MolBlock,
  molblockV3000: DG.chem.Notation.V3KMolBlock,
  smarts: DG.chem.Notation.Smarts,
};

/** Crux Sketch (`vendor/crux-sketch/`), bundled with Chem like any dependency: a chunk of its own, loaded the first
 * time a Crux sketcher opens, its WebAssembly emitted beside Chem's chunks. */
export function loadCruxSketch(): Promise<typeof CruxSketch> {
  return import('../../vendor/crux-sketch/crux-sketch.js');
}

/** What the Crux sketchers share on the page, kept on the page, not in this module: Chem's package tests run in a
 * bundle of their own (package-test.js) beside package.js, whose functions make the sketchers the platform opens, and a
 * test's setter reaches those sketchers only so. With one bundle (production) it is this module's state. `loader`
 * (tests): what the sketchers load Crux Sketch with instead (null: {@link loadCruxSketch}), to hand a sketcher a handle
 * as an engine that has stopped gives one; it is given the sketchers' own loader, as the test bundle's copy of Crux
 * Sketch cannot make a sketcher once package.js's has defined its element. */
interface PageState {
  loader: CruxSketchLoaderForTests | null;
  /** The hosts whose Crux opens in query mode ({@link openCruxInQueryMode}); made on first use. */
  queryHosts?: WeakSet<DG.chem.Sketcher>;
}
const PAGE_STATE = Symbol.for('datagrok.chem.crux-sketcher.page');
const pageState: PageState = (globalThis as unknown as {[key: symbol]: PageState | undefined})[PAGE_STATE] ??=
  {loader: null};

/** Makes the Crux sketcher `host` brings up open in query mode, its query tools offered from the start, as the
 * substructure filter's does (Apps | Chem | Sketcher: the platform's sketcher host has no say in a sketcher's mode).
 * Call it before the host's sketcher starts; the gear's Query mode still switches it. */
export function openCruxInQueryMode(host: DG.chem.Sketcher): void {
  (pageState.queryHosts ??= new WeakSet()).add(host);
}

/** Tests: what the Crux sketchers load Crux Sketch with, given Chem's own loader (the one to wrap). */
export type CruxSketchLoaderForTests = (load: () => Promise<typeof CruxSketch>) => Promise<typeof CruxSketch>;

/** Tests: sets the loader every Crux sketcher created from now on uses (null: Chem's own). */
export function setCruxSketchLoaderForTests(loader: CruxSketchLoaderForTests | null): void {
  pageState.loader = loader;
}

/** Crux Sketch as the reaction sketchers load it (the Reactions app's, the reaction cell editor's): Chem's own
 * {@link loadCruxSketch}, or the tests' loader ({@link setCruxSketchLoaderForTests}), as the molecule sketcher does. */
export function loadCruxSketchForReactions(): Promise<typeof CruxSketch> {
  return pageState.loader?.(loadCruxSketch) ?? loadCruxSketch();
}

/** The name Crux's hit areas and readings go under in the sketcher's status (`getWidgetStatus()`). */
const STATUS = 'crux-sketch';

/** The sketcher's u2 name (`data-u2-name` on its root, HOST-043): what the platform's automation finds it by. */
const U2_NAME = 'cruxSketch';

/** What the sketcher's area says when the sketcher could not start (a user sees it: no engine's name and no error's
 * text; the error is in the console). The words of the reaction editors' (reactions-app.ts EDITOR_FAILED).
 * */
export const SKETCHER_NOT_STARTED = 'The drawing area could not start.';

/** Crux's toolbars by their place, as the u2 parts the status names them (HOST-044): the top bar holds the actions (undo,
 * clear, zoom, layout, the clipboard), the left one the drawing tools, the right one the elements, the bottom one the
 * templates (the rings and the Structure Library). */
const PARTS: {[place: string]: string} = {top: 'actions', left: 'tools', right: 'elements', bottom: 'templates'};

type Rect = {x: number, y: number, width: number, height: number};

function isEmptyValue(value: string | null | undefined): boolean {
  return !value?.trim() || grok.chem.Sketcher.isEmptyMolfile(value);
}

/** The notations Crux's getters write, by Chem's names. */
export type CruxNotation = Notation;

/** Why one of Crux's getters threw (the engine-poison foot-gun): `engine`, the engine stopped (it
 * trapped, or never loaded: crux-sketch's SketchError of format `engine`, or any throw that is no SketchError of a
 * notation, a raw RuntimeError among them); else the drawing cannot be written in `notation`. `message` is the
 * throw's. */
export interface CruxReadFailure {
  notation: Notation;
  engine: boolean;
  message: string;
}

/** A getter's throw, as a {@link CruxReadFailure}. */
export function readFailureOf(e: unknown, notation: Notation): CruxReadFailure {
  const error = e !== null && typeof e === 'object' ? e as {format?: unknown, message?: unknown} : null;
  const format = typeof error?.format === 'string' ? error.format : null;
  const engine = e instanceof Error || format === null || format === 'engine';
  const message = typeof error?.message === 'string' ? error.message : String(e);
  return {notation, engine, message};
}

/** The keys Crux Sketch acts on in its drawing that Datagrok's shortcuts act on too, at the window (core command.dart:
 * one keydown listener in the bubble phase, which reads no defaultPrevented): Select all (Ctrl+A) on the table behind a
 * dialog, undo and redo (Ctrl+Z, Ctrl+Y, Ctrl+Shift+Z), Delete and Backspace, the arrows. Never Escape, Enter or Tab:
 * a dialog's Cancel and OK, and focus moving on, are the page's. */
export function sketcherOwnsKey(e: KeyboardEvent): boolean {
  if (e.altKey)
    return false;
  if (e.ctrlKey || e.metaKey)
    return ['a', 'z', 'y'].includes(e.key.toLowerCase());
  return ['Delete', 'Backspace', 'ArrowUp', 'ArrowDown', 'ArrowLeft', 'ArrowRight'].includes(e.key);
}

/** Keeps the keys a Crux sketcher in `host` took (cancelled) or owns ({@link sketcherOwnsKey}) from the page's
 * shortcuts (Ctrl+A over the reaction cell editor's canvas also selected every row of the grid behind the dialog): a
 * keydown from inside `<crux-sketch>` stops at `host`, its keyup with it. Crux Sketch keeps the keys its handlers
 * cancel in its shadow root (its containment); this keeps any that get past it. `anywhere`: a dialog that holds a
 * sketcher keeps its keys pressed anywhere in it but in a text field (the dialog itself focused, as it opens), so that
 * the table behind never acts on them. Returns the listeners' removal. */
export function keepSketcherKeys(host: HTMLElement, anywhere: boolean = false): () => void {
  const kept = new Set<string>();
  const inSketcher = (e: Event) => e.composedPath().some((t) => (t as Element).localName === 'crux-sketch');
  const inText = (e: Event) => {
    const t = e.composedPath()[0] as HTMLElement | undefined;
    return !!t && (t.localName === 'input' || t.localName === 'textarea' || t.localName === 'select' ||
      t.isContentEditable === true);
  };
  const down = (e: KeyboardEvent) => {
    if (!(inSketcher(e) || (anywhere && !inText(e))) || !(e.defaultPrevented || sketcherOwnsKey(e)))
      return;
    e.stopPropagation();
    kept.add(e.code || e.key);
  };
  const up = (e: KeyboardEvent) => {
    if (kept.delete(e.code || e.key))
      e.stopPropagation();
  };
  host.addEventListener('keydown', down);
  host.addEventListener('keyup', up);
  return () => {
    host.removeEventListener('keydown', down);
    host.removeEventListener('keyup', up);
  };
}

/** Gives the keyboard to the drawing of the Crux sketcher in `host`, as a click on it does (its canvas, focusable):
 * a dialog that opens on a sketcher otherwise keeps the focus itself, and the keys go to the page (Ctrl+A among them).
 * Whether a canvas took it. */
export function focusDrawing(host: HTMLElement): boolean {
  const canvas = host.querySelector('crux-sketch')?.shadowRoot?.querySelector('svg[tabindex="0"]') as
    SVGElement | null | undefined;
  canvas?.focus({preventScroll: true});
  return !!canvas && canvas.getRootNode() instanceof ShadowRoot &&
    (canvas.getRootNode() as ShadowRoot).activeElement === canvas;
}

/** {@link focusDrawing} as a dialog on a sketcher opens, and again over its first second: the sketcher replaces its
 * canvas as it draws what it read (a ChemDraw page, a route arranged), and the focus fell back to the dialog; the page
 * behind (the grid the dialog opened from) never keeps it. Only a focus a person gave something else in the dialog (its
 * OK, a field) is left alone. */
export function focusDrawingSoon(host: HTMLElement, dialog: HTMLElement): void {
  const wait = [0, 100, 250, 500, 1000];
  const go = (i: number) => {
    if (!dialog.isConnected)
      return;
    const a = document.activeElement;
    if (a !== null && a !== dialog && dialog.contains(a) && !host.contains(a))
      return;
    focusDrawing(host);
    if (i + 1 < wait.length)
      setTimeout(() => go(i + 1), wait[i + 1] - wait[i]);
  };
  go(0);
}

/** Crux Sketch as a Datagrok molecule sketcher. Crux's getters are synchronous and stay readable after
 * `destroy()`, and its `change` event says whether a change is the user's edit or a value written to it, so
 * the host's rules map onto it directly: a setter is one write and one `onChanged`, a user edit is one
 * `onChanged` and ends `explicitMol`. Its root is named `cruxSketch` (`data-u2-name`), and its status
 * (`getWidgetStatus()`) shows the platform's tests what Crux draws and holds: a hit area for every atom and bond, its
 * parts and its controls, readings, and `isRenderPending` and `onRendered` to settle on. */
export class CruxSketcher extends grok.chem.SketcherBase {
  private sketch: SketcherHandle | null = null;
  // A value set while no live Crux could show it: before `init` resolved, or after `detach`.
  private unshown: {notation: Notation, value: string} | null = null;
  // Crux is being created: from `init` until it is ready or failed to start.
  private loading = false;
  // The host's change events since this sketcher was made (`onChanged`): the status's `changes` reading.
  private changes = 0;
  // The getters that threw, by notation, until each answers again (the engine-poison foot-gun): `read` answers
  // '' then, which no host may write over a value; the cell editor blocks OK and writes nothing (utils/footguns).
  private failures = new Map<Notation, CruxReadFailure>();
  // The sketcher started on an engine that had failed (crux-sketch's lastError of format `engine` as it started: an
  // engine that could make no fresh instance, the last resort of API-069, or a start the engine refused): it holds
  // nothing, and its getters may answer '' without a throw, so OK writes nothing for as long as it is open.
  private engineStop: CruxReadFailure | null = null;

  constructor() {
    super();
    this.root.classList.add('crux-sketcher');
    this.name = U2_NAME;
    DG.UndoService.ownScope(this.root);
    keepSketcherKeys(this.root); // Ctrl+A, Delete, the arrows: the drawing's, never the table's behind
    this.onChanged.subscribe(() => this.changes++);
    this.addStatusProvider(STATUS, () => this.status());
  }

  get type(): string {
    return 'Crux sketcher';
  }

  /** Whether Crux has something on its way (its `isPending`: a press held, a preview or a render waiting for a frame),
   * or is still loading: what the platform's tests settle on before they act or read. */
  get isRenderPending(): boolean {
    return this.loading || (this.sketch?.isPending ?? false);
  }

  /** Whether the drawing is a query no SMILES writes (an atom list, a query bond, a ring query): the `smiles` the host
   * reads of it is then empty, with the platform's warning (`read`), and a SMILES cell's OK would empty the cell. The
   * cell editor's foot-gun check asks this first, so that it never raises that warning itself (utils/footguns). */
  get smilesIsQuery(): boolean {
    if (this.explicitMol?.notation === 'smiles' || this.unshown !== null || this.sketch === null)
      return false;
    try {
      return holdsQuery(this.sketch);
    } catch {
      return false;
    }
  }

  /** What the last value write, paste or drop left out on reading, by crux-sketch's codes (`arrow`, `text`: a ChemDraw
   * reaction pasted in becomes one molecule; `sd-data`, `coordinates-3d`, ...): the cell editor's foot-gun check names
   * the arrow (utils/footguns, MOL-17). */
  get readLosses(): string[] {
    try {
      return (this.sketch?.readWarnings ?? []).map((w) => w.code ?? '').filter((c) => c);
    } catch {
      return [];
    }
  }

  async init(host: DG.chem.Sketcher): Promise<void> {
    this.host = host;
    let sketch: SketcherHandle;
    this.loading = true;
    try {
      // Where Datagrok asks for a query (its substructure filter) or a host asked for it (openCruxInQueryMode: the
      // Sketcher app), Crux opens in query mode, offering the query tools; every other host opens it in molecule mode.
      // The user can switch it in Crux's gear menu either way.
      const query = host.isSubstructureFilter || (pageState.queryHosts?.has(host) ?? false);
      sketch = await (await (pageState.loader?.(loadCruxSketch) ?? loadCruxSketch())).createSketcher(this.root,
        query ? {mode: 'query'} : undefined);
    } catch (e) {
      console.error(e);
      this.root.append(ui.divText(SKETCHER_NOT_STARTED, 'crux-sketcher-error'));
      return;
    } finally {
      this.loading = false;
    }
    // A sketcher started on an engine that had failed says so while it starts, before createSketcher resolves: no
    // listener made now hears it, but lastError holds it; such a sketcher holds nothing and its getters may answer ''
    // without a throw, so it is noted here (the never-write-empty guard; no block and no balloon: the engine recovers)
    let startError: {format?: string, message?: string} | null = null;
    try {
      startError = sketch.lastError;
    } catch (e) {
      startError = {message: String((e as {message?: string})?.message ?? e)};
    }
    if (startError !== null && (startError.format === 'engine' || startError.format === undefined))
      this.engineStop = readFailureOf(startError, 'smiles');
    if (this.isDetached) {
      // The host's dialog closed while Crux loaded (OK pressed at once): Crux lets go of the page, and answers as a
      // sketcher detached after it was ready does. Its getters read what it holds, and the value the host stored
      // meanwhile, which the host sets on it now, is kept and said once (`write`), so a host listening for changes
      // (the substructure filter, "Filter as you draw" on) still gets the SMILES typed into its field.
      sketch.destroy();
      this.sketch = sketch;
      return;
    }
    // Each drawing Crux shows (its `render`): the platform's `onRendered`, which tests settle on with `isRenderPending`
    // (a platform whose js-api predates it has none)
    sketch.on('render', () => this.onRendered?.next());
    sketch.on('change', ({source}) => {
      if (source === 'user')
        this.explicitMol = null;
      this.onChanged.next(null);
    });
    this.sketch = sketch;
    if (this.unshown !== null) {
      const {notation, value} = this.unshown;
      this.unshown = null;
      this.show(notation, value);
    }
  }

  get supportedExportFormats(): string[] {
    return ['smiles', 'mol', 'molV3000', 'smarts'];
  }

  /** No minimum size: Crux lays itself out in any box, and a dialog it is in can be shrunk to anything. Its size
   * as a dialog opens, 500 x 400, is in css/crux-sketcher.css. */
  get width(): number {
    return 0;
  }

  get height(): number {
    return 0;
  }

  get isInitialized(): boolean {
    return this.sketch !== null;
  }

  /** Crux Sketch's own handle once the sketcher is up (null before, and when it did not start): what the Sketcher app
   * reads the drawing's other formats and pictures from (CXSMILES, InChI, CDXML, `svgWith`, `png`), which the
   * platform's sketcher API does not give. */
  get handle(): SketcherHandle | null {
    return this.sketch;
  }

  get smiles(): string {
    return this.read('smiles');
  }

  set smiles(value: string) {
    this.write('smiles', value);
  }

  get molFile(): string {
    return this.read('molblock');
  }

  set molFile(value: string) {
    this.write('molblock', value);
  }

  get molV3000(): string {
    return this.read('molblockV3000');
  }

  set molV3000(value: string) {
    this.write('molblockV3000', value);
  }

  async getSmarts(): Promise<string> {
    return this.read('smarts');
  }

  set smarts(value: string) {
    this.write('smarts', value);
  }

  detach(): void {
    this.sketch?.destroy();
    super.detach();
  }

  // What the platform's tests read of Crux (`getWidgetStatus()`, the BDD library's areas and readings), in px of this
  // sketcher's root (the widget's hit-area system):
  // - a hit area for every atom and bond Crux draws, `atom 0`, `bond 0` (by index, as its molfile numbers them), centred
  //   where a press lands on it (Crux's `positions`) and half a bond wide;
  // - its parts (HOST-044; Crux's `layout`): `canvas`, the toolbars as `actions` (top), `tools` (left), `elements` (right)
  //   and `templates` (bottom), and `label editor` while it is open; and `tool <name>` for each control a toolbar shows,
  //   named as Crux's toolbar configuration names it (`tool bond.single`, `tool ring.benzene`, `tool undo`) (HOST-045);
  // - the readings `ready`, `pending`, `smiles` and `smarts` (Crux's own, whatever the host was given), `atoms`, `bonds`, `mode`,
  //   `selected atoms` and `selected bonds` (their indices in molfile order, `0, 2`, empty for none), `tool` (Crux's),
  //   `query` (whether the drawing has a query feature), `empty` and `changes` (the host's change events since this
  //   sketcher was made) (HOST-045).
  private status(): {hitAreas: {[name: string]: Rect}, values: {[name: string]: string | number | boolean}} {
    const sketch = this.sketch;
    const element = this.root.querySelector<HTMLElement>(':scope > crux-sketch');
    if (sketch === null || element === null || this.isDetached)
      return {hitAreas: {}, values: {ready: false, pending: this.isRenderPending, changes: this.changes}};
    const {atoms, bonds, bondLength} = sketch.positions;
    const at = element.getBoundingClientRect();
    const root = this.root.getBoundingClientRect();
    const [ox, oy] = [at.left - root.left, at.top - root.top];
    const half = Math.max(4, bondLength / 4);
    const hitAreas: {[name: string]: Rect} = {};
    const add = (name: string, p: {x: number, y: number} | null) => {
      if (p !== null)
        hitAreas[name] = {x: ox + p.x - half, y: oy + p.y - half, width: 2 * half, height: 2 * half};
    };
    atoms.forEach((p, i) => add(`atom ${i}`, p));
    bonds.forEach((p, i) => add(`bond ${i}`, p));
    const box = (name: string, b: Rect | null | undefined) => {
      if (b)
        hitAreas[name] = {x: ox + b.x, y: oy + b.y, width: b.width, height: b.height};
    };
    const layout = sketch.layout;
    box('canvas', layout.canvas);
    for (const [place, part] of Object.entries(PARTS))
      box(part, layout.toolbars[place as keyof typeof layout.toolbars]);
    box('label editor', layout.labelEditor);
    for (const [name, b] of Object.entries(layout.controls))
      box(`tool ${name}`, b);
    let smiles = '';
    try {
      smiles = sketch.smiles;
    } catch {
      // a drawing no SMILES holds (a query): its reading is empty
    }
    let smarts = '';
    try {
      smarts = sketch.smarts;
    } catch {
      // the engine cannot write it: the reading is empty
    }
    const selection = sketch.selection;
    return {hitAreas, values: {ready: true, pending: sketch.isPending, smiles, smarts, atoms: atoms.length, bonds: bonds.length,
      mode: sketch.mode, 'selected atoms': selection.atoms.join(', '), 'selected bonds': selection.bonds.join(', '),
      tool: sketch.tool, query: sketch.hasQuery, empty: sketch.isEmpty, changes: this.changes}};
  }

  /** Why the last read of one of `notations` threw, an engine failure first (any notation); null when each of them
   * last answered. What a host about to write one of them must not write (the engine-poison foot-gun). Each is
   * forgotten when its getter answers again (the engine recovers from a trap by itself: crux-sketch API-069); a start
   * on a failed engine is kept while the sketcher is open. */
  readFailure(notations: readonly Notation[]): CruxReadFailure | null {
    const all = [...this.failures.values()];
    return this.engineStop ?? all.find((f) => f.engine) ?? all.find((f) => notations.includes(f.notation)) ?? null;
  }

  private read(notation: Notation): string {
    if (this.explicitMol?.notation === notation)
      return this.explicitMol.value;
    if (this.unshown !== null) {
      const {notation: from, value} = this.unshown;
      return isEmptyValue(value) ? emptyValue(notation) :
        DG.chem.convert(value, DG_NOTATION[from], DG_NOTATION[notation]);
    }
    const sketch = this.sketch;
    if (sketch === null)
      return emptyValue(notation);
    try {
      let value = '';
      switch (notation) {
      case 'smiles': {
        const smiles = sketch.smiles;
        value = holdsQuery(sketch) ? DG.chem.smilesFromSmartsWarning() : smiles;
        break;
      }
      case 'molblock':
        value = sketch.molfile || emptyValue(notation);
        break;
      case 'molblockV3000':
        value = sketch.molV3000;
        break;
      case 'smarts':
        value = sketch.smarts;
        break;
      }
      this.failures.delete(notation);
      return value;
    } catch (e) {
      // the engine cannot write this drawing in that notation (a SketchError, never thrown for an empty drawing), or
      // it failed on this read twice, on a fresh engine too (crux-sketch API-069: that read alone). The platform's
      // getters answer '', and the failure is kept for the hosts that write (`readFailure`) until the getter answers
      // again: the cell editor writes nothing meanwhile
      const failure = readFailureOf(e, notation);
      this.failures.set(notation, failure);
      console.warn(`Crux Sketch: ${failure.message}`);
      return '';
    }
  }

  private write(notation: Notation, value: string): void {
    // A molblock given as SMILES is read as the molblock it is, as Ketcher and OpenChemLib read it. Core's Sketch action
    // (the context panel's Actions on a molecule) sets the current molecule's value as SMILES, and once a sketcher has
    // synced the current object (js-api's Sketcher.syncCurrentObject) that value is a molblock: Crux opened empty.
    if (notation === 'smiles' && DG.chem.isMolBlock(value))
      notation = value.includes('V3000') ? 'molblockV3000' : 'molblock';
    // and the other way round: a SMILES given as a molblock (SequenceTranslator's Markush enumerator sets its cores'
    // SMILES with setMolFile) is read as the SMILES it is
    else if ((notation === 'molblock' || notation === 'molblockV3000') && value && !DG.chem.isMolBlock(value))
      notation = 'smiles';
    // A SMARTS that also reads as SMILES (`[NH2]c1ccccc1`, which js-api's substructure host sets as SMARTS for its bracket
    // H count; run20 STAND S-004) is not kept as the caller's string: the substructure filter tells SMARTS from SMILES by
    // the text, took it for a molecule and dropped the H count (264 of mol1K's rows where it finds 11). The drawing's own
    // query answers instead.
    this.explicitMol = notation === 'smarts' && !DG.chem.isSmarts(value) ? null : {notation, value};
    if (this.sketch === null || this.isDetached) {
      this.unshown = {notation, value};
      if (this.isDetached)
        this.onChanged.next(null);
      return;
    }
    this.show(notation, value);
  }

  // Replaces the drawing with one write, and so one change event. A value Crux cannot read changes nothing and
  // says nothing, as with the other sketchers: the host has shown it malformed, and keeps showing it so.
  private show(notation: Notation, value: string): void {
    if (isEmptyValue(value))
      this.sketch!.setValue('');
    else
      this.sketch!.setValue(value, {format: CRUX_FORMAT[notation]});
  }
}

/** Whether Crux's drawing is a query its SMILES can't write (Crux warns so). */
function holdsQuery(sketch: SketcherHandle): boolean {
  return sketch.warnings.some((w) => w.kind === 'format' && w.format === 'smiles' && w.code === 'query');
}

function emptyValue(notation: Notation): string {
  return notation === 'molblock' ? DG.WHITE_MOLBLOCK : '';
}
