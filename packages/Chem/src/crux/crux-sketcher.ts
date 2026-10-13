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
function loadCruxSketch(): Promise<typeof CruxSketch> {
  return import('../../vendor/crux-sketch/crux-sketch.js');
}

/** The name Crux's hit areas and readings go under in the sketcher's status (`getWidgetStatus()`). */
const STATUS = 'crux-sketch';

/** The sketcher's u2 name (`data-u2-name` on its root, HOST-043): what the platform's automation finds it by. */
const U2_NAME = 'cruxSketch';

/** Crux's toolbars by their place, as the u2 parts the status names them (HOST-044): the top bar holds the actions (undo,
 * clear, zoom, layout, the clipboard), the left one the drawing tools, the right one the elements, the bottom one the
 * templates (the rings and the Structure Library). */
const PARTS: {[place: string]: string} = {top: 'actions', left: 'tools', right: 'elements', bottom: 'templates'};

type Rect = {x: number, y: number, width: number, height: number};

function isEmptyValue(value: string | null | undefined): boolean {
  return !value?.trim() || grok.chem.Sketcher.isEmptyMolfile(value);
}

/** Crux Sketch as a Datagrok molecule sketcher. Crux's getters are synchronous and stay readable after
 * `destroy()`, and its `change` event says whether a change is the user's edit or a value written to it, so
 * the host's rules map onto it directly: a setter is one write and one `onChanged`, a user edit is one
 * `onChanged` and ends `explicitMol`. Its root is named `cruxSketch` (`data-u2-name`), and its status
 * (`getWidgetStatus()`) shows the platform's tests what Crux draws and holds: a hit area for every atom and bond, its
 * parts and its controls, readings, and `isRenderPending` and `onRendered` to settle on. */
export class CruxSketcher extends grok.chem.SketcherBase {
  private sketch: SketcherHandle | null = null;
  /** A value set while no live Crux could show it: before `init` resolved, or after `detach`. */
  private unshown: {notation: Notation, value: string} | null = null;
  /** Crux is being created: from `init` until it is ready or failed to start. */
  private loading = false;
  /** The host's change events since this sketcher was made (`onChanged`): the status's `changes` reading. */
  private changes = 0;

  constructor() {
    super();
    this.root.classList.add('crux-sketcher');
    this.name = U2_NAME;
    DG.UndoService.ownScope(this.root);
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

  async init(host: DG.chem.Sketcher): Promise<void> {
    this.host = host;
    let sketch: SketcherHandle;
    this.loading = true;
    try {
      // Where Datagrok asks for a query (its substructure filter), Crux opens in query mode, offering the query tools;
      // every other host opens it in molecule mode. The user can switch it in Crux's gear menu either way.
      sketch = await (await loadCruxSketch()).createSketcher(this.root,
        host.isSubstructureFilter ? {mode: 'query'} : undefined);
    } catch (e) {
      console.error(e);
      this.root.append(ui.divText(`Crux Sketch could not start: ${(e as {message?: string})?.message ?? e}`,
        'crux-sketcher-error'));
      return;
    } finally {
      this.loading = false;
    }
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

  /**
   * What the platform's tests read of Crux (`getWidgetStatus()`, the BDD library's areas and readings), in px of this
   * sketcher's root (the widget's hit-area system):
   * - a hit area for every atom and bond Crux draws, `atom 0`, `bond 0` (by index, as its molfile numbers them), centred
   *   where a press lands on it (Crux's `positions`) and half a bond wide;
   * - its parts (HOST-044; Crux's `layout`): `canvas`, the toolbars as `actions` (top), `tools` (left), `elements` (right)
   *   and `templates` (bottom), and `label editor` while it is open; and `tool <name>` for each control a toolbar shows,
   *   named as Crux's toolbar configuration names it (`tool bond.single`, `tool ring.benzene`, `tool undo`) (HOST-045);
   * - the readings `ready`, `pending`, `smiles` and `smarts` (Crux's own, whatever the host was given), `atoms`, `bonds`, `mode`,
   *   `selected atoms` and `selected bonds` (their indices in molfile order, `0, 2`, empty for none), `tool` (Crux's),
   *   `query` (whether the drawing has a query feature), `empty` and `changes` (the host's change events since this
   *   sketcher was made) (HOST-045).
   */
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
      switch (notation) {
      case 'smiles': {
        const smiles = sketch.smiles;
        return sketch.warnings.some((w) => w.kind === 'format' && w.format === 'smiles' && w.code === 'query') ?
          DG.chem.smilesFromSmartsWarning() : smiles;
      }
      case 'molblock':
        return sketch.molfile || emptyValue(notation);
      case 'molblockV3000':
        return sketch.molV3000;
      case 'smarts':
        return sketch.smarts;
      }
    } catch (e) {
      // the engine cannot write this drawing in that notation (a SketchError, never thrown for an empty drawing)
      console.warn(`Crux Sketch: ${(e as {message?: string})?.message ?? e}`);
      return '';
    }
  }

  private write(notation: Notation, value: string): void {
    this.explicitMol = {notation, value};
    if (this.sketch === null || this.isDetached) {
      this.unshown = {notation, value};
      if (this.isDetached)
        this.onChanged.next(null);
      return;
    }
    this.show(notation, value);
  }

  /** Replaces the drawing with one write, and so one change event. A value Crux cannot read changes nothing and
   * says nothing, as with the other sketchers: the host has shown it malformed, and keeps showing it so. */
  private show(notation: Notation, value: string): void {
    if (isEmptyValue(value))
      this.sketch!.setValue('');
    else
      this.sketch!.setValue(value, {format: CRUX_FORMAT[notation]});
  }
}

function emptyValue(notation: Notation): string {
  return notation === 'molblock' ? DG.WHITE_MOLBLOCK : '';
}
