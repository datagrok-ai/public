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

function isEmptyValue(value: string | null | undefined): boolean {
  return !value?.trim() || grok.chem.Sketcher.isEmptyMolfile(value);
}

/** Crux Sketch as a Datagrok molecule sketcher. Crux's getters are synchronous and stay readable after
 * `destroy()`, and its `change` event says whether a change is the user's edit or a value written to it, so
 * the host's rules map onto it directly: a setter is one write and one `onChanged`, a user edit is one
 * `onChanged` and ends `explicitMol`. */
export class CruxSketcher extends grok.chem.SketcherBase {
  private sketch: SketcherHandle | null = null;
  /** A value set while no live Crux could show it: before `init` resolved, or after `detach`. */
  private unshown: {notation: Notation, value: string} | null = null;

  constructor() {
    super();
    this.root.classList.add('crux-sketcher');
    DG.UndoService.ownScope(this.root);
  }

  async init(host: DG.chem.Sketcher): Promise<void> {
    this.host = host;
    let sketch: SketcherHandle;
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
    }
    if (this.isDetached) {
      sketch.destroy();
      return;
    }
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
