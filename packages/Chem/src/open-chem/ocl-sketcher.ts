import * as grok from 'datagrok-api/grok';
import {chem} from 'datagrok-api/grok';
import * as OCL from 'openchemlib/full';
import {PackageFunctions} from '../package';
import * as DG from 'datagrok-api/dg';

let sketcherId = 0;

export class OpenChemLibSketcher extends grok.chem.SketcherBase {
  _sketcher: OCL.StructureEditor | null = null;
  mouseDown = false;
  savedMolecule: string | null = null;
  setSavedMolecule = false;

  constructor() {
    super();
  }

  async init(host: chem.Sketcher): Promise<void> {
    this.host = host;
    const id = `ocl-sketcher-${sketcherId++}`;
    this.root.id = id;
    this._sketcher = OCL.StructureEditor.createSVGEditor(id, 1);
    this._sketcher.setChangeListenerCallback((_) => {
      this.explicitMol = null;
      if (this.setSavedMolecule) {
        this.setSavedMolecule = false;
        this.molFile = this.savedMolecule!;
        this.savedMolecule = null;
      } else
        this.onChanged.next(null);
    });

    /* workaround for situations when mouse was down outside sketcher editor but up inside sketcher editor
    (in this case the action selected on the toolbox is applyed unintendedly. For instance, molecule is removed).
    stopImmediatePropagation doesn't prevent OCL editor from deleting molecule.
    Using pointerup instead of mouseup here since at the moment of mouseup event molecule is already removed
    TODO!!!:  create issue on openChemLib github*/

    const canvas = this.root.querySelectorAll('canvas')[1];
    canvas!.addEventListener('pointerdown', (_) => {
      this.mouseDown = true;
    });
    canvas!.addEventListener('pointerup', (_) => {
      if (this.mouseDown)
        this.mouseDown = false;
      else {
        this.setSavedMolecule = true;
        this.savedMolecule = this.molFile;
      }
    });
  }

  get supportedExportFormats() {
    return ['smiles', 'mol'];
  }

  get smiles() {
    if (this.explicitMol?.notation === 'smiles')
      return this.explicitMol.value;
    return this._sketcher ? this._sketcher.getSmiles() : '';
  }

  set smiles(s) {
    this._sketcher?.setSmiles(s);
    this.explicitMol = {notation: 'smiles', value: s};
  }

  get molFile() {
    if (this.explicitMol?.notation === 'molblock')
      return this.explicitMol.value;
    return this._sketcher ? this._sketcher.getMolFile() : '';
  }

  set molFile(s) {
    this._sketcher?.setMolFile(s);
    this.explicitMol = {notation: 'molblock', value: s};
  }

  get molV3000() {
    if (this.explicitMol?.notation === 'molblockV3000')
      return this.explicitMol.value;
    return this._sketcher ? this._sketcher.getMolFileV3() : '';
  }

  set molV3000(s) {
    this._sketcher?.setMolFile(s);
    this.explicitMol = {notation: 'molblockV3000', value: s};
  }

  async getSmarts(): Promise<string> {
    return PackageFunctions.convertMolNotation(this.molFile, DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts);
  }

  set smarts(s: string) {
    this.convertAndSetSmarts(s);
  }

  get isInitialized() {
    return this._sketcher !== null;
  }

  detach(): void {
    console.log('OCL sketcher detached');
    super.detach();
  }

  private async convertAndSetSmarts(s: string) {
    const molfile = PackageFunctions.convertMolNotation(s, DG.chem.Notation.Smarts, DG.chem.Notation.MolBlock);
    this._sketcher?.setMolFile(molfile);
  }
}
