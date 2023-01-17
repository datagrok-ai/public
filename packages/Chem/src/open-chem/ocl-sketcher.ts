import * as grok from 'datagrok-api/grok';
import {chem} from 'datagrok-api/grok';
import * as OCL from 'openchemlib/full';
import { convertMolNotation } from '../package';
import * as DG from 'datagrok-api/dg';

let sketcherId = 0;

export class OpenChemLibSketcher extends grok.chem.SketcherBase {

  _sketcher: OCL.StructureEditor | null = null;

  constructor() {
    super();
  }

  async init(host: chem.Sketcher): Promise<void> {
    this.host = host;
    const id = `ocl-sketcher-${sketcherId++}`;
    this.root.id = id;
    this._sketcher = OCL.StructureEditor.createSVGEditor(id, 1);
    this._sketcher.setChangeListenerCallback((id: any, molecule: any) => {
      this.onChanged.next(null);
    });
  }

  get supportedExportFormats() {
    return ['smiles', 'mol'];
  }

  get smiles() {
    return this._sketcher ? this._sketcher.getSmiles() : '';
  }

  set smiles(s) {
    this._sketcher?.setSmiles(s);
  }

  get molFile() {
      return this._sketcher ? this._sketcher.getMolFile() : '';
  }

  set molFile(s) {
    this._sketcher?.setMolFile(s);
  }

  get molV3000() {
      return this._sketcher ? this._sketcher.getMolFileV3() : '';
  }

  set molV3000(s) {
    this._sketcher?.setMolFile(s);
  }

  async getSmarts(): Promise<string> {
    return convertMolNotation(this.molFile, DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts);
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
    const molfile = convertMolNotation(s, DG.chem.Notation.Smarts, DG.chem.Notation.MolBlock);
    this._sketcher?.setMolFile(molfile);
  }
}
