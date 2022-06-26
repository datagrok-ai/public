import * as grok from 'datagrok-api/grok';
import * as OCL from 'openchemlib/full';

let sketcherId = 0;

export class OpenChemLibSketcher extends grok.chem.SketcherBase {
  // @ts-ignore
  _sketcher: OCL.StructureEditor;

  constructor() {
    super();
  }

  async init(): Promise<void> {
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

  async getSmarts(): Promise<string> {
    return this.smiles;
  }

  setSmarts(s: string): void {
    this.smiles = s;
  }

  detach(): void {
    console.log('OCL sketcher detached');
    super.detach();
  }
}
