import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';

let sketcherId = 0;

export class OpenChemLibSketcher extends DG.chem.SketcherBase {
  // @ts-ignore
  _sketcher: OCL.StructureEditor;

  constructor() {
    super();
  }

  async init() {
    let id = `ocl-sketcher-${sketcherId++}`;
    this.root.id = id;

    this._sketcher = OCL.StructureEditor.createSVGEditor(id, 1);
    this._sketcher.setChangeListenerCallback((id, molecule) => {
      this.onChanged.next(null);
    });
  }

  get smiles() { return this._sketcher.getSmiles(); }
  set smiles(s) { this._sketcher.setSmiles(s); }

  get molFile() { return this._sketcher.getMolFile(); }
  set molFile(s) { this._sketcher.setMolFile(s); }
}
