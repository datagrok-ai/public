import * as grok from 'datagrok-api/grok';
import {chem} from 'datagrok-api/grok';
import * as OCL from 'openchemlib/full';
import {getRdKitModule} from '../utils/chem-common-rdkit';

let sketcherId = 0;

export class OpenChemLibSketcher extends grok.chem.SketcherBase {
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
    return this._sketcher ? this._sketcher.getSmiles() : this.host?.getSmiles();
  }

  set smiles(s) {
    this._sketcher.setSmiles(s);
  }

  get molFile() {
    return this._sketcher ? this._sketcher.getMolFile() : this.host?.getMolFile();
  }

  set molFile(s) {
    this._sketcher.setMolFile(s);
  }

  get molV3000() {
    return this._sketcher ? this._sketcher.getMolFileV3() : this.host?.getMolFile();
  }

  set molV3000(s) {
    this._sketcher.setMolFile(s);
  }

  async getSmarts(): Promise<string> {
    if (this._sketcher) {
      const mol = getRdKitModule().get_mol(this.molFile);
      const smarts = mol.get_smarts();
      mol?.delete();
      return smarts;
    } else
      return await this.host!.getSmarts() as string;
  }

  set smarts(s: string) {
    this.convertAndSetSmarts(s);
  }

  detach(): void {
    console.log('OCL sketcher detached');
    super.detach();
  }

  private async convertAndSetSmarts(s: string) {
    const mol = getRdKitModule().get_mol(s);
    this._sketcher?.setMolFile(mol.get_molblock());
    mol?.delete();
  }
}
