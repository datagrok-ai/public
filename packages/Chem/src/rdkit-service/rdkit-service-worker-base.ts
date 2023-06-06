import {RDModule, RDMol, RDSubstructLibrary} from '@datagrok-libraries/chem-meta/src/rdkit-api';

export class RdKitServiceWorkerBase {
  _rdKitModule: RDModule;
  _webRoot: string;
  _rdKitMols: (RDMol | null)[] | null = null;
  _substructLibrary: RDSubstructLibrary | null = null;
  constructor(module: RDModule, webRoot: string) {
    this._rdKitModule = module;
    this._webRoot = webRoot;
  }
}
