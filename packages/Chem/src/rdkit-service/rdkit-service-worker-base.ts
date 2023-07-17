import {RDModule, RDMol, RDSubstructLibrary} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import BitArray from '@datagrok-libraries/utils/src/bit-array';

export class RdKitServiceWorkerBase {
  _rdKitModule: RDModule;
  _webRoot: string;
  _rdKitMols: (RDMol | null)[] | null = null;
  _substructLibrary: RDSubstructLibrary | null = null;
  _malformedIdxs: BitArray | null = null;
  _terminateFlag: boolean = false;
  constructor(module: RDModule, webRoot: string) {
    this._rdKitModule = module;
    this._webRoot = webRoot;
  }
}
