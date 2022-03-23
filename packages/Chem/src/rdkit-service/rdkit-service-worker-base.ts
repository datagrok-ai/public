import {RDModule, RDMol} from "../rdkit-api";

export class RdKitServiceWorkerBase {
  _rdKitModule: RDModule;
  _webRoot: string;
  _rdKitMols: RDMol[] | null = null;
  constructor(module: RDModule, webRoot: string) {
    this._rdKitModule = module;
    this._webRoot = webRoot;
  }
}