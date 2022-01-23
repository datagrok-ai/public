import {RDMol} from "../rdkit-api";

export class RdKitServiceWorkerBase {
  _rdKitModule: any;
  _webRoot: string;
  _rdKitMols: RDMol[] | null = null;
  constructor(module: Object, webRoot: string) {
    this._rdKitModule = module;
    this._webRoot = webRoot;
  }
}