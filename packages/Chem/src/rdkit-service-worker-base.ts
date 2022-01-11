export class RdKitServiceWorkerBase {
  _rdKitModule: any | null;
  _webRoot: string | null;
  constructor(module: Object, webRoot: string) {
    this._rdKitModule = module;
    this._webRoot = webRoot;
  }
}
