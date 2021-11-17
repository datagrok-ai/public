export class RdKitServiceWorkerBase {
  rdKitModule: any | null;
  webRoot: string | null;
  constructor(module: Object, webRoot: string) {
    this.rdKitModule = module;
    this.webRoot = webRoot;
  }
}