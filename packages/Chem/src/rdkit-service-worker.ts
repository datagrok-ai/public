import {RdKitServiceWorkerSubstructure} from './rdkit-service-worker-substructure';

export class RdKitServiceWorker extends RdKitServiceWorkerSubstructure {
  constructor(module: Object, webRoot: string) {
    super(module, webRoot);
  }
}
