import {RdkitServiceWorkerSubstructure} from './rdkit_service_worker_substructure';

export class RdKitServiceWorker extends RdkitServiceWorkerSubstructure {
  constructor(module: Object, webRoot: string) {
    super(module, webRoot);
  }
}
