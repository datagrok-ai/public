import {RdKitServiceWorkerSimilarity} from './rdkit_service_worker_similarity';

export class RdKitServiceWorker extends RdKitServiceWorkerSimilarity {
  constructor(module: Object, webRoot: string) {
    super(module, webRoot);
  }
}