import {RdKitServiceWorkerSubstructure} from './rdkit-service-worker-substructure';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

/**
 * RDKit's worker code, written against RDKit.js's module: the RDKit pool's workers run it on RDKit (rdkit.worker.ts),
 * the Crux workers on the Crux module (crux/crux.worker.ts), so it stays when RDKit goes (engine/README.md).
 */
export class RdKitServiceWorker extends RdKitServiceWorkerSubstructure {
  constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }
}
