import {RdKitServiceWorkerSimilarity} from './rdkit_service_worker_similarity';
import { loadAlertsCollection, getStructuralAlerts, setStructuralAlertsRdKitModule } from './widgets/structural-alerts';

export class RdKitServiceWorkerMisc extends RdKitServiceWorkerSimilarity {

  constructor(module: Object, webRoot: string) {
    super(module, webRoot);
    setStructuralAlertsRdKitModule(module, webRoot);
  }

  initStructuralAlerts(smarts: string[]) {
    loadAlertsCollection(smarts);
  }

  getStructuralAlerts(smiles: string) {
    return getStructuralAlerts(smiles);
  }

}