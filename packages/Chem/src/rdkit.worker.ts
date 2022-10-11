import {WORKER_CALL} from './rdkit-service/rdkit-service-worker-api';
import {RdKitServiceWorker as ServiceWorkerClass} from './rdkit-service/rdkit-service-worker';
// @ts-ignore
import initRDKitModule from './RDKit_minimal.js';
//@ts-ignore
import rdKitLibVersion from './rdkit_lib_version';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

const ctx: Worker = self as any;
let _rdKitModule: RDModule;
let _rdKitServiceWorker: ServiceWorkerClass | null = null;

ctx.addEventListener('message', async (e: any) => {
  const {op, args} = e.data;
  const port = e.ports[0];
  if (op === 'module::init') {
    const webRoot = args[0];
    _rdKitModule = await initRDKitModule({locateFile: () => `${webRoot}/dist/${rdKitLibVersion}.wasm`});
    console.log('RDKit (worker) initialized');
    _rdKitServiceWorker = new ServiceWorkerClass(_rdKitModule, webRoot);
    port.postMessage({op: op, retval: null});
  } else if (op === WORKER_CALL.INIT_MOLECULES_STRUCTURES) {
    const result = _rdKitServiceWorker!.initMoleculesStructures(args[0]);
    port.postMessage({op: op, retval: result});
  } else if (op === WORKER_CALL.SEARCH_SUBSTRUCTURE) {
    const result = _rdKitServiceWorker!.searchSubstructure(args[0], args[1], args[2]);
    port.postMessage({op: op, retval: result});
  } else if (op === WORKER_CALL.FREE_MOLECULES_STRUCTURES) {
    _rdKitServiceWorker!.freeMoleculesStructures();
    _rdKitServiceWorker = null;
    port.postMessage({op: op, retval: null});
  } else if (op === WORKER_CALL.GET_FINGERPRINTS) {
    const result = _rdKitServiceWorker!.getFingerprints(args[0]);
    port.postMessage({op: op, retval: result});
  }
});
