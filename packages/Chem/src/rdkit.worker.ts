import {WORKER_CALL} from './rdkit_service_worker_api';
import {RdKitServiceWorker as ServiceWorkerClass} from './rdkit_service_worker';
// @ts-ignore
import {createRDKit} from './RDKit_minimal_2021.03_18.js';

const ctx: Worker = self as any;
let _rdKitModule: any | null = null;
let _rdkitServiceWorker: ServiceWorkerClass | null = null;

ctx.addEventListener("message", async (e: any) => {
  const {op, args} = e.data;
  let port = e.ports[0];
  if (op === 'module::init') {
    const webRoot = args[0];
    _rdKitModule = await createRDKit(webRoot);
    console.log("RDKit (worker) initialized");
    _rdkitServiceWorker = new ServiceWorkerClass(_rdKitModule, webRoot);
    port.postMessage({op: op, retval: null});
  } else if (op === WORKER_CALL.INIT_MOLECULES_STRUCTURES) {
    const result = _rdkitServiceWorker!.initMoleculesStructures(args[0]);
    port.postMessage({op: op, retval: result});
  } else if (op === WORKER_CALL.SEARCH_SUBSTRUCTURE) {
    const result = _rdkitServiceWorker!.searchSubstructure(args[0], args[1]);
    port.postMessage({op: op, retval: result});
  } else if (op === WORKER_CALL.FREE_MOLECULES_STRUCTURES) {
    _rdkitServiceWorker!.freeMoleculesStructures();
    _rdkitServiceWorker = null;
    port.postMessage({op: op, retval: null});
  } else if (op === WORKER_CALL.INIT_MORGAN_FINGERPRINTS) {
    _rdkitServiceWorker!.initMorganFingerprints();
    port.postMessage({op: op, retval: null});
  } else if (op === WORKER_CALL.GET_MORGAN_FINGERPRINTS) {
    const result = _rdkitServiceWorker!.getMorganFingerprints();
    port.postMessage({op: op, retval: result});
  }
});