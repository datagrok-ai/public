import {WORKER_CALL} from './rdkit_service_worker_api';
import {RdKitServiceWorker as ServiceWorkerClass} from './rdkit_service_worker';
// @ts-ignore
import {createRDKit} from './RDKit_minimal_2021.03_18.js';

const ctx: Worker = self as any;

let handler: any = {}; // Own for each worker

ctx.addEventListener("message", async (e: any) => {
  const {op, args} = e.data;
  let port = e.ports[0];
  if (op === 'module::init') {
    const webRoot = args[0];
    handler._rdKitModule = await createRDKit(webRoot);
    console.log("RDKit (worker) initialized");
    handler._rdkitServiceWorker = new ServiceWorkerClass(handler._rdKitModule, webRoot);
    port.postMessage({op: op, retval: null});
  } else if (op === WORKER_CALL.INIT_MOLECULES_STRUCTURES) {
    const result = handler._rdkitServiceWorker.initMoleculesStructures(args[0]);
    port.postMessage({op: op, retval: result});
  } else if (op === WORKER_CALL.SEARCH_SUBSTRUCTURE) {
    const result = handler._rdkitServiceWorker.searchSubstructure(args[0], args[1]);
    port.postMessage({op: op, retval: result});
  } else if (op === WORKER_CALL.FREE_MOLECULES_STRUCTURES) {
    handler._rdkitServiceWorker.freeMoleculesStructures();
    handler._rdkitServiceWorker = null;
    port.postMessage({op: op, retval: null});
  } else if (op === WORKER_CALL.INIT_MORGAN_FINGERPRINTS) {
    handler._rdkitServiceWorker.initMorganFingerprints();
    port.postMessage({op: op, retval: null});
  } else if (op === WORKER_CALL.GET_MORGAN_FINGERPRINTS) {
    const result = handler._rdkitServiceWorker.getMorganFingerprints();
    port.postMessage({op: op, retval: result});
  }
});