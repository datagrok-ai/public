import {createRDKit} from './RDKit_minimal_2021.03_18.js';
import {RdKitServiceWorker as ServiceClass} from './rdkit_service_worker';
// import {RdKitFingerprintSearcher as SearcherClass} from './rdkit_fingerprint_searcher';

const ctx: Worker = self as any;

let handler: any = {}; // Own for each worker

ctx.addEventListener("message", async (e: any) => {
  const {op, args} = e.data;
  let port = e.ports[0];
  if (op === 'module::init') {
    const webRoot = args[0];
    handler._rdKitModule = await createRDKit(webRoot);
    console.log("RDKit (worker) initialized");
    handler._rdkitServiceWorker = new ServiceClass(handler._rdKitModule, webRoot);
    port.postMessage({op: op});
  } else if (op === 'substructLibrary::init') {
    const result = handler._rdkitServiceWorker.init(args[0]);
    port.postMessage({op: op, retval: result});
  } else if (op === 'substructLibrary::search') {
    const result = handler._rdkitServiceWorker.search(args[0], args[1]);
    port.postMessage({op: op, retval: result});
  } else if (op === 'substructLibrary::deinit') {
    handler._rdkitServiceWorker.deinit();
    handler._rdkitServiceWorker = null;
    port.postMessage({op: op});
  } else if (op === 'structuralAlerts::init') {
    const result = handler._rdkitServiceWorker.initStructuralAlerts(args[0]);
    port.postMessage({op: op});
  } else if (op === 'structuralAlerts::get') {
    const result = handler._rdkitServiceWorker.getStructuralAlerts(args[0]);
    port.postMessage({op: op, retval: result});
  }
});