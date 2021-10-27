import {createRDKit} from './RDKit_minimal_2021.03_17.js'
import {RdKitSubstructLibrary} from './rdkit_substruct_library.js'

// var _rdKitModule = null;
// var _substructLibrary = null;

async function handler(e) {
  const {op, args} = e.data;
  let port = e.ports[0];
  if (op === 'module::init') {
    const webRoot = args[0];
    handler._rdKitModule = await createRDKit(webRoot);
    console.log("RDKit (worker) initialized");
    handler._substructLibrary = new RdKitSubstructLibrary(handler._rdKitModule);
    port.postMessage({op: op});
  } else if (op === 'substructLibrary::init') {
    const result = handler._substructLibrary.init(args[0]);
    port.postMessage({op: op, retval: result});
  } else if (op === 'substructLibrary::search') {
    const result = handler._substructLibrary.search(args[0], args[1]);
    port.postMessage({op: op, retval: result});
  } else if (op === 'substructLibrary::deinit') {
    handler._substructLibrary.deinit();
    handler._substructLibrary = null;
    port.postMessage({op: op});
  }
}

onmessage = handler;
