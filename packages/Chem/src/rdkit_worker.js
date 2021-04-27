importScripts('RDKit_minimal_2021.03_13.js');
importScripts('rdkit_substruct_library.js');

let _rdKitModule = null;
let _substructLibrary = null;

onmessage = async function (e) {

  const {op, args} = e.data;
  let port = e.ports[0];

  if (op === 'module::init') {
    _rdKitModule = await initRDKitModule();
    console.log("RDKit (worker) initialized");
    _substructLibrary = new RdKitSubstructLibrary(_rdKitModule);
    port.postMessage({op: op});
  } else if (op === 'substructLibrary::init') {
    _substructLibrary.init(args[0]);
    port.postMessage({op: op});
  } else if (op === 'substructLibrary::search') {
    const result = _substructLibrary.search(args[0]);
    port.postMessage({op: op, retval: result});
  } else if (op === 'substructLibrary::deinit') {
    _substructLibrary.deinit();
    _substructLibrary = null;
    port.postMessage({op: op});
  }

}
