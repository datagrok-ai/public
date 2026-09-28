import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

/**
 * Initializes the RDKit build from a module compiled once, on the main thread. Emscripten reports an asynchronous
 * `instantiateWasm` failure nowhere, so the instantiation's rejection is passed on here.
 */
export function initRdKitFrom(init: (options: object) => Promise<RDModule>,
  wasm: WebAssembly.Module): Promise<RDModule> {
  return new Promise((resolve, reject) => {
    init({instantiateWasm: (imports: WebAssembly.Imports,
      receive: (instance: WebAssembly.Instance, module: WebAssembly.Module) => void) => {
      WebAssembly.instantiate(wasm, imports).then((instance) => receive(instance, wasm), reject);
      return {};
    }}).then(resolve, reject);
  });
}

export enum WORKER_CALL {
  INIT_MOLECULES_STRUCTURES = 'initRdKitMolecules',
  FREE_MOLECULES_STRUCTURES = 'freeRdKitMolecules',
  GET_FINGERPRINTS = 'getFingerprints',
  SEARCH_SUBSTRUCTURE = 'searchSubstructure',
  CONVERT_MOL_NOTATION = 'convertMolNotation',
  FLATTEN_MOLECULES = 'flattenMolecules',
  GET_INCHIS = 'getInchis',
  GET_STRUCTURAL_ALERTS = 'getStructuralAlerts',
  INVALIDATE_CACHE = 'invalidateCache',
  SET_TERMINATE_FLAG = 'setTerminateFlag',
  MMP_GET_FRAGMENTS = 'mmpGetFragments',
  MMP_LINK_FRAGMENTS = 'mmpLinkFragments',
  LINK_R_GROUP_FRAGMENTS = 'linkRGroupFragments',
  MMP_GET_MCS = 'mmpGetMcs',
  MOST_COMMON_STRUCTURE = 'mostCommonStructure',
  R_GROUP_ANALYSIS = 'rGroupAnalysis',
  BEAUTIFY_MOLS = 'beautifyMols',
  GET_COORDGEN_COORDS = 'getCoordGenCoords',
}
