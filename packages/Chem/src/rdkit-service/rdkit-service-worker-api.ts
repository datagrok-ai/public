export enum WORKER_CALL {
  INIT_MOLECULES_STRUCTURES = 'initRdKitMolecules',
  FREE_MOLECULES_STRUCTURES = 'freeRdKitMolecules',
  GET_FINGERPRINTS = 'getFingerprints',
  SEARCH_SUBSTRUCTURE = 'searchSubstructure',
  CONVERT_MOL_NOTATION = 'convertMolNotation',
  GET_STRUCTURAL_ALERTS = 'getStructuralAlerts',
  INVALIDATE_CACHE = 'invalidateCache',
  SET_TERMINATE_FLAG = 'setTerminateFlag',
  MOST_COMMON_STRUCTURE = 'mostCommonStructure',
}
