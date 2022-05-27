export enum WORKER_CALL {
  INIT_MOLECULES_STRUCTURES = 'initRdKitMolecules',
  FREE_MOLECULES_STRUCTURES = 'freeRdKitMolecules',
  GET_MORGAN_FINGERPRINTS = 'getMorganFingerprints',
  GET_PATTERN_FINGERPRINTS = 'getPatternFingerprints',
  SEARCH_SUBSTRUCTURE = 'searchSubstructure',
}
