export enum WORKER_CALL {
  INIT_MOLECULES_STRUCTURES = 'initRdKitMolecules',
  FREE_MOLECULES_STRUCTURES = 'freeRdKitMolecules',
  INIT_MORGAN_FINGERPRINTS = 'initMorganFingerprints',
  GET_MORGAN_FINGERPRINTS = 'getMorganFingerprints',
  FREE_MORGAN_FINGERPRINTS = 'freeMorganFingerprints',
  SEARCH_SUBSTRUCTURE = 'searchSubstructure',
}
