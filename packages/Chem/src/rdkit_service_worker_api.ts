export enum WORKER_CALL {
  INIT_MORGAN_FINGERPRINTS = 'initMorganFingerprints',
  FREE_MORGAN_FINGERPRINTS = 'freeMorganFingerprints',
  GET_MORGAN_FINGERPRINTS = 'getTanimotoFingerprints',
  INIT_MOLECULES_STRUCTURES = 'initRdKitMolecules',
  FREE_MOLECULES_STRUCTURES = 'freeRdKitMolecules',
  SEARCH_SUBSTRUCTURE = 'searchSubstructure',
  GET_SIMILARITIES = 'getSimilarities',
  INIT_STRUCTURAL_ALERTS = 'initStructuralAlerts',
  GET_STRUCTURAL_ALERTS = 'getStructuralAlerts'
}