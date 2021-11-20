export enum WORKER_CALL {
  INIT_TANIMOTO_FINGERPRINTS = 'initTanimotoFingerprints',
  FREE_TANIMOTO_FINGERPRINTS = 'freeTanimotoFingerprints',
  INIT_MOLECULES_STRUCTURES = 'initRdKitMolecules',
  FREE_MOLECULES_STRUCTURES = 'freeRdKitMolecules',
  SEARCH_SUBSTRUCTURE = 'searchSubstructure',
  GET_TANIMOTO_FINGERPRINTS = 'getTanimotoFingerprints',
  GET_SIMILARITIES = 'getSimilarities',
  INIT_STRUCTURAL_ALERTS = 'initStructuralAlerts',
  GET_STRUCTURAL_ALERTS = 'getStructuralAlerts'
}