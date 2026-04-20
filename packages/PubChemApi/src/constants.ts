export const pubChemBaseURL = 'https://pubchem.ncbi.nlm.nih.gov';
export const pubChemRest = `${pubChemBaseURL}/rest`;
export const pubChemPug = `${pubChemRest}/pug`;

export enum COLUMN_NAMES {
  CANONICAL_SMILES = 'CanonicalSMILES',
  CONNECTIVITY_SMILES = 'ConnectivitySMILES',
  CID = 'CID',
  SCORE = 'score',
  INDEX = 'index',
}