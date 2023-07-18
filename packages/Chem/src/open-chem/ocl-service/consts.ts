export const enum OCLServiceCall {
  TOXICITY = 'toxicity',
  CHEM_PROPERTIES = 'chem_properties',
  DRUG_LIKENESS = 'drug_likeness',
};

export const riskTypes: {[index: number]: string} = {
  0: 'Mutagenicity',
  1: 'Tumorigenicity',
  2: 'Irritating effects',
  3: 'Reproductive effects',
} as const;

export const riskLevels: {[index: number]: string} = {
  0: 'Unknown',
  1: 'None',
  2: 'Low',
  3: 'High',
} as const;

export const enum MolNotationType {
  SMILES = 'smiles',
  MOLBLOCK = 'molblock',
}
