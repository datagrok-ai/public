export const SEQUENCE_TYPES = {
  SENSE_STRAND: 'SS',
  ANTISENSE_STRAND: 'AS',
  DUPLEX: 'Duplex',
  TRIPLEX: 'Triplex',
  DIMER: 'Dimer',
};

export const COL_NAMES = {
  CHEMISTRY: 'Chemistry',
  NUMBER: 'Number',
  TYPE: 'Type',
  CHEMISTRY_NAME: 'Chemistry Name',
  INTERNAL_COMPOUND_ID: 'Internal compound ID',
  IDP: 'IDP',
  SEQUENCE: 'Sequence',
  COMPOUND_NAME: 'Compound Name',
  COMPOUND_COMMENTS: 'Compound Comments',
  SALT: 'Salt',
  EQUIVALENTS: 'Equivalents',
  PURITY: 'Purity',
  COMPOUND_MOL_WEIGHT: 'Cpd MW',
  SALT_MOL_WEIGHT: 'Salt MW',
  SALT_MASS: 'Salt mass',
  BATCH_MOL_WEIGHT: 'Batch MW',
  SOURCE: 'Source',
  ICD: 'ICD',
  OWNER: 'Owner',
};

export const GENERATED_COL_NAMES = [
  COL_NAMES.COMPOUND_NAME,
  COL_NAMES.COMPOUND_COMMENTS,
  COL_NAMES.COMPOUND_MOL_WEIGHT,
  COL_NAMES.SALT_MASS,
  COL_NAMES.BATCH_MOL_WEIGHT,
];
