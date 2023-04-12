// todo: trim duplicates if any

export const enum PREFIXES {
  AS = 'AS',
  SS = 'SS',
  AS1 = 'AS1',
  AS2 = 'AS2'
}

export const enum SEQ_TYPE {
  AS = 'AS',
  SS = 'SS',
  DUPLEX = 'Duplex',
  DIMER = 'Dimer',
}

/** Computable categories of sequence types */
export const enum SEQ_TYPE_CATEGORY {
  AS_OR_SS,
  DUPLEX,
  DIMER,
}

/** Map between types and their categories inferrable from 'Sequence' column */
export const seqTypeToCategoryDict = {
  [SEQ_TYPE.AS]: SEQ_TYPE_CATEGORY.AS_OR_SS,
  [SEQ_TYPE.SS]: SEQ_TYPE_CATEGORY.AS_OR_SS,
  [SEQ_TYPE.DIMER]: SEQ_TYPE_CATEGORY.DIMER,
  [SEQ_TYPE.DUPLEX]: SEQ_TYPE_CATEGORY.DUPLEX,
};

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
