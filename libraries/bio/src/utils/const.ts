import {MonomerType, PolymerType} from '../helm/types';
import {PolymerTypes} from '../helm/consts';
import {RGroup} from '../types/monomer-library';

// hwe exports the HELM type vocabulary as `as const` objects (value-only), not
// TS enums, so we re-create the value+type duality the old `import X = Enum`
// alias provided — consumers use it both as a value (`.RNA`) and as a type.
const HELM_POLYMER_TYPE = PolymerTypes;
type HELM_POLYMER_TYPE = PolymerType;

export {HELM_POLYMER_TYPE};

/** Required HELM library monomer fields:
 * https://github.com/PistoiaHELM/HELMMonomerSets/blob/master/HELMmonomerSchema.json */
export const enum HELM_REQUIRED_FIELD {
  SYMBOL = 'symbol',
  NAME = 'name',
  MOLFILE = 'molfile',
  AUTHOR = 'author',
  ID = 'id',
  RGROUPS = 'rgroups',
  SMILES = 'smiles',
  POLYMER_TYPE = 'polymerType',
  MONOMER_TYPE = 'monomerType',
  CREATE_DATE = 'createDate',
}

// fields of "rgroups" sub-object in HELM library
export const enum HELM_RGROUP_FIELDS {
  CAP_GROUP_SMILES = 'capGroupSmiles',
  // WARNING: both capitalization variants coexist
  CAP_GROUP_SMILES_UPPERCASE = 'capGroupSMILES',
  ALTERNATE_ID = 'alternateId',
  CAP_GROUP_NAME = 'capGroupName',
  LABEL = 'label',
}

export const enum HELM_OPTIONAL_FIELDS {
  NATURAL_ANALOG = 'naturalAnalog',
  META = 'meta', // for SequenceTranslator
}

// todo: remove
export const enum HELM_FIELDS {
  MONOMER_TYPE = 'monomerType',
  SMILES = 'smiles',
  NAME = 'name',
  AUTHOR = 'author',
  MOLFILE = 'molfile',
  NATURAL_ANALOG = 'naturalAnalog',
  RGROUPS = 'rgroups',
  CREATE_DATE = 'createDate',
  ID = 'id',
  POLYMER_TYPE = 'polymerType',
  SYMBOL = 'symbol'
}

export const SDF_MONOMER_NAME = 'MonomerName';

// todo: ideally, keys should be expressed via constants
export const jsonSdfMonomerLibDict = {
  'monomerType': null, // -> Backbone
  'smiles': null,
  'name': 'Name',
  'author': null,
  'molfile': 'molecule',
  'naturalAnalog': 'MonomerNaturalAnalogCode',
  'rgroups': 'MonomerCaps',
  'createDate': null,
  'id': null,
  'polymerType': 'MonomerType',
  'symbol': 'MonomerName'
};

export const DUMMY_MONOMER = {
  'monomerType': 'Backbone' as MonomerType,
  'smiles': '',
  'name': '',
  'author': 'Datagrok',
  'molfile': '',
  'naturalAnalog': '',
  'rgroups': [] as RGroup[],
  'createDate': null,
  'id': 0,
  'polymerType': 'PEPTIDE' as PolymerType,
  'symbol': ''
} as const;

export const RIBOSE_SYMBOL = 'r';
export const DEOXYRIBOSE_SYMBOL = 'd';
export const PHOSPHATE_SYMBOL = 'p';
export const HELM_WRAPPERS_REGEXP = new RegExp(
  `[${RIBOSE_SYMBOL}${DEOXYRIBOSE_SYMBOL}]\\((\\w)\\)${PHOSPHATE_SYMBOL}?`,
  'g'
);

