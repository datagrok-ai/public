import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export const HELM_CORE_LIB_FILENAME = '/data/HELMCoreLibrary.json';

/** Required HELM library monomer fields:
 * https://github.com/PistoiaHELM/HELMMonomerSets/blob/master/HELMmonomerSchema.json */
export const enum HELM_REQUIRED_FIELDS {
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

export const enum HELM_OPTIONAL_FIELDS {
  NATURAL_ANALOG = 'naturalAnalog',
  META = 'meta', // for SequenceTranslator
}

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

// fields of "rgroups" sub-object in HELM library
export const enum HELM_RGROUP_FIELDS {
  CAP_GROUP_SMILES = 'capGroupSmiles',
  CAP_GROUP_SMILES_UPPERCASE = 'capGroupSMILES', // alas, both variants coexist
  ALTERNATE_ID = 'alternateId',
  CAP_GROUP_NAME = 'capGroupName',
  LABEL = 'label',
}

// possible values of polymers
export const enum HELM_POLYMER_TYPE {
  PEPTIDE = 'PEPTIDE',
  RNA = 'RNA',
}

export const enum HELM_MONOMER_TYPE {
  BACKBONE = 'Backbone',
  TERMINAL = 'Terminal',
  BRANCH = 'Branch',
}

// core fields of HELM library object used in toAtomicLevel function
export const HELM_CORE_FIELDS = [
  HELM_FIELDS.SYMBOL,
  HELM_FIELDS.MOLFILE,
  HELM_FIELDS.RGROUPS,
  HELM_FIELDS.NAME,
  // HELM_FIELDS.MONOMER_TYPE, // add if terminal monomers for PEPTIDEs to be
  // supported
];

export const SDF_MONOMER_NAME = 'MonomerName';

/** For Enumerator  */
export const helmFieldsToEnumeratorInputFields = {
  [HELM_REQUIRED_FIELDS.SYMBOL]: 'Short Name',
  [HELM_REQUIRED_FIELDS.NAME]: 'Medium Name',
  [HELM_REQUIRED_FIELDS.SMILES]: 'SMILES',
};

/** For Enumerator  */
export const rGroupsDummy = [
  {
    'capGroupSmiles': '[*:1][H]',
    'alternateId': 'R1-H',
    'capGroupName': 'H',
    'label': 'R1'
  },
  {
    'capGroupSmiles': 'O[*:2]',
    'alternateId': 'R2-OH',
    'capGroupName': 'OH',
    'label': 'R2'
  }
];

// todo: ideally, keys should be expressed via constants
export const jsonSdfMonomerLibDict = {
  'monomerType': null,
  'smiles': null,
  'name': 'MonomerName',
  'author': null,
  'molfile': 'molecule',
  'naturalAnalog': 'MonomerNaturalAnalogCode',
  'rgroups': 'MonomerCaps',
  'createDate': null,
  'id': null,
  'polymerType': 'MonomerType',
  'symbol': 'MonomerCode'
};

export const dummyMonomer = {
  'monomerType': 'Backbone',
  'smiles': '',
  'name': '',
  'author': 'Datagrok',
  'molfile': '',
  'naturalAnalog': '',
  'rgroups': [],
  'createDate': null,
  'id': 0,
  'polymerType': 'PEPTIDE',
  'symbol': ''
};
// range of hex nubers used in PepSea library to endode monomers
export const MONOMER_ENCODE_MIN = 0x100;
export const MONOMER_ENCODE_MAX = 0x40A;
