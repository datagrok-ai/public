import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export const HELM_CORE_LIB_FILENAME = '/data/HELMCoreLibrary.json';

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
export const enum RGROUP_FIELDS {
  CAP_GROUP_SMILES = 'capGroupSmiles',
  ALTER_ID = 'alternateId',
  CAP_GROUP_NAME = 'capGroupName',
  LABEL = 'label',
}

// possible values of polymers
export const enum HELM_POLYMER_TYPE {
  PEPTIDE = 'PEPTIDE',
  RNA = 'RNA',
}

// core fields of HELM library object used in toAtomicLevel function
export const HELM_CORE_FIELDS = [
  HELM_FIELDS.SYMBOL,
  HELM_FIELDS.MOLFILE,
  HELM_FIELDS.RGROUPS,
  HELM_FIELDS.NAME,
];

export const SDF_MONOMER_NAME = 'MonomerName';

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

// range of hex nubers used in PepSea library to endode monomers
export const MONOMER_ENCODE_MIN = 0x100;
export const MONOMER_ENCODE_MAX = 0x40A;
