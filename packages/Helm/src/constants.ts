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

export const SMILES = 'smiles';
export const RGROUPS = 'rgroups';
export const MONOMER_SYMBOL = 'symbol';
export const RGROUP_CAP_GROUP_SMILES = 'capGroupSmiles';
export const RGROUP_ALTER_ID = 'alternateId';
export const RGROUP_CAP_GROUP_NAME = 'capGroupName';
export const RGROUP_LABEL = 'label';
export const SDF_MONOMER_NAME = 'MonomerName';

export const enum TAGS {
  cellRendererRenderError = '.cell-renderer.render.error',
}

/** Names of package properties declared in the `properties` section of `package.json`. */
export const enum HelmPackagePropertiesNames {
  MonomersPerRow = 'MonomersPerRow',
}

/**
 * Fallback for `MonomersPerRow` when the package property is missing or
 * unparseable (a package loaded before its settings, an old install). Kept in
 * step with the `defaultValue` in package.json and with hwe's own default.
 */
export const DEFAULT_MONOMERS_PER_ROW = 20;
