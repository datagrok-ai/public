import {DBExplorerConfig} from '@datagrok-libraries/db-explorer/src/types';

export const explorerConfig: DBExplorerConfig = {
  'connectionName': 'CHEMBL',
  'schemaName': 'public',
  'dataSourceName': 'postgres',
  'nqName': 'Chembl:Chembl',
  'entryPoints': {
    'CHEMBL_ID': {
      'table': 'molecule_dictionary',
      'column': 'chembl_id',
      'regexpExample': {
        'example': 'CHEMBL1234',
        'nonVariablePart': 'CHEMBL',
        'regexpMarkup': 'CHEMBL[0-9]+'
      },
      'matchRegexp': 'CHEMBL\\d+'
    },
    'molregno': {
      'table': 'molecule_dictionary',
      'column': 'molregno'
    }
  },
  'joinOptions': [
    {
      'fromTable': 'molecule_dictionary',
      'columnName': 'molregno',
      'tableName': 'compound_structures',
      'onColumn': 'molregno',
      'select': ['canonical_smiles', 'standard_inchi']
    },
    // example of join with fromSchema and onSchema
    // {
    //   'fromSchema': 'public',
    //   'fromTable': 'molecule_dictionary',
    //   'columnName': 'molregno',
    //   'onSchema': 'rdk',
    //   'tableName': 'mols',
    //   'onColumn': 'molregno',
    //   'select': ['m']
    // },
    {
      'fromTable': 'compound_structural_alerts',
      'columnName': 'alert_id',
      'tableName': 'structural_alerts',
      'onColumn': 'alert_id',
      'select': ['alert_name']
    }
  ],
  'headerNames': {
    'action_type': 'action_type',
    'activities': 'type',
    'activity_properties': 'standard_type',
    'assays': 'tid',
    'chembl_id_lookup': 'chembl_id',
    'drug_indication': 'efo_term',
    'drug_mechanism': 'action_type',
    'drug_warning': 'warning_type',
    'frac_classification': 'active_ingredient',
    'molecule_synonyms': 'synonyms',
    'compound_records': 'record_id',
    'compound_structural_alerts': 'alert_name',
  },
  explicitReferences: [
    {
      schema: 'rdk',
      table: 'mols',
      column: 'molregno',
      refSchema: 'public',
      refTable: 'molecule_dictionary',
      refColumn: 'molregno'
    },
    {
      schema: 'rdk',
      table: 'fps',
      column: 'molregno',
      refSchema: 'public',
      refTable: 'molecule_dictionary',
      refColumn: 'molregno'
    }
  ],

  'uniqueColumns': {
    'activities': 'activity_id',
    'compound_records': 'record_id',
    'molecule_synonyms': 'molsyn_id'
  },
  'customSelectedColumns': {
    'molecule_dictionary': [
      'canonical_smiles',
      'standard_inchi',
      'molregno',
      'pref_name',
      'chembl_id',
      'structure_type',
      'molecule_type',
      'first_approval',
      'indication_class',
      'chirality',
    ]
  }
};
