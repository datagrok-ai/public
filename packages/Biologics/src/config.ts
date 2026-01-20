import {DBExplorerConfig} from '@datagrok-libraries/db-explorer/src/types';

export const biologicsConfig: DBExplorerConfig = {
  'connectionName': 'Biologics',
  'schemaName': 'biologics',
  'nqName': 'Biologics:biologics',
  'dataSourceName': 'postgresDart',
  'entryPoints': {
    'DG_BIOLOGICS_DRUG_ID': {
      'table': 'drugs',
      'column': 'identifier',
      'matchRegexp': 'GROKMOL-\\d{6}',
      'regexpExample': {
        'example': 'GROKMOL-000001',
        'nonVariablePart': 'GROKMOL-',
        'regexpMarkup': 'GROKMOL-{d6}'
      }
    },
    'DG_BIOLOGICS_PEPTIDE_ID': {
      'table': 'peptides',
      'column': 'identifier',
      'matchRegexp': 'GROKPEP-\\d{6}',
      'regexpExample': {
        'example': 'GROKPEP-000001',
        'nonVariablePart': 'GROKPEP-',
        'regexpMarkup': 'GROKPEP-{d6}'
      }
    },
    'DG_BIOLOGICS_SEQUENCE_ID': {
      'table': 'sequences',
      'column': 'identifier',
      'matchRegexp': 'GROKSEQ-\\d{6}',
      'regexpExample': {
        'example': 'GROKSEQ-000001',
        'nonVariablePart': 'GROKSEQ-',
        'regexpMarkup': 'GROKSEQ-{d6}'
      }
    },
    'DG_BIOLOGICS_LINKER_ID': {
      'table': 'linkers',
      'column': 'identifier',
      'matchRegexp': 'GROKLINKER-\\d{6}',
      'regexpExample': {
        'example': 'GROKLINKER-000001',
        'nonVariablePart': 'GROKLINKER-',
        'regexpMarkup': 'GROKLINKER-{d6}'
      }
    },
    'DG_BIOLOGICS_ADC_ID': {
      'table': 'adc',
      'column': 'identifier',
      'matchRegexp': 'GROKADC-\\d{6}',
      'regexpExample': {
        'example': 'GROKADC-000001',
        'nonVariablePart': 'GROKADC-',
        'regexpMarkup': 'GROKADC-{d6}'
      }
    },
    'DG_BIOLOGICS_ORGANISM_ID': {
      'table': 'target_organisms',
      'column': 'identifier',
      'matchRegexp': 'GROKORG-\\d{6}',
      'regexpExample': {
        'example': 'GROKORG-000001',
        'nonVariablePart': 'GROKORG-',
        'regexpMarkup': 'GROKORG-{d6}'
      }
    },
    'DG_BIOLOGICS_PURIFICATION_ID': {
      'table': 'purification_batches',
      'column': 'identifier',
      'matchRegexp': 'GROKPUR-\\d{6}',
      'regexpExample': {
        'example': 'GROKPUR-000001',
        'nonVariablePart': 'GROKPUR-',
        'regexpMarkup': 'GROKPUR-{d6}'
      }
    },
    'DG_BIOLOGICS_EXPRESSION_ID': {
      'table': 'expression_batches',
      'column': 'identifier',
      'matchRegexp': 'GROKEXP-\\d{6}',
      'regexpExample': {
        'example': 'GROKEXP-000001',
        'nonVariablePart': 'GROKEXP-',
        'regexpMarkup': 'GROKEXP-{d6}'
      }
    }
  },
  'joinOptions': [
    {
      'fromTable': 'adc',
      'columnName': 'drug_id',
      'tableName': 'drugs',
      'onColumn': 'id',
      'select': ['smiles as compound_structure', 'name as compound_name']
    },
    {
      'fromTable': 'adc',
      'columnName': 'antibody_id',
      'tableName': 'sequences',
      'onColumn': 'id',
      'select': ['sequence as antibody_sequence', 'name as sequence_name']
    },
    {
      'fromTable': 'assay_results',
      'columnName': 'assay_id',
      'tableName': 'assay_types',
      'onColumn': 'id',
      'select': ['name']
    }
  ],
  'headerNames': {
    'linkers': 'linker_type',
  },
  'uniqueColumns': {
    'adc': 'identifier',
    'drugs': 'identifier',
    'sequences': 'identifier',
    'linkers': 'identifier',
    'target_organisms': 'identifier',
    'purification_batches': 'identifier',
    'expression_batches': 'identifier',
    'peptides': 'identifier',
  },
  'customSelectedColumns': {
    'adc': [
      'name',
      'identifier',
      'antibody_sequence',
      'antibody_id',
      'drug_id',
      'linker_id',
      'compound_structure',
      'glyph'
    ],
    'peptides': [
      'name',
      'identifier',
      'helm'
    ]
  }
};
