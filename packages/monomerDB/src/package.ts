/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as MonomerLib from '@datagrok-libraries/bio/src/types/monomer-library';
import {DBLibraryProvider} from './libraryProvider';
import {DBExplorerConfig} from '@datagrok-libraries/db-explorer/src/types';
import {DBExplorer} from '@datagrok-libraries/db-explorer/src/db-explorer';
import {moleculeRenderer} from '@datagrok-libraries/db-explorer/src/renderer';
export * from './package.g';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

export class PackageFunctions {
  @grok.decorators.func({
    name: 'Monomers DB Provider',
    outputs: [{type: 'object', name: 'result'}],
    meta: {role: 'monomerLibProvider'},
  })
  static async getMonomerDBProvider(): Promise<MonomerLib.IMonomerLibProvider> {
    return new DBLibraryProvider();
  }

  @grok.decorators.init()
  static async init(): Promise<void> {
    const config: DBExplorerConfig = {
      'connectionName': 'Mdb1',
      'schemaName': 'mdb1',
      'dataSourceName': 'postgresDart',
      'entryPoints': {
        'DG_MONOMERDB_ID': {
          'table': 'monomers',
          'column': 'grok_identifier',
          'regexpExample': {
            'example': 'GROKMONO-1',
            'nonVariablePart': 'GROKMONO-',
            'regexpMarkup': 'GROKMONO-{d}'
          }
        },
      },
      'joinOptions': [
        {
          'fromTable': 'monomers',
          'columnName': 'library_id',
          'tableName': 'monomer_libraries',
          'onColumn': 'id',
          'select': ['friendly_name as library_name']
        },
      ],
      'headerNames': {
        'smiles': 'Compound'
      },
      'customSelectedColumns': {
        'monomers': [
          'name',
          'symbol',
          'library_name',
          'smiles',
        ]
      }
    };

    const exp = DBExplorer.initFromConfig(config);
    if (!exp) {
      grok.shell.error('Failed to initialize Monomer DB Explorer');
      return;
    }
    exp.addCustomRenderer((_, colName, value) => {
      const lc = colName?.toLowerCase() || '';
      return (lc === 'structure' || lc.includes('smiles') || lc.includes('compound_structure')) && typeof value === 'string' && grok.chem.checkSmiles(value);
    }, (value) => moleculeRenderer(value as string));

    // exp.addDefaultHeaderReplacerColumns(['units']);
    console.log('MonomerDB object handlers registered');
  }
}
