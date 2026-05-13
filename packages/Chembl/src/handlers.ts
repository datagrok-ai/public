import * as grok from 'datagrok-api/grok';
import {DBExplorer} from '@datagrok-libraries/db-explorer/src/db-explorer';
import {moleculeRenderer} from '@datagrok-libraries/db-explorer/src/renderer';
import {explorerConfig} from './explorer-config';


export function registerChemblIdHandler() {
  const exp = DBExplorer.initFromConfig(explorerConfig);
  if (!exp) {
    grok.shell.error('Failed to load db-explorer config');
    return;
  }
  exp.addCustomRenderer((_, colName, value) => {
    const lc = colName?.toLowerCase() || '';
    return (lc === 'structure' || lc.includes('smiles')) && typeof value === 'string' && grok.chem.checkSmiles(value);
  }, (value) => moleculeRenderer(value as string));

  console.log('CHEMBL object handlers registered');
}
