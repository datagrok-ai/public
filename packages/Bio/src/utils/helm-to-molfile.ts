/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export async function helmToMolfile(helmCol: DG.Column<string>): Promise<DG.Column> {
  const polymerMolfile: DG.Column<string> = await grok.functions.call('HELM:getMolfiles', {col: helmCol});
  return polymerMolfile;
}
