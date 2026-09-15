/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import {domains} from '@datagrok-libraries/u2/src/dg/index.js';
export * from './package.g';

// the skins of everything the u2 domain stack renders
import '@datagrok-libraries/u2/src/dg/domain/styles.js';

export const _package = new DG.Package();

//name: Stockroom
//description: Chemical stockroom on the GHS classification — the zero-code app over databases/stockroom/schema.json
//tags: app
//meta.icon: images/flask.svg
//output: view result
export async function stockroomApp(): Promise<DG.ViewBase> {
  return (await domains.table('stockroom.substance')).app();
}
