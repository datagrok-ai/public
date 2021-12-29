/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PDBViewerHelper} from './pdb-viewer-helper';

export const _package = new DG.Package();

//name: testPDB
export async function test() {
  const h = new PDBViewerHelper('2v0a');

  await h.fetchInfo();

  console.log(h.body);
  console.log(h.entities);
}
