import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import $ from 'cash-dom';

import {awaitCheck, category, test, testViewer} from '@datagrok-libraries/utils/src/test';

import {TreeHelper} from '../utils/tree-helper';

import {_package} from '../package-test';


// category('viewers', () => {
//   const viewers = DG.Func.find({package: 'Dendrogram', meta: {role: 'viewer'}}).map((f) => f.friendlyName);
//   for (const v of viewers) {
//     test(v, async () => {
//       const df = await (async () => {
//         const newickStr: string = await _package.files.readAsText('data/tree95.nwk');
//         const treeHelper = new TreeHelper();
//         return treeHelper.newickToDf(newickStr, 'tree95');
//       })();
//       await testViewer(v, df, {detectSemanticTypes: true});
//     });
//   }
// });
