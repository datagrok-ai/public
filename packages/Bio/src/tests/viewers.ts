import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';

import {category, test, testViewer} from '@datagrok-libraries/test/src/test';
import {readDataframe} from './utils';


// category('viewers', () => {
//   const viewers = DG.Func.find({package: 'Bio', tags: ['viewer']}).map((f) => f.friendlyName);
//   for (const v of viewers) {
//     test(v, async () => {
//       const df = await readDataframe('samples/FASTA_DNA.csv');
//       await df.meta.detectSemanticTypes();
//       await grok.data.detectSemanticTypes(df);
//       await testViewer(v, df, {detectSemanticTypes: true});
//     });
//   }
// });
