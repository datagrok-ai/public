import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';

import {category, test, testViewer} from '@datagrok-libraries/utils/src/test';
import {aligned1} from './test-data';


category('Viewers', () => {
  const df = DG.DataFrame.fromCsv(aligned1);
  const viewers = DG.Func.find({package: 'Peptides', tags: ['viewer']}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      await testViewer(v, df.clone(), true);
    }, {skipReason: 'GROK-11534'});
  }
});