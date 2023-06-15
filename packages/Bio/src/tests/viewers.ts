import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';

import {category, delay, test} from '@datagrok-libraries/utils/src/test';
import {readDataframe} from './utils';


category('viewers', () => {
  const viewers = DG.Func.find({package: 'Bio', tags: ['viewer']}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      const df = await readDataframe('data/sample_FASTA_DNA.csv');
      const tv = grok.shell.addTableView(df);
      await grok.data.detectSemanticTypes(df);
      tv.addViewer(v);
      await delay(2000);
      // await testViewer(v, df, {detectSemanticTypes: true});
    });
  }
});
