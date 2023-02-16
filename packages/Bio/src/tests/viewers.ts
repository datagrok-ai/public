import * as DG from 'datagrok-api/dg';
// import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';

import {category, test, testViewer} from '@datagrok-libraries/utils/src/test';
import {readDataframe} from './utils';


category('viewers', () => {
  const viewers = DG.Func.find({package: 'Bio', tags: ['viewer']}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      await testViewer(v, await readDataframe('data/sample_FASTA_DNA.csv'), true);
    });
  }
});
