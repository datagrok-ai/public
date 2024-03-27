import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';

import {category, test, testViewer} from '@datagrok-libraries/utils/src/test';
import {readDataframe} from './utils';


category('viewers', () => {
  const viewers = DG.Func.find({package: 'Bio', tags: ['viewer']}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      const df = await readDataframe('samples/FASTA_DNA.csv');
      await testViewer(v, df, {detectSemanticTypes: true});
    }, {
      skipReason: {
        'Sequence Similarity Search': 'GROK-13162',
        'Sequence Diversity Search': 'GROK-13162',
        'WebLogo': 'GROK-13162',
        'VdRegions': 'GROK-13162',
      }[v],
    });
  }
});
