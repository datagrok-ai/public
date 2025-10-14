// import * as ui from 'datagrok-api/ui';
// import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, expectArray} from '@datagrok-libraries/utils/src/test';

category('DataFrame: Methods: Join', () => {
  test('join.inPlace.fullTargetColumnList', async () => {
    const tgtDf1: DG.DataFrame = DG.DataFrame.fromCsv(`id,length
id1,48
id2,47
id2,50`);
    const tgtDf1ColList = tgtDf1.columns.names();
    const srcDf2: DG.DataFrame = DG.DataFrame.fromCsv(`id,start,end,seq
id2,1,23,ACGT
id3,1,23,TCGA
id1,1,23,CGTA
`);
    // In case of arg inPlace = true, the value of valueColumns1 arg can't be null,
    // to keep all columns in target use [] (the empty col list) or full list of target data frame columns.
    tgtDf1.join(srcDf2, ['id'], ['id'], tgtDf1ColList, ['seq'], DG.JOIN_TYPE.LEFT, true);
    expectArray(tgtDf1.columns.names(), [...tgtDf1ColList, 'seq']);
  }, {});

  test('join.inPlace.emptyTargetColumnList', async () => {
    const tgtDf1: DG.DataFrame = DG.DataFrame.fromCsv(`id,length
id1,48
id2,47
id2,50`);
    const tgtDf1ColList = tgtDf1.columns.names();
    const srcDf2: DG.DataFrame = DG.DataFrame.fromCsv(`id,start,end,seq
id2,1,23,ACGT
id3,1,23,TCGA
id1,1,23,CGTA
`);
    // In case of arg inPlace = true, the value of valueColumns1 arg can't be null,
    // to keep all columns in target use [] (the empty col list) or full list of target data frame columns.
    tgtDf1.join(srcDf2, ['id'], ['id'], [], ['seq'], DG.JOIN_TYPE.LEFT, true);
    expectArray(tgtDf1.columns.names(), [...tgtDf1ColList, 'seq']);
  }, {});
}, { owner: 'aparamonov@datagrok.ai' });
