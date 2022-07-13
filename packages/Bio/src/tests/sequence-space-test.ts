import {after, before, category, test, expect} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import {sequenceSpace} from '../utils/sequence-space';
import {readDataframe} from './utils';
import * as grok from 'datagrok-api/grok';

category('sequenceSpace', async () => {
  let testFastaDf: DG.DataFrame;

  before(async () => {
    testFastaDf = await readDataframe('samples/sample_FASTA.csv');
    // await grok.data.detectSemanticTypes(testFastaDf);
  });

  after(async () => {
    grok.shell.closeTable(testFastaDf);
  });

  test('sequenceSpaceOpens', async () => {
    const sequenceSpaceParams = {
      seqCol: testFastaDf.col('Sequence')!,
      methodName: 't-SNE',
      similarityMetric: 'Levenshtein',
      embedAxesNames: ['Embed_X', 'Embed_Y']
    };
    const res = await sequenceSpace(sequenceSpaceParams);
    expect(res.coordinates != undefined, true);
    expect(res.distance != undefined, true);
  });
});
