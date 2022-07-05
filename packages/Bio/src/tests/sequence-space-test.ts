import {before, category, test, expect} from '@datagrok-libraries/utils/src/test';
import * as DG from "datagrok-api/dg";
import { sequenceSpace } from '../utils/sequence-space';
import { readDataframe } from './utils';
//import * as grok from 'datagrok-api/grok';

category('sequenceSpace', async () => {

    let testFastaDf: DG.DataFrame;
  
    before(async () => {
        //@ts-ignore
        testFastaDf = await readDataframe('sample_FASTA.csv');
    });
  
  
    test('sequenceSpaceOpens', async () => {
         //@ts-ignore
        const res = await sequenceSpace(testFastaDf.col('Sequence')!, 't-SNE', 'Levenshtein', ['Embed_X', 'Embed_Y']);
        expect(res.coordinates != undefined, true); 
        expect(res.distance != undefined, true); 
    });
  
  });