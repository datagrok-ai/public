import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import { category, test, expect, expectArray, expectTable } from '@datagrok-libraries/utils/src/test';

import {
    dna_nucleotide,
    fuzzyJoin
} from '../package';

category('stanislav-sequence-counter', () => {
    test('count-subsequence', async () => {
        const sequence = "GACTACA";
        const subsequence = "AC";

        const resultValue = await grok.functions.call('StanislavSequence:CountSubsequencePython', { sequence, subsequence });
        expect(resultValue, 2);
    })

    test('count-subsequence-dataframe', async () => {
        const testCsv = `Sequence,id
        GATTACA,1997
        AGGCGGA,1984
        TTTAGGC,2021`
        const columnName = "Sequence";
        const subsequence = "GG";
        const df = DG.DataFrame.fromCsv(testCsv);

        const resultValue = await grok.functions.call('StanislavSequence:CountSubsequencePythonDataframe', { df, columnName, subsequence });
        const expectedValues = [0, 2, 1];
        console.log(resultValue);
        // const gainedValues = resultValue.getCol("N("+subsequence+")").getRawData();

        if(resultValue.length !== expectedValues.length){
            expect(1,2);
        }

        for (let i = 0; i < expectedValues.length; i++) {
            expect(expectedValues[i],resultValue[i]);
        }
    })
    test('count-subsequence-dataframes', async () => {
        const testCsv = `Sequence
        GATTACA` 
        const testCsv2 = `Sequence
        GATTACA
        AGGCGGA
        GATTACC` 
        const expectedRes = `Sequence, Subsequence Count
        GATTACA,3
        GATTACA,2
        AGGCGGA,0
        GATTACC,1` 

        const df1 = DG.DataFrame.fromCsv(testCsv);

        const result1 = fuzzyJoin(df1, DG.DataFrame.fromCsv(testCsv2), df1.col('Sequence')!, 6)

        expectTable(result1, DG.DataFrame.fromCsv(expectedRes));
    })
});
