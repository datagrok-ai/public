import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import { category, expect, test } from '@datagrok-libraries/utils/src/test';

category('detectors', () => {
    test('detect-dna_nucleotide', async () => {
        const testCsv = `sequence,id
        GATTACA,1997
        ATTCGGA,1984
        TTTAGGC,2021`

        const df = DG.DataFrame.fromCsv(testCsv);
        const col1 = df.columns.byName('sequence');
        await grok.data.detectSemanticTypes(df);

        expect(col1.semType, "dna_nucleotide");
    });
});