import * as grok from 'datagrok-api/grok';  

import { category, test, expect, expectArray } from '@datagrok-libraries/utils/src/test';
 
category('stanislav-sequence-dna', () => {
	test('dna_nucleotide', async () => {

		const sequence = "GACTACA"; 
		const subsequence = "AC"; 

		const resultValue = await grok.functions.call('StanislavSequence:CountSubsequencePython', { sequence, subsequence });
		expect(resultValue, 2);
	})
});
