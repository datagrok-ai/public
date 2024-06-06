import * as grok from 'datagrok-api/grok';  

import { category, test, expect, expectArray } from '@datagrok-libraries/utils/src/test';
import { complement, dna_nucleotide } from '../package';
 
category('stanislav-sequence-dna', () => {
	test('dna_nucleotide', async () => { 
		const value = "GATTACA";
		const resultValue = dna_nucleotide(value);

		expect(resultValue, value);

	})

	test('dna_complement', async () => { 
		const value = "GATTACA";
		const resultValue = complement(value);
		
		expect(resultValue, "CTAATGT");
	})
});
