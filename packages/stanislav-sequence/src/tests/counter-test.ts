
import { category, test, expect, expectArray } from '@datagrok-libraries/utils/src/test';

import {
	dna_nucleotide
} from '../package';

category('stanislav-sequence-counter', () => {
	test('count-subsequence', async () => {

		const value = "GATTACA";
		const resultValue = dna_nucleotide(value);

		expect(resultValue, value);
	})
});
