import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import {complement} from '../package';

category('complement function', () => {
  test('characters replacement', async () => {
    const expected = 'TACG';
    const input = 'ATGC';

    const actual = complement(input);
    console.log({actual, expected});

    const errMsg = `extected to have "${expected}" by passing "${input}", but got "${actual}"`;
    expect(actual, expected, errMsg);
  });
});
