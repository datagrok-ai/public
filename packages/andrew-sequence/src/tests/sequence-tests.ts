import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import {complement} from '../package';

category('complement function', (): void => {
  const checkComplementResult = (input: string, expected: string): void => {
    const actual = complement(input);

    const errMsg = `Extected to have "${expected}" by passing "${input}", but got "${actual}".`;
    expect(actual, expected, errMsg);
  };

  test('upper case characters replacement', async (): Promise<void> => {
    const expected = 'TACG';
    const input = 'ATGC';

    checkComplementResult(input, expected);
  });

  test('lower case characters replacement', async (): Promise<void> => {
    const expected = 'tacg';
    const input = 'atgc';

    checkComplementResult(input, expected);
  });

  test('lower & upper case characters replacement including spaces', async (): Promise<void> => {
    const expected = 'tacg AaAt GCGC';
    const input = 'atgc TtTa CGCG';

    checkComplementResult(input, expected);
  });
});
