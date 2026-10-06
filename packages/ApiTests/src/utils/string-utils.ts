import * as DG from 'datagrok-api/dg';
import {category, expect, expectFloat, test} from '@datagrok-libraries/test/src/test';

// npm datagrok-api@1.27.9 typings predate the 1.28 StringUtils distance statics
const strings = DG.StringUtils as any;

category('Utils: StringUtils', () => {
  test('levenshteinDistance', async () => {
    expect(strings.levenshteinDistance('', ''), 0);
    expect(strings.levenshteinDistance('age', 'age'), 0);
    expect(strings.levenshteinDistance('a', 'xyz'), 1);
    expectFloat(strings.levenshteinDistance('kitten', 'sitting'), 3 / 7, 0.0001);
  });

  test('jaroWinklerDistance', async () => {
    expect(strings.jaroWinklerDistance('weight', 'weight'), 0);
    expect(strings.jaroWinklerDistance('abc', 'xyz'), 1);
    expectFloat(strings.jaroWinklerDistance('martha', 'marhta'), 0.0389, 0.0001);
  });
}, {owner: 'lbankurova@datagrok.ai'});
