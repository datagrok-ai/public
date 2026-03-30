import { category, expect, test } from '@datagrok-libraries/test/src/test';
import { countSubsequencePythonPackageTS } from '../package';

category('CountSubsequencePythonPackageTS', () => {
  test('basic count', async () => {
    expect(await countSubsequencePythonPackageTS('ATGATC', 'A'), 2);
  });

  test('overlapping count', async () => {
    expect(await countSubsequencePythonPackageTS('AAAA', 'AA'), 3);
  });

  test('sequence example', async () => {
    expect(await countSubsequencePythonPackageTS('fasta: gttctctacc', 'acc'), 1);
  });
});