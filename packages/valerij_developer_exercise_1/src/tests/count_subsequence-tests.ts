import { category, expect, test } from '@datagrok-libraries/test/src/test';
import { countSubsequencePythonPackageTS } from '../package';

category('CountSubsequencePythonPackageTS', () => {
  test('basic count', async () => {
    const result = await countSubsequencePythonPackageTS('ATGATC', 'A');
    expect(result, 2);
  });

  test('overlapping count', async () => {
    const result = await countSubsequencePythonPackageTS('AAAA', 'AA');
    expect(result, 3);
  });

  test('sequence example', async () => {
    const result = await countSubsequencePythonPackageTS('fasta: gttctctacc', 'acc');
    expect(result, 1);
  });
});