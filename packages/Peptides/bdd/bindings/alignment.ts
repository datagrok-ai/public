import {expect, type Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';

declare const grok: any;

export const alignedPositions = Then('the split monomer columns of row {int} should match its peptide sequence',
  async (page: Page, row: number) => {
    const result = await page.evaluate((number) => {
      const table = grok.shell.t;
      if (number < 1 || number > table.rowCount)
        throw new Error(`row ${number} is outside the peptide table`);
      const sequence = table.getCol('AlignedSequence');
      const separator = sequence.getTag('separator');
      if (sequence.meta.units !== 'separator' || !separator)
        throw new Error('the alignment oracle requires separator-notation peptides');
      const tokens = String(sequence.get(number - 1)).split(separator);
      const positions = table.columns.toList().filter((c: any) => c.getTag('isPositionCol') === 'true');
      return {
        expected: tokens,
        actual: positions.sort((a: any, b: any) => Number(a.name) - Number(b.name))
          .map((column: any) => String(column.get(number - 1) ?? '')),
      };
    }, row);
    expect(result.actual, `all split monomers of row ${row}`).toEqual(result.expected);
  }, {description: 'all tagged position columns in numerical order checked against the stored separator sequence'});
