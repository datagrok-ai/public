import {expect, type Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';

declare const grok: any;

export const invariantCounts = Then('the invariant-map export should match the monomer counts of table {string}',
  async (page: Page, source: string) => {
    const result = await page.evaluate((name) => {
      const table = grok.shell.tables.find((t: any) => t.name === name);
      if (!table || table.rowCount === 0)
        throw new Error(`no populated source table "${name}"`);
      const sequences = table.getCol('AlignedSequence');
      const separator = sequences.getTag('separator');
      if (sequences.meta.units !== 'separator' || !separator)
        throw new Error('the count oracle requires separator-notation peptides');
      const split = sequences.toList().map((s: string) => s.split(separator)) as string[][];
      const width = Math.max(...split.map((s) => s.length));
      const monomers = [...new Set(split.flat().filter(Boolean))].sort();
      const exported = grok.shell.t;
      const columns = ['AAR', ...Array.from({length: width}, (_, i) => String(i + 1))];
      const mismatches: string[] = [];
      if (exported.columns.names().join(',') !== columns.join(','))
        mismatches.push(`columns: ${exported.columns.names().join(',')}`);
      if (exported.getCol('AAR').toList().join(',') !== monomers.join(','))
        mismatches.push(`monomers: ${exported.getCol('AAR').toList().join(',')}`);
      for (let position = 0; position < width; position++) {
        const column = exported.getCol(String(position + 1));
        if (column.type !== 'int')
          mismatches.push(`${column.name} has type ${column.type}`);
        for (let row = 0; row < monomers.length; row++) {
          const expected = split.filter((s) => s[position] === monomers[row]).length;
          if (column.get(row) !== expected)
            mismatches.push(`${monomers[row]} at ${position + 1}: ${column.get(row)}, expected ${expected}`);
        }
      }
      return mismatches;
    }, source);
    expect(result, 'every exported monomer-position count, including zeros').toEqual([]);
  }, {description: 'all cells and integer column types checked against the original separator sequences, excluding gaps'});

export const mutationPairs = Then('the mutation-cliff export should contain every single-mutation pair from table {string}',
  async (page: Page, source: string) => {
    const result = await page.evaluate((name) => {
      const table = grok.shell.tables.find((t: any) => t.name === name);
      if (!table || table.rowCount === 0)
        throw new Error(`no populated source table "${name}"`);
      const sequenceCol = table.getCol('AlignedSequence');
      const separator = sequenceCol.getTag('separator');
      if (sequenceCol.meta.units !== 'separator' || !separator)
        throw new Error('the mutation oracle requires separator-notation peptides');
      const sequences = sequenceCol.toList() as string[];
      const activity = table.getCol('IC50').toList() as number[];
      const byActivity = new Map(activity.map((value, row) => [value, row]));
      if (byActivity.size !== table.rowCount)
        throw new Error('the fixture must have unique IC50 values to identify exported source rows');
      const split = sequences.map((s) => s.split(separator));
      const expected = new Set<string>();
      for (let a = 0; a < split.length; a++) {
        for (let b = a + 1; b < split.length; b++) {
          if (split[a].length !== split[b].length)
            throw new Error('the mutation oracle requires equally aligned sequences');
          if (split[a].filter((monomer, position) => monomer !== split[b][position]).length === 1)
            expected.add(`${a}:${b}`);
        }
      }
      if (!expected.size)
        throw new Error('the source fixture contains no single-mutation pairs');
      const exported = grok.shell.t;
      const actual = new Set<string>();
      const mismatches: string[] = [];
      for (let row = 0; row < exported.rowCount; row++) {
        const activity1 = exported.get('Seq 1 IC50', row);
        const activity2 = exported.get('Seq 2 IC50', row);
        const a = byActivity.get(activity1);
        const b = byActivity.get(activity2);
        if (a === undefined || b === undefined) {
          mismatches.push(`row ${row + 1}: activity absent from the source`);
          continue;
        }
        const key = `${Math.min(a, b)}:${Math.max(a, b)}`;
        if (actual.has(key) || !expected.has(key))
          mismatches.push(`row ${row + 1}: duplicate or non-cliff pair ${key}`);
        actual.add(key);
        if (exported.get('Seq 1', row) !== sequences[a] || exported.get('Seq 2', row) !== sequences[b])
          mismatches.push(`row ${row + 1}: sequence does not match its source activity`);
        const delta = exported.get('Delta', row);
        if (!Number.isFinite(delta) || Math.abs(delta - (activity1 - activity2)) > 1e-10)
          mismatches.push(`row ${row + 1}: wrong activity difference ${delta}`);
      }
      return {missing: [...expected].filter((key) => !actual.has(key)).length, mismatches};
    }, source);
    expect(result, 'unique single-mutation source pairs, sequence identity and activity differences')
      .toEqual({missing: 0, mismatches: []});
  }, {description: 'all Hamming-distance-one row pairs of the unfiltered separator fixture; IC50 identifies source rows'});

export const mutationExtraValues = Then('the mutation-cliff export should preserve {string} values from table {string}',
  async (page: Page, extra: string, source: string) => {
    const mismatches = await page.evaluate(([column, name]) => {
      const table = grok.shell.tables.find((t: any) => t.name === name);
      if (!table || table.rowCount === 0)
        throw new Error(`no populated source table "${name}"`);
      const activity = table.getCol('IC50').toList() as number[];
      const values = table.getCol(column);
      const byActivity = new Map(activity.map((value, row) => [value, values.get(row)]));
      if (byActivity.size !== table.rowCount)
        throw new Error('the fixture must have unique IC50 values to identify exported source rows');
      const exported = grok.shell.t;
      if (!exported.rowCount)
        throw new Error('the mutation-cliff export is empty');
      const mismatches: string[] = [];
      for (const side of [1, 2]) {
        const extras = exported.getCol(`Seq ${side} ${column}`);
        for (let row = 0; row < exported.rowCount; row++) {
          const value = exported.get(`Seq ${side} IC50`, row);
          if (!byActivity.has(value) || extras.get(row) !== byActivity.get(value))
            mismatches.push(`row ${row + 1}, Seq ${side}: ${extras.get(row)}`);
        }
      }
      return mismatches;
    }, [extra, source]);
    expect(mismatches, `exported ${extra} for each member of every mutation pair`).toEqual([]);
  }, {description: 'both extra-column copies checked against the corresponding source rows identified by IC50'});
