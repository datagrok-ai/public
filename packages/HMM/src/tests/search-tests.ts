import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';

/** Fixtures from the Rusty-HMMER corpora (synced): inputs and canonical Linux
 * SSE C results, printed values (%.1f scores, %.2g E-values). */
async function sequenceTable(file: string): Promise<{df: DG.DataFrame; col: DG.Column<string>}> {
  const df = DG.DataFrame.fromCsv(await _package.files.readAsText(`tests/${file}`));
  const col = df.col('sequence') as DG.Column<string>;
  col.semType = DG.SEMTYPE.MACROMOLECULE;
  col.meta.units = 'fasta';
  col.setTag('aligned', 'SEQ');
  col.setTag('alphabet', 'PT');
  return {df, col};
}

/** Expected rows as strings (DG float columns are 32-bit: E-values would underflow). */
async function expectedRows(file: string): Promise<Record<string, string>[]> {
  const [header, ...lines] = (await _package.files.readAsText(`tests/${file}`)).trim().split('\n');
  const names = header.split(',');
  return lines.map((line) => Object.fromEntries(line.split(',').map((v, i) => [names[i], v])));
}

const g2 = (x: number) => Number(x.toPrecision(2));
const f1 = (x: number) => Number(x.toFixed(1));

category('HMMER search and domains', () => {
  test('profile search matches C hmmsearch', async () => {
    const {df, col} = await sequenceTable('search-targets.csv');
    const model = await _package.files.readAsText('tests/search-model.hmm');
    const expected = await expectedRows('search-expected.csv');
    const result: DG.DataFrame = await grok.functions.call('HMM:searchWithHmm',
      {table: df, sequence: col, model, evalue: 10, gathering: false});
    const names = df.col('name')!.toList() as string[];
    const found = new Set<string>();
    for (const e of expected) {
      const row = names.indexOf(e.target);
      found.add(e.target);
      expect(f1(result.columns.byIndex(0).get(row)), Number(e.score), `${e.target} score`);
      expect(g2(result.columns.byIndex(1).get(row)), Number(e.evalue), `${e.target} E-value`);
      expect(result.columns.byIndex(2).get(row), Number(e.reported_domains), `${e.target} domains`);
    }
    for (let row = 0; row < df.rowCount; row++)
      if (!found.has(names[row])) expect(result.columns.byIndex(0).isNone(row), true, `${names[row]} not reported`);
  });

  test('domain annotation matches C hmmscan --cut_ga', async () => {
    const {df, col} = await sequenceTable('scan-queries.csv');
    const library = await _package.files.readAsText('tests/scan-library.hmm');
    const expected = await expectedRows('scan-expected.csv');
    const result: DG.DataFrame = await grok.functions.call('HMM:findDomains',
      {table: df, sequence: col, library, evalue: 0.01, gathering: true});
    const names = df.col('name')!.toList() as string[];
    const key = (query: string, model: string, from: number, to: number, score: number, evalue: number) =>
      [query, model, from, to, f1(score), g2(evalue)].join('|');
    const actual = new Set<string>();
    for (let i = 0; i < result.rowCount; i++) {
      actual.add(key(names[result.get('row', i)], result.get('model', i), result.get('from', i), result.get('to', i),
        result.get('score', i), result.get('evalue', i)));
    }
    const wanted = new Set<string>();
    for (const e of expected)
      wanted.add(key(e.query, e.model, Number(e.env_from), Number(e.env_to), Number(e.score), Number(e.i_evalue)));

    expect(actual.size, wanted.size, 'domain count');
    for (const k of wanted) expect(actual.has(k), true, `missing ${k}`);
  });
});
