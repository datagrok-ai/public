import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {parseMaxQuantText} from '../parsers/maxquant-parser';
import {SEMTYPE} from '../utils/proteomics-types';

// Helper: build a minimal proteinGroups TSV for testing
function makeTsv(rows: string[][], headers: string[]): string {
  return [headers.join('\t'), ...rows.map((r) => r.join('\t'))].join('\n');
}

const HEADERS = [
  'Protein IDs', 'Majority protein IDs', 'Gene names',
  'Potential contaminant', 'Reverse', 'Only identified by site',
  'LFQ intensity Sample1', 'LFQ intensity Sample2', 'iBAQ',
];

category('Parsers', () => {
  test('MaxQuant: filters contaminant rows', async () => {
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '1000', '2000', '500'],
      ['P99999', 'P99999', 'TP53', '+', '', '', '3000', '4000', '600'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    expect(df.rowCount, 1);
  });

  test('MaxQuant: filters reverse rows', async () => {
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '1000', '2000', '500'],
      ['P99999', 'P99999', 'TP53', '', '+', '', '3000', '4000', '600'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    expect(df.rowCount, 1);
  });

  test('MaxQuant: filters only-by-site rows', async () => {
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '1000', '2000', '500'],
      ['P99999', 'P99999', 'TP53', '', '', '+', '3000', '4000', '600'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    expect(df.rowCount, 1);
  });

  test('MaxQuant: filters CON__ prefix rows', async () => {
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '1000', '2000', '500'],
      ['CON__P99999', 'CON__P99999', 'KRT1', '', '', '', '3000', '4000', '600'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    expect(df.rowCount, 1);
  });

  test('MaxQuant: filters REV__ prefix rows', async () => {
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '1000', '2000', '500'],
      ['REV__P99999', 'REV__P99999', 'TP53', '', '', '', '3000', '4000', '600'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    expect(df.rowCount, 1);
  });

  test('MaxQuant: detects LFQ intensity columns', async () => {
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '1000', '2000', '500'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    expect(df.col('LFQ intensity Sample1') !== null, true);
    expect(df.col('log2(LFQ intensity Sample1)') !== null, true);
    expect(df.col('log2(LFQ intensity Sample2)') !== null, true);
  });

  test('MaxQuant: log2 transforms intensity values correctly', async () => {
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '1024', '2048', '500'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    const log2Col = df.col('log2(LFQ intensity Sample1)');
    expect(log2Col !== null, true);
    if (log2Col)
      expect(Math.abs(log2Col.get(0)! - 10.0) < 0.001, true);
  });

  test('MaxQuant: zero intensity produces FLOAT_NULL in log2', async () => {
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '0', '2048', '500'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    const log2Col = df.col('log2(LFQ intensity Sample1)');
    expect(log2Col !== null, true);
    if (log2Col)
      expect(log2Col.isNone(0), true);
  });

  test('MaxQuant: assigns semantic types', async () => {
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '1000', '2000', '500'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    expect(df.col('Protein IDs')?.semType, SEMTYPE.PROTEIN_ID);
    expect(df.col('Gene names')?.semType, SEMTYPE.GENE_SYMBOL);
    expect(df.col('LFQ intensity Sample1')?.semType, SEMTYPE.INTENSITY);
    expect(df.col('log2(LFQ intensity Sample1)')?.semType, SEMTYPE.INTENSITY);
  });

  test('MaxQuant: parses primary protein ID from semicolon-delimited', async () => {
    const tsv = makeTsv([
      ['P12345;Q67890;R11111', 'P12345', 'BRCA1;BRCA2', '', '', '', '1000', '2000', '500'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    const primaryId = df.col('Primary Protein ID');
    expect(primaryId !== null, true);
    if (primaryId)
      expect(primaryId.get(0), 'P12345');
    const primaryGene = df.col('Primary Gene Name');
    expect(primaryGene !== null, true);
    if (primaryGene)
      expect(primaryGene.get(0), 'BRCA1');
  });

  test('MaxQuant: no contaminant IDs remain after filtering', async () => {
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '1000', '2000', '500'],
      ['P22222', 'P22222', 'TP53', '', '', '', '3000', '4000', '600'],
      ['CON__P99999', 'CON__P99999', 'KRT1', '+', '', '', '5000', '6000', '700'],
      ['REV__P88888', 'REV__P88888', 'FAKE1', '', '+', '', '7000', '8000', '800'],
      ['P77777', 'P77777', 'FAKE2', '', '', '+', '9000', '1000', '900'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    expect(df.rowCount, 2);
    const idCol = df.col('Protein IDs')!;
    for (let i = 0; i < df.rowCount; i++) {
      const val = idCol.get(i) as string;
      if (val.startsWith('CON__') || val.startsWith('REV__'))
        throw new Error(`Found filtered ID in output: ${val}`);
    }
  });

  test('MaxQuant: assigns intensity semantic type to log2 columns', async () => {
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '1000', '2000', '500'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    expect(df.col('log2(LFQ intensity Sample1)')?.semType, SEMTYPE.INTENSITY);
    expect(df.col('log2(LFQ intensity Sample2)')?.semType, SEMTYPE.INTENSITY);
    expect(df.col('log2(iBAQ)')?.semType, SEMTYPE.INTENSITY);
  });

  test('MaxQuant: handles MQ 2.x dot column names', async () => {
    const headers2x = [
      'Protein IDs', 'Majority protein IDs', 'Gene names',
      'Potential.contaminant', 'Reverse', 'Only.identified.by.site',
      'LFQ intensity Sample1',
    ];
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '1000'],
      ['P99999', 'P99999', 'TP53', '+', '', '', '3000'],
    ], headers2x);
    const df = await parseMaxQuantText(tsv);
    expect(df.rowCount, 1);
  });

  test('MaxQuant: parser tags proteomics.source', async () => {
    // parseMaxQuantText does NOT set df.name — that's the caller's job (importMaxQuant
    // derives it from the filename). The parser DOES tag proteomics.source so
    // downstream pipeline steps can identify the import origin.
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '1000', '2000', '500'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    expect(df.getTag('proteomics.source'), 'maxquant');
  });

  test('MaxQuant: iBAQ columns detected as intensity', async () => {
    const tsv = makeTsv([
      ['P12345', 'P12345', 'BRCA1', '', '', '', '1000', '2000', '500'],
    ], HEADERS);
    const df = await parseMaxQuantText(tsv);
    expect(df.col('iBAQ')?.semType, SEMTYPE.INTENSITY);
    expect(df.col('log2(iBAQ)') !== null, true);
  });
});
