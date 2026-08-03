import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {parseFragPipeText} from '../parsers/fragpipe-parser';
import {SEMTYPE} from '../utils/proteomics-types';

/** Builds a minimal FragPipe combined_protein.tsv string. */
function makeTsv(rows: string[][], headers: string[]): string {
  return [headers.join('\t'), ...rows.map((r) => r.join('\t'))].join('\n');
}

const BASE_HEADERS = [
  'Protein', 'Protein ID', 'Entry Name', 'Gene', 'Length', 'Organism', 'Protein Description',
  'Sample1 Spectral Count', 'Sample1 Intensity', 'Sample1 MaxLFQ Intensity',
  'Sample2 Spectral Count', 'Sample2 Intensity', 'Sample2 MaxLFQ Intensity',
];

function baseRow(proteinId: string, gene: string,
  s1Int: string, s1Lfq: string, s2Int: string, s2Lfq: string): string[] {
  return [
    `sp|${proteinId}|${gene}_HUMAN`, proteinId, `${gene}_HUMAN`, gene, '500', 'Homo sapiens', `${gene} protein`,
    '12', s1Int, s1Lfq, '15', s2Int, s2Lfq,
  ];
}

category('FragPipe', () => {
  test('parses standard rows', async () => {
    const tsv = makeTsv([
      baseRow('P12345', 'BRCA1', '10000', '12000', '20000', '24000'),
      baseRow('Q67890', 'TP53', '30000', '36000', '40000', '48000'),
    ], BASE_HEADERS);
    const df = await parseFragPipeText(tsv);
    expect(df.rowCount, 2);
  });

  test('filters contam_ rows', async () => {
    const tsv = makeTsv([
      baseRow('P12345', 'BRCA1', '10000', '12000', '20000', '24000'),
      baseRow('contam_P99999', 'KRT1', '30000', '36000', '40000', '48000'),
    ], BASE_HEADERS);
    const df = await parseFragPipeText(tsv);
    expect(df.rowCount, 1);
    expect(df.col('Protein ID')!.get(0), 'P12345');
  });

  test('filters rev_ rows', async () => {
    const tsv = makeTsv([
      baseRow('P12345', 'BRCA1', '10000', '12000', '20000', '24000'),
      baseRow('rev_P99999', 'DECOY', '30000', '36000', '40000', '48000'),
    ], BASE_HEADERS);
    const df = await parseFragPipeText(tsv);
    expect(df.rowCount, 1);
  });

  test('filters CON__/REV__ rows for MaxQuant-style FASTAs', async () => {
    const tsv = makeTsv([
      baseRow('P12345', 'BRCA1', '10000', '12000', '20000', '24000'),
      baseRow('CON__P99999', 'KRT1', '30000', '36000', '40000', '48000'),
      baseRow('REV__P88888', 'FAKE1', '50000', '60000', '70000', '80000'),
    ], BASE_HEADERS);
    const df = await parseFragPipeText(tsv);
    expect(df.rowCount, 1);
  });

  test('assigns semantic types', async () => {
    const tsv = makeTsv([
      baseRow('P12345', 'BRCA1', '10000', '12000', '20000', '24000'),
    ], BASE_HEADERS);
    const df = await parseFragPipeText(tsv);
    expect(df.col('Protein ID')?.semType, SEMTYPE.PROTEIN_ID);
    expect(df.col('Gene')?.semType, SEMTYPE.GENE_SYMBOL);
    // Only the chosen primary quant (MaxLFQ) gets SEMTYPE.INTENSITY. Redundant
    // bare-Intensity / Razor columns are intentionally NOT tagged so they don't
    // show up in downstream intensity-column pickers.
    expect(df.col('Sample1 MaxLFQ Intensity')?.semType, SEMTYPE.INTENSITY);
  });

  test('prefers MaxLFQ for log2 transform', async () => {
    const tsv = makeTsv([
      baseRow('P12345', 'BRCA1', '1024', '2048', '4096', '8192'),
    ], BASE_HEADERS);
    const df = await parseFragPipeText(tsv);
    // MaxLFQ columns get log2-transformed; bare Intensity ones do not
    expect(df.col('log2(Sample1 MaxLFQ Intensity)') !== null, true);
    expect(df.col('log2(Sample2 MaxLFQ Intensity)') !== null, true);
    expect(df.col('log2(Sample1 Intensity)'), null);
    const log2Col = df.col('log2(Sample1 MaxLFQ Intensity)')!;
    expect(Math.abs(log2Col.get(0)! - 11.0) < 0.001, true);
  });

  test('falls back to bare Intensity when no MaxLFQ present', async () => {
    const headers = [
      'Protein', 'Protein ID', 'Gene',
      'Sample1 Spectral Count', 'Sample1 Intensity',
      'Sample2 Spectral Count', 'Sample2 Intensity',
    ];
    const tsv = makeTsv([
      ['sp|P12345|BRCA1_HUMAN', 'P12345', 'BRCA1', '12', '1024', '15', '4096'],
    ], headers);
    const df = await parseFragPipeText(tsv);
    expect(df.col('log2(Sample1 Intensity)') !== null, true);
    expect(df.col('log2(Sample2 Intensity)') !== null, true);
    const log2Col = df.col('log2(Sample1 Intensity)')!;
    expect(Math.abs(log2Col.get(0)! - 10.0) < 0.001, true);
  });

  test('zero intensity produces FLOAT_NULL in log2', async () => {
    const tsv = makeTsv([
      baseRow('P12345', 'BRCA1', '0', '0', '4096', '8192'),
    ], BASE_HEADERS);
    const df = await parseFragPipeText(tsv);
    const log2Col = df.col('log2(Sample1 MaxLFQ Intensity)')!;
    expect(log2Col.isNone(0), true);
  });

  test('primary protein ID extracted from semicolon-delimited', async () => {
    const tsv = makeTsv([
      ['sp|P12345|BRCA1_HUMAN', 'P12345;Q67890', 'BRCA1_HUMAN', 'BRCA1;BRCA2', '500',
        'Homo sapiens', 'BRCA1 protein',
        '12', '10000', '12000', '15', '20000', '24000'],
    ], BASE_HEADERS);
    const df = await parseFragPipeText(tsv);
    const primaryId = df.col('Primary Protein ID');
    expect(primaryId !== null, true);
    if (primaryId)
      expect(primaryId.get(0), 'P12345');
    const primaryGene = df.col('Primary Gene Name');
    expect(primaryGene !== null, true);
    if (primaryGene)
      expect(primaryGene.get(0), 'BRCA1');
  });

  test('proteomics.source tag set to fragpipe', async () => {
    const tsv = makeTsv([
      baseRow('P12345', 'BRCA1', '10000', '12000', '20000', '24000'),
    ], BASE_HEADERS);
    const df = await parseFragPipeText(tsv);
    expect(df.getTag('proteomics.source'), 'fragpipe');
  });

  test('intensity columns survive as raw plus log2 pair', async () => {
    const tsv = makeTsv([
      baseRow('P12345', 'BRCA1', '10000', '12000', '20000', '24000'),
    ], BASE_HEADERS);
    const df = await parseFragPipeText(tsv);
    expect(df.col('Sample1 MaxLFQ Intensity') !== null, true);
    expect(df.col('log2(Sample1 MaxLFQ Intensity)') !== null, true);
  });
});
