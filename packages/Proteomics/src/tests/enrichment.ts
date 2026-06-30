import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {buildEnrichmentDf, countSignificantProteins, splitGenesByDirection,
  ORGANISM_LIST, GostResult} from '../analysis/enrichment';
import {wireEnrichmentToVolcano} from '../viewers/enrichment-viewers';
import {SEMTYPE} from '../utils/proteomics-types';

// --- Helpers ---

function makeMockResults(): GostResult[] {
  return [
    {
      native: 'GO:0006915',
      name: 'apoptotic process',
      source: 'GO:BP',
      p_value: 0.001,
      significant: true,
      term_size: 500,
      query_size: 50,
      intersection_size: 10,
      effective_domain_size: 20000,
      precision: 0.2,
      recall: 0.02,
      intersections: [['IEA'], ['IDA'], [], ['TAS']],
    },
    {
      native: 'KEGG:04210',
      name: 'Apoptosis',
      source: 'KEGG',
      p_value: 0.03,
      significant: true,
      term_size: 120,
      query_size: 50,
      intersection_size: 5,
      effective_domain_size: 20000,
      precision: 0.1,
      recall: 0.042,
      intersections: [[], ['evidence'], [], []],
    },
    {
      native: 'REAC:R-HSA-109581',
      name: 'Apoptosis',
      source: 'REAC',
      p_value: 0.12,
      significant: false,
      term_size: 80,
      query_size: 50,
      intersection_size: 3,
      effective_domain_size: 20000,
      precision: 0.06,
      recall: 0.0375,
      intersections: [['IEA'], [], [], []],
    },
  ];
}

// --- Tests ---

category('Enrichment', () => {
  test('buildEnrichmentDf creates correct schema', async () => {
    const results = makeMockResults();
    const queryGenes = ['TP53', 'BRCA1', 'EGFR', 'MYC'];
    const df = buildEnrichmentDf(results, queryGenes);

    // 7 data cols (Source/Term ID/Term Name/FDR/Gene Count/Gene Ratio/Intersection)
    // + 1 Significant. g:GOSt with significance_threshold_method='fdr' returns the
    // FDR-corrected p-value with no separate raw p-value, so the schema exposes
    // only FDR — see buildEnrichmentDf in src/analysis/enrichment.ts.
    expect(df.columns.length, 8);
    expect(df.col('Source') !== null, true);
    expect(df.col('Term ID') !== null, true);
    expect(df.col('Term Name') !== null, true);
    expect(df.col('FDR') !== null, true);
    expect(df.col('Gene Count') !== null, true);
    expect(df.col('Gene Ratio') !== null, true);
    expect(df.col('Intersection') !== null, true);
    expect(df.col('Significant') !== null, true);
    expect(df.rowCount, 3);
  });

  test('buildEnrichmentDf sets enrichment tag', async () => {
    const results = makeMockResults();
    const df = buildEnrichmentDf(results, ['TP53', 'BRCA1', 'EGFR', 'MYC']);
    expect(df.getTag('proteomics.enrichment'), 'true');
  });

  test('buildEnrichmentDf marks significant terms', async () => {
    const results = makeMockResults(); // first two < 0.05, third > 0.05
    const df = buildEnrichmentDf(results, ['TP53', 'BRCA1', 'EGFR', 'MYC']);
    const sigCol = df.col('Significant')!;
    expect(sigCol.get(0), true);  // p=0.001
    expect(sigCol.get(1), true);  // p=0.03
    expect(sigCol.get(2), false); // p=0.12
  });

  test('buildEnrichmentDf extracts intersection genes', async () => {
    const queryGenes = ['TP53', 'BRCA1', 'EGFR', 'MYC'];
    const results: GostResult[] = [{
      native: 'GO:0006915',
      name: 'apoptotic process',
      source: 'GO:BP',
      p_value: 0.01,
      significant: true,
      term_size: 500,
      query_size: 4,
      intersection_size: 2,
      effective_domain_size: 20000,
      precision: 0.5,
      recall: 0.004,
      intersections: [['IEA'], [], ['TAS'], []],
    }];
    const df = buildEnrichmentDf(results, queryGenes);
    const intersectionVal = df.col('Intersection')!.get(0) as string;
    expect(intersectionVal, 'TP53, EGFR');
  });

  test('countSignificantProteins counts correctly', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat32Array('log2FC', new Float32Array([2.0, 0.5, -1.5, 3.0])),
      DG.Column.fromFloat32Array('adjP', new Float32Array([0.01, 0.001, 0.1, 0.03])),
    ]);
    const result = countSignificantProteins(
      df, 1.0, 0.05, df.col('log2FC')!, df.col('adjP')!,
    );
    // |2.0|>=1 && 0.01<=0.05 => yes; |0.5|<1 => no; |-1.5|>=1 && 0.1>0.05 => no; |3.0|>=1 && 0.03<=0.05 => yes
    expect(result.significant, 2);
    expect(result.total, 4);
  });

  test('countSignificantProteins handles null values', async () => {
    const fcCol = DG.Column.fromFloat32Array('log2FC', new Float32Array([2.0, 0, -1.5]));
    const pCol = DG.Column.fromFloat32Array('adjP', new Float32Array([0.01, 0, 0.03]));
    // Set index 1 to null
    fcCol.set(1, DG.FLOAT_NULL);
    pCol.set(1, DG.FLOAT_NULL);

    const df = DG.DataFrame.fromColumns([fcCol, pCol]);
    const result = countSignificantProteins(df, 1.0, 0.05, fcCol, pCol);
    // row 0: |2.0|>=1 && 0.01<=0.05 => yes; row 1: null => skip; row 2: |-1.5|>=1 && 0.03<=0.05 => yes
    expect(result.significant, 2);
    expect(result.total, 3); // all rows counted in total
  });

  // --- R2: directional split (13-02) ---

  test('splitGenesByDirection splits significant genes by fc sign, shared background', async () => {
    const geneForRow = new Map<number, string>([[0, 'A'], [1, 'B'], [2, 'C'], [3, 'D']]);
    const fc = new Float32Array([2.0, -2.0, 0.5, -3.0]);
    const p = new Float32Array([0.01, 0.01, 0.001, 0.2]);
    const {upGenes, downGenes, background} = splitGenesByDirection(geneForRow, fc, p, 1.0, 0.05);
    // A: |2.0|>=1 & .01<=.05 → up. B: |-2.0|>=1 & .01<=.05 → down.
    // C: |0.5|<1 → neither. D: .2>.05 → neither.
    expect(upGenes.length, 1);
    expect(upGenes[0], 'A');
    expect(downGenes.length, 1);
    expect(downGenes[0], 'B');
    // Background = all detected genes (both directional queries share it).
    expect(background.length, 4);
    expect(background.includes('A') && background.includes('B') &&
      background.includes('C') && background.includes('D'), true);
    // Up/Down are disjoint.
    expect(upGenes.some((g) => downGenes.includes(g)), false);
  });

  test('splitGenesByDirection skips null fc/p rows', async () => {
    const geneForRow = new Map<number, string>([[0, 'A'], [1, 'B']]);
    const fc = new Float32Array([2.0, 0]);
    const p = new Float32Array([0.01, 0]);
    fc[1] = DG.FLOAT_NULL;
    p[1] = DG.FLOAT_NULL;
    const {upGenes, downGenes, background} = splitGenesByDirection(geneForRow, fc, p, 1.0, 0.05);
    expect(upGenes.length, 1);
    expect(downGenes.length, 0);
    expect(background.length, 2); // both rows still counted as detected background
  });

  test('buildEnrichmentDf adds Direction column when label given', async () => {
    const df = buildEnrichmentDf(makeMockResults(), ['TP53', 'BRCA1', 'EGFR', 'MYC'], 0.05, 'Up');
    const dir = df.col('Direction');
    expect(dir !== null, true);
    expect(dir!.get(0), 'Up');
    expect(dir!.get(1), 'Up');
    // Phase-9 schema preserved + Direction = 9 columns.
    expect(df.columns.length, 9);
    expect(df.col('Intersection') !== null, true);
    expect(df.getTag('proteomics.enrichment'), 'true');
  });

  test('buildEnrichmentDf omits Direction column when no label', async () => {
    const df = buildEnrichmentDf(makeMockResults(), ['TP53', 'BRCA1', 'EGFR', 'MYC']);
    expect(df.col('Direction'), null);
    expect(df.columns.length, 8);
  });

  test('wireEnrichmentToVolcano selects matching proteins, source DF unmutated', async () => {
    const enrichDf = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Intersection', ['TP53, EGFR']),
      DG.Column.fromStrings('Direction', ['Up']),
    ]);
    const geneCol = DG.Column.fromStrings('Gene Symbol', ['TP53', 'EGFR', 'MYC']);
    geneCol.semType = SEMTYPE.GENE_SYMBOL;
    const proteinDf = DG.DataFrame.fromColumns([geneCol]);

    const sub = wireEnrichmentToVolcano(enrichDf, proteinDf);
    enrichDf.currentRowIdx = 0;
    // TP53 + EGFR selected (from Intersection); MYC not.
    expect(proteinDf.selection.trueCount, 2);
    expect(proteinDf.selection.get(2), false);
    // Cross-link must not mutate either source frame's shape.
    expect(enrichDf.rowCount, 1);
    expect(proteinDf.rowCount, 3);
    sub.unsubscribe();
  });

  test('ORGANISM_LIST contains 9 entries with valid codes', async () => {
    expect(ORGANISM_LIST.length, 9);
    expect(ORGANISM_LIST[0].display, 'Homo sapiens (Human)');
    expect(ORGANISM_LIST[0].code, 'hsapiens');
    for (const org of ORGANISM_LIST) {
      expect(org.display.length > 0, true);
      expect(org.code.length > 0, true);
    }
  });
});
