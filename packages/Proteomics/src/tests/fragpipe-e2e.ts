import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {parseFragPipeText} from '../parsers/fragpipe-parser';
import {setGroups} from '../analysis/experiment-setup';
import {medianNormalize} from '../analysis/normalization';
import {imputeMean} from '../analysis/imputation';
import {runDifferentialExpression} from '../analysis/differential-expression';
import {SEMTYPE} from '../utils/proteomics-types';
import {_package} from '../package-test';

const FIXTURE_PATH = 'demo/fragpipe-smoke-test.tsv';

const GROUPS = {
  group1: {name: 'Control', columns: ['log2(Ctrl_1 MaxLFQ Intensity)', 'log2(Ctrl_2 MaxLFQ Intensity)']},
  group2: {name: 'Treatment', columns: ['log2(Treat_1 MaxLFQ Intensity)', 'log2(Treat_2 MaxLFQ Intensity)']},
};

const UP_GENES = ['TP53', 'AKT1', 'F2'];
const DOWN_GENES = ['BRCA1', 'PRKCD'];
const UNCHANGED_GENES = ['EGFR', 'MTOR', 'HRAS'];

category('FragPipe E2E', () => {
  test('full pipeline (parse → normalize → impute → DE) on fragpipe-smoke-test.tsv', async () => {
    const text = await _package.files.readAsText(FIXTURE_PATH);
    const df = await parseFragPipeText(text);

    expect(df.rowCount, 8);
    expect(df.getTag('proteomics.source'), 'fragpipe');
    expect(df.col('Protein ID')!.semType, SEMTYPE.PROTEIN_ID);
    expect(df.col('Gene')!.semType, SEMTYPE.GENE_SYMBOL);
    expect(df.col('log2(Ctrl_1 MaxLFQ Intensity)') !== null, true);
    expect(df.col('log2(Treat_2 MaxLFQ Intensity)') !== null, true);

    setGroups(df, GROUPS);
    expect(df.getTag('proteomics.groups') !== null, true);

    const log2Cols = [...GROUPS.group1.columns, ...GROUPS.group2.columns];
    medianNormalize(df, log2Cols);
    expect(df.getTag('proteomics.normalized'), 'true');
    for (const name of log2Cols)
      expect(Math.abs(df.col(name)!.stats.med) < 1e-4, true);

    imputeMean(df, log2Cols);
    expect(df.getTag('proteomics.imputed'), 'true');

    runDifferentialExpression(df, GROUPS.group1.columns, GROUPS.group2.columns,
      'Control', 'Treatment', 1.0, 0.05);
    expect(df.getTag('proteomics.de_complete'), 'true');
    expect(df.col('log2FC') !== null, true);
    expect(df.col('p-value') !== null, true);
    expect(df.col('adj.p-value') !== null, true);
    expect(df.col('significant') !== null, true);

    const geneRow: Record<string, number> = {};
    const geneCol = df.col('Gene')!;
    for (let i = 0; i < df.rowCount; i++)
      geneRow[geneCol.get(i) as string] = i;

    const fc = df.col('log2FC')!;
    const sig = df.col('significant')!;

    for (const g of UP_GENES) {
      expect(fc.get(geneRow[g]) > 1.0, true);
      expect(sig.get(geneRow[g]), true);
    }
    for (const g of DOWN_GENES) {
      expect(fc.get(geneRow[g]) < -1.0, true);
      expect(sig.get(geneRow[g]), true);
    }
    for (const g of UNCHANGED_GENES) {
      expect(Math.abs(fc.get(geneRow[g])) < 0.5, true);
      expect(sig.get(geneRow[g]), false);
    }
  });
});
