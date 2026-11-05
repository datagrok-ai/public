import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {AnnotationCategory} from '@datagrok-libraries/bio/src/utils/macromolecule/annotations';

import {
  BUILTIN_LIABILITY_RULES, LiabilityRule, ScanLiabilitiesResult,
  applyLiabilityScanResults, createLiabilitySummaryColumn, scanLiabilities,
} from '../utils/annotations/liability-scanner';
import {
  getAnnotationColumnName, getColumnAnnotations, getRowAnnotations,
} from '../utils/annotations/annotation-manager';
import {readDataframe} from './utils';

category('liabilities', () => {
  let seqHelper: ISeqHelper;

  before(async () => {
    seqHelper = await getSeqHelper();
  });

  function rules(enabled?: (rule: LiabilityRule) => boolean): LiabilityRule[] {
    return BUILTIN_LIABILITY_RULES.map((r) => ({...r, pattern: new RegExp(r.pattern.source, 'g'),
      enabled: enabled ? enabled(r) : r.enabled}));
  }

  function fastaTable(seqs: string[]): DG.DataFrame {
    const col = DG.Column.fromStrings('seq', seqs);
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = NOTATION.FASTA;
    return DG.DataFrame.fromColumns([col]);
  }

  function scan(df: DG.DataFrame, colName: string, ruleSet: LiabilityRule[]): ScanLiabilitiesResult {
    const col = df.getCol(colName);
    return scanLiabilities(col, seqHelper.getSeqHandler(col), ruleSet);
  }

  function hitList(result: ScanLiabilitiesResult, rowIdx: number): string[] {
    return result.rowData[rowIdx].map((h) => `${h.annotationId}@${h.positionIndex}:${h.matchedMonomers}`);
  }

  test('defaultRules', async () => {
    const df = fastaTable(['QNGMW', 'ANSTDSC', 'PPPP']);
    const result = scan(df, 'seq', rules());
    expectArray(hitList(result, 0), ['deamid-ng@1:NG', 'oxid-m@3:M', 'oxid-w@4:W']);
    expectArray(hitList(result, 1), ['deamid-ns@1:NS', 'isom-ds@4:DS', 'glyco-nxst@1:NST']);
    expectArray(hitList(result, 2), []);
    expect(result.totalHits, 6);
    expectArray(result.annotations.map((a) => a.id),
      ['deamid-ng', 'deamid-ns', 'isom-ds', 'oxid-m', 'oxid-w', 'glyco-nxst']);
  });

  test('freeCysteineRule', async () => {
    const df = fastaTable(['ANSTDSC']);
    const result = scan(df, 'seq', rules((r) => r.id === 'free-cys'));
    expectArray(hitList(result, 0), ['free-cys@6:C']);
  });

  test('hitsAtReportedPositions', async () => {
    const src = await readDataframe('samples/antibodies.csv');
    const df = src.clone(DG.BitSet.create(src.rowCount, (i) => i < 40));
    await grok.data.detectSemanticTypes(df);
    const seqCol = df.getCol('AntibodyHC');
    const result = scan(df, 'AntibodyHC', rules());
    expect(result.totalHits > 0, true, 'no liabilities found in the first 40 heavy chains');
    for (let i = 0; i < df.rowCount; i++) {
      for (const hit of result.rowData[i]) {
        expect(seqCol.get(i)!.substr(hit.positionIndex, hit.matchedMonomers!.length), hit.matchedMonomers,
          `row ${i}: ${hit.annotationId} is not at ${hit.positionIndex}`);
      }
    }
  });

  test('annotationColumn', async () => {
    const df = fastaTable(['QNGMW', 'ANSTDSC', 'PPPP']);
    const seqCol = df.getCol('seq');
    const result = scan(df, 'seq', rules());
    applyLiabilityScanResults(df, seqCol, result);
    const annotCol = df.getCol(getAnnotationColumnName('seq'));
    expect(annotCol.name, '~seq_annotations');
    expect(getColumnAnnotations(seqCol).filter((a) => a.category === AnnotationCategory.Liability).length, 6);
    expect(getRowAnnotations(annotCol, 0)!.length, 3);
    expect(getRowAnnotations(annotCol, 1)!.length, 3);
    expect(getRowAnnotations(annotCol, 2), null);
    expect(df.col('seq_liability_count'), null);
  });

  test('summaryCountOxidationOnly', async () => {
    const df = fastaTable(['QNGMW', 'MWMDSW', 'PPPP']);
    const seqCol = df.getCol('seq');
    const all = scan(df, 'seq', rules());
    const oxidation = scan(df, 'seq', rules((r) => r.ruleCategory === 'oxidation'));
    const countCol = createLiabilitySummaryColumn(df, seqCol, oxidation);
    expect(countCol.name, 'seq_liability_count');
    expect(countCol.type, DG.TYPE.INT);
    expect(countCol.stats.missingValueCount, 0);
    expectArray(countCol.toList(), [2, 4, 0]);
    expect(countCol.stats.sum < all.totalHits, true, `${countCol.stats.sum} oxidation hits of ${all.totalHits}`);
  });
});
