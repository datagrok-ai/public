import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {before, category, test, expect} from '@datagrok-libraries/utils/src/test';

import {_package} from '../package-test';
import {startAnalysis} from '../widgets/analyze-peptides';
import {PeptidesModel} from '../model';
import * as C from '../utils/constants';

category('Core', () => {
  let simpleTable: DG.DataFrame;
  let simpleActivityCol: DG.Column<number>;
  let simpleAlignedSeqCol: DG.Column<string>;
  let simpleScaledDf: DG.DataFrame;
  let simpleScaledColName: string;

  let complexTable: DG.DataFrame;
  let complexActivityCol: DG.Column<number>;
  let complexAlignedSeqCol: DG.Column<string>;
  let complexScaledDf: DG.DataFrame;
  let complexScaledColName: string;
  const alignedSeuqnceCol = 'AlignedSequence';

  let model: PeptidesModel | null = null;

  test('Start analysis: simple', async () => {
    const simpleActivityColName = 'IC50';
    simpleTable = DG.DataFrame.fromCsv(await _package.files.readAsText('aligned.csv'));
    simpleActivityCol = simpleTable.getCol(simpleActivityColName);
    simpleAlignedSeqCol = simpleTable.getCol(alignedSeuqnceCol);
    simpleAlignedSeqCol.semType = C.SEM_TYPES.ALIGNED_SEQUENCE;
    [simpleScaledDf, simpleScaledColName] =
      await PeptidesModel.scaleActivity('-lg', simpleTable, simpleActivityColName, true);

    model = await startAnalysis(
      simpleActivityCol, simpleAlignedSeqCol, simpleTable, simpleScaledDf, simpleScaledColName, _package);
    expect(model instanceof PeptidesModel, true);
    model?.setSARGridCellAt('D', '11');
    if (model != null)
      grok.shell.closeTable(model.dataFrame);
  });

  test('Start analysis: сomplex', async () => {
    const measureCategory = 'interleukin 23 receptor_h_wt_DC_HEK293_SPR_biotin_BIND Dose Response (Kd (uM):koff)';
    const complexActivityColName = 'Value';
    complexTable = DG.DataFrame.fromCsv(await _package.files.readAsText('aligned_2.csv'));
    const measuredCol: DG.Column<string> = complexTable.getCol('Measured');
    complexTable.filter.init((idx) => measuredCol.get(idx) == measureCategory);
    complexActivityCol = complexTable.getCol(complexActivityColName);
    complexAlignedSeqCol = complexTable.getCol(alignedSeuqnceCol);
    complexAlignedSeqCol.semType = C.SEM_TYPES.ALIGNED_SEQUENCE;
    [complexScaledDf, complexScaledColName] =
      await PeptidesModel.scaleActivity('-lg', complexTable, complexActivityColName, true);

    model = await startAnalysis(
      complexActivityCol, complexAlignedSeqCol, complexTable, complexScaledDf, complexScaledColName, _package);
    expect(model instanceof PeptidesModel, true);
    model?.setSARGridCellAt('-', '13');
    if (model != null)
      grok.shell.closeTable(model.dataFrame);
  });
});
