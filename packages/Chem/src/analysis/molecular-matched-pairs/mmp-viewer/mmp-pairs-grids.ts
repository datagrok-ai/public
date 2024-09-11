import * as DG from 'datagrok-api/dg';
import {ILineSeries} from '@datagrok-libraries/utils/src/render-lines-on-sp';
import {SUBSTRUCT_COL} from '../../../constants';
import {MMP_NAMES, columnsDescriptions} from './mmp-constants';
import {PaletteCodes} from './palette';
import {createColWithDescription} from './mmp-generations';
import {MMPA} from '../mmp-analysis/mmpa';

export function getMmpPairsGrids(mmpa: MMPA, palette: PaletteCodes): {
  activityMeanNames: Array<string>,
  fpGrid: DG.Grid,
  mmpGrid: DG.Grid,
  lines: ILineSeries,
  linesIdxs: Uint32Array,
  linesActivityCorrespondance: Uint32Array
} {
  const activityMeanNames = getActivitMeanNames(mmpa);
  const fpGrid = getFragmetsPairsGrid(activityMeanNames, mmpa);
  const mmpGrid = getMatchedPairsGrid(mmpa);
  const lines = createLines(mmpa, palette);

  return {activityMeanNames, linesIdxs: lines.linesIdxs, fpGrid, mmpGrid,
    lines: lines.lines, linesActivityCorrespondance: lines.linesActivityCorrespondance};
}

function getActivitMeanNames(mmpa: MMPA) {
  const activityMeanNames = Array<string>(mmpa.initData.activitiesCount);
  for (let i = 0; i < mmpa.initData.activitiesCount; i++) {
    const name = MMP_NAMES.MEANDIFF + ' ' + mmpa.initData.activitiesNames[i];
    activityMeanNames[i] = name;
  }

  return activityMeanNames;
}

function getFragmetsPairsGrid(activityMeanNames: string[], mmpa: MMPA) : DG.Grid {
  const fromCol = createColWithDescription('string', MMP_NAMES.FROM, mmpa.rulesBased.fromFrag);
  const toCol = createColWithDescription('string', MMP_NAMES.TO, mmpa.rulesBased.toFrag);
  const occasionsCol = DG.Column.fromInt32Array(MMP_NAMES.PAIRS, mmpa.rulesBased.occasions);
  occasionsCol.setTag('description', columnsDescriptions[MMP_NAMES.PAIRS]);
  fromCol.semType = DG.SEMTYPE.MOLECULE;
  toCol.semType = DG.SEMTYPE.MOLECULE;
  occasionsCol.semType = DG.TYPE.INT;

  const fpCols = [fromCol, toCol, occasionsCol];

  for (let i = 0; i < activityMeanNames.length; i++)
    fpCols.push(DG.Column.fromFloat32Array(activityMeanNames[i], mmpa.rulesBased.meanDiffs[i]));

  const colorsPal = new Int32Array(fromCol.length);
  for (let i = 0; i < fromCol.length; i++)
    colorsPal[i] = i < fromCol.length ? 0 : 1;

  const colorCol = DG.Column.fromInt32Array(MMP_NAMES.COLOR, colorsPal);
  fpCols.push(colorCol);

  const fpDf = DG.DataFrame.fromColumns(fpCols);
  const fpGrid = fpDf.plot.grid();
  fpGrid.col(MMP_NAMES.COLOR)!.visible = false;
  return fpGrid;
}

function getMatchedPairsGrid(mmpa: MMPA) : DG.Grid {
  const pairsFromCol = createColWithDescription('string', MMP_NAMES.FROM, mmpa.allCasesBased.molFrom);
  const pairsToCol = createColWithDescription('string', MMP_NAMES.TO, mmpa.allCasesBased.molTo);
  const structureDiffFromCol =
    DG.Column.fromType('object', MMP_NAMES.STRUCT_DIFF_FROM_NAME, mmpa.allCasesBased.molFrom.length);
  const structureDiffToCol =
    DG.Column.fromType('object', MMP_NAMES.STRUCT_DIFF_TO_NAME, mmpa.allCasesBased.molFrom.length);
  const pairNumberCol = DG.Column.fromInt32Array(MMP_NAMES.PAIRNUM, mmpa.allCasesBased.pairNum);
  const pairNumberFromCol = DG.Column.fromInt32Array(MMP_NAMES.PAIRNUM_FROM, mmpa.allCasesBased.molNumFrom);
  const pairNumberToCol = DG.Column.fromInt32Array(MMP_NAMES.PAIRNUM_TO, mmpa.allCasesBased.molNumTo);

  const pairsFromSmilesCol = DG.Column.fromStrings(MMP_NAMES.SMI1, mmpa.allCasesBased.pairsFromSmiles);
  const pairsToSmilesCol = DG.Column.fromStrings(MMP_NAMES.SMI2, mmpa.allCasesBased.pairsToSmiles);
  const ruleNumCol = DG.Column.fromInt32Array(MMP_NAMES.RULENUM, mmpa.allCasesBased.ruleNum);

  pairsFromSmilesCol.semType = DG.SEMTYPE.MOLECULE;
  pairsToSmilesCol.semType = DG.SEMTYPE.MOLECULE;
  pairsFromCol.semType = DG.SEMTYPE.MOLECULE;
  pairsToCol.semType = DG.SEMTYPE.MOLECULE;
  pairsFromCol.temp[SUBSTRUCT_COL] = MMP_NAMES.STRUCT_DIFF_FROM_NAME;
  pairsToCol.temp[SUBSTRUCT_COL] = MMP_NAMES.STRUCT_DIFF_TO_NAME;

  const allTransformationsCols = [pairsFromCol, pairsToCol,
    structureDiffFromCol, structureDiffToCol,
    pairNumberCol, pairNumberFromCol, pairNumberToCol,
    pairsFromSmilesCol, pairsToSmilesCol, ruleNumCol];

  for (let i = 0; i < mmpa.initData.activitiesCount; i++) {
    const name = MMP_NAMES.DIFF + ' ' + mmpa.initData.activitiesNames[i];
    allTransformationsCols.push(DG.Column.fromFloat32Array(name, mmpa.allCasesBased.diffs[i]));
  }

  const pairedTransformations = DG.DataFrame.fromColumns(allTransformationsCols);
  return pairedTransformations.plot.grid();
}

function createLines(mmpa: MMPA, palette: PaletteCodes) : {
    linesIdxs: Uint32Array,
    lines: ILineSeries,
    linesActivityCorrespondance: Uint32Array
  } {
  let allCount = 0;
  for (let i = 0; i < mmpa.initData.activitiesCount; i++)
    allCount += mmpa.allCasesBased.activityPairsIdxs[i].trueCount();

  const pointsFrom = new Uint32Array(allCount);
  const pointsTo = new Uint32Array(allCount);
  const linesIdxs = new Uint32Array(allCount);
  const colors = new Array<string>(allCount);
  const linesActivityCorrespondance = new Uint32Array(allCount);

  let pairsCounter = 0;
  for (let j = 0; j < mmpa.initData.activitiesCount; j++) {
    for (let i = -1; (i = mmpa.allCasesBased.activityPairsIdxs[j].findNext(i)) !== -1;) {
      pointsFrom[pairsCounter] = mmpa.allCasesBased.molNumFrom[i];
      pointsTo[pairsCounter] = mmpa.allCasesBased.molNumTo[i];
      linesIdxs[pairsCounter] = i;
      colors[pairsCounter] = palette.rgbCut[j];
      linesActivityCorrespondance[pairsCounter] = j;
      pairsCounter++;
    }
  }

  const lines: ILineSeries = {
    from: pointsFrom,
    to: pointsTo,
    drawArrows: true,
    colors: colors,
    arrowSize: 10,
    width: 0.5,
  };

  return {linesIdxs, lines, linesActivityCorrespondance};
}
