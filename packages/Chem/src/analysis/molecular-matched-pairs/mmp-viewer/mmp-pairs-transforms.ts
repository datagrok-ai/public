import * as DG from 'datagrok-api/dg';
import {ILineSeries} from '@datagrok-libraries/utils/src/render-lines-on-sp';
import {MmpInput} from './mmp-viewer';
import {SUBSTRUCT_COL} from '../../../constants';
import {MMP_NAMES, columnsDescriptions} from './mmp-constants';
import {PaletteCodes} from './palette';
import {createColWithDescription} from '../mmp-generations';
import {MMPA} from '../mmp-analysis/mmpa';

function createLines(mmpa: MMPA, palette: PaletteCodes) : {
    linesIdxs: Uint32Array,
    lines: ILineSeries,
    linesActivityCorrespondance: Uint32Array
  } {
  let allCount = 0;
  for (let i = 0; i < mmpa.allCasesBased.variates; i++)
    allCount += mmpa.allCasesBased.activityPairsIdxs[i].trueCount();

  const pointsFrom = new Uint32Array(allCount);
  const pointsTo = new Uint32Array(allCount);
  const linesIdxs = new Uint32Array(allCount);
  const colors = new Array<string>(allCount);
  const linesActivityCorrespondance = new Uint32Array(allCount);

  let pairsCounter = 0;
  for (let j = 0; j < mmpa.allCasesBased.variates; j++) {
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

function getTransFragmetsGrid(activities: DG.ColumnList, mmpa: MMPA) : [string[], DG.Grid] {
  const fromCol = createColWithDescription('string', MMP_NAMES.FROM, mmpa.rulesBased.fromFrag);
  const toCol = createColWithDescription('string', MMP_NAMES.TO, mmpa.rulesBased.toFrag);
  const occasionsCol = DG.Column.fromInt32Array(MMP_NAMES.PAIRS, mmpa.rulesBased.occasions);
  occasionsCol.setTag('description', columnsDescriptions[MMP_NAMES.PAIRS]);
  fromCol.semType = DG.SEMTYPE.MOLECULE;
  toCol.semType = DG.SEMTYPE.MOLECULE;
  occasionsCol.semType = DG.TYPE.INT;

  const allPairsCols = [fromCol, toCol, occasionsCol];
  const activityMeanNames = Array<string>(mmpa.allCasesBased.variates);
  for (let i = 0; i < mmpa.allCasesBased.variates; i++) {
    const name = MMP_NAMES.MEANDIFF + ' ' + activities.byIndex(i).name;
    activityMeanNames[i] = name;
    allPairsCols.push(DG.Column.fromFloat32Array(name, mmpa.rulesBased.meanDiffs[i]));
  }

  //TODO: make compatible with trellis
  const colorsPal = new Int32Array(fromCol.length);
  for (let i = 0; i < fromCol.length; i++)
    colorsPal[i] = i < fromCol.length ? 0 : 1;

  const colorCol = DG.Column.fromInt32Array(MMP_NAMES.COLOR, colorsPal);
  allPairsCols.push(colorCol);

  const dfAllPairs = DG.DataFrame.fromColumns(allPairsCols);
  const grid = dfAllPairs.plot.grid();
  grid.col(MMP_NAMES.COLOR)!.visible = false;
  return [activityMeanNames, grid];
}

function getTransPairsGrid(activities: DG.ColumnList, mmpa: MMPA) : DG.Grid {
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

  for (let i = 0; i < mmpa.allCasesBased.variates; i++) {
    const name = MMP_NAMES.DIFF + ' ' + activities.byIndex(i).name;
    allTransformationsCols.push(DG.Column.fromFloat32Array(name, mmpa.allCasesBased.diffs[i]));
  }

  const pairedTransformations = DG.DataFrame.fromColumns(allTransformationsCols);
  return pairedTransformations.plot.grid();
}

export function getMmpActivityPairsAndTransforms(mmpInput: MmpInput, mmpa: MMPA, palette: PaletteCodes): {
  activityMeanNames: Array<string>,
  transFragmentsGrid: DG.Grid,
  transPairsGrid: DG.Grid,
  lines: ILineSeries,
  linesIdxs: Uint32Array,
  linesActivityCorrespondance: Uint32Array
} {
  const [activityMeanNames, transFragmentsGrid] = getTransFragmetsGrid(mmpInput.activities, mmpa);
  const transPairsGrid = getTransPairsGrid( mmpInput.activities, mmpa);
  const lines = createLines(mmpa, palette);

  return {activityMeanNames, linesIdxs: lines.linesIdxs, transFragmentsGrid, transPairsGrid,
    lines: lines.lines, linesActivityCorrespondance: lines.linesActivityCorrespondance};
}
