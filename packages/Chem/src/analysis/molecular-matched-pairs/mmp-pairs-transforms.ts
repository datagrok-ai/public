import * as DG from 'datagrok-api/dg';
import {ILineSeries} from '@datagrok-libraries/utils/src/render-lines-on-sp';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {MmpRules} from './mmp-fragments';
import {SUBSTRUCT_COL} from '../../constants';
import {MMP_COLNAME_FROM, MMP_COLNAME_TO, MMP_COLNAME_PAIRS, MMP_COL_PAIRNUM, MMP_COL_PAIRNUM_FROM,
  MMP_COL_PAIRNUM_TO, MMP_COLNAME_MEANDIFF, MMP_COLNAME_DIFF, MMP_STRUCT_DIFF_FROM_NAME, MMP_STRUCT_DIFF_TO_NAME,
  MMP_COLOR} from './mmp-constants';

export function getMmpActivityPairsAndTransforms(molecules: DG.Column, activities: DG.ColumnList, mmpRules: MmpRules,
  allCasesNumber: number, palette: Array<string>): {
  maxActs: number[],
  diffs: Array<Float32Array>,
  activityMeanNames: Array<string>,
  linesIdxs: Uint32Array,
  allPairsGrid: DG.Grid,
  casesGrid: DG.Grid,
  lines: ILineSeries,
  linesActivityCorrespondance: Uint32Array
} {
  const variates = activities.length;

  const maxActs = new Array<number>(variates);
  for (let i = 0; i < variates; i++)
    maxActs[i] = 0;
  const fromFrag = new Array<string>(mmpRules.rules.length);
  const toFrag = new Array<string>(mmpRules.rules.length);
  const occasions = new Int32Array(mmpRules.rules.length);
  const meanDiff = new Array<Float32Array>(variates);
  for (let i = 0; i < variates; i++)
    meanDiff[i] = new Float32Array(mmpRules.rules.length);

  const molFrom = new Array<string>(allCasesNumber);
  const molTo = new Array<string>(allCasesNumber);
  const pairNum = new Int32Array(allCasesNumber);
  const molNumFrom = new Int32Array(allCasesNumber);
  const molNumTo = new Int32Array(allCasesNumber);
  const diffs = new Array<Float32Array>(variates);
  for (let i = 0; i < variates; i++)
    diffs[i] = new Float32Array(allCasesNumber);

  const pairsFromSmiles = new Array<string>(allCasesNumber);
  const pairsToSmiles = new Array<string>(allCasesNumber);
  const ruleNum = new Int32Array(allCasesNumber);

  const activityPairsIdxs = new Array<BitArray>(variates);
  for (let i = 0; i < variates; i++)
    activityPairsIdxs[i] = new BitArray(allCasesNumber);

  //set activity differences
  let pairIdx = 0;
  for (let i = 0; i < mmpRules.rules.length; i++) {
    fromFrag[i] = mmpRules.smilesFrags[mmpRules.rules[i].smilesRule1];
    toFrag[i] = mmpRules.smilesFrags[mmpRules.rules[i].smilesRule2];
    occasions[i] = mmpRules.rules[i].pairs.length;

    const mean = new Float32Array(variates);
    for (let j = 0; j < occasions[i]; j++) {
      const idx1 = mmpRules.rules[i].pairs[j].firstStructure;
      const idx2 = mmpRules.rules[i].pairs[j].secondStructure;

      molFrom[pairIdx] = molecules.get(idx1);
      molTo[pairIdx] = molecules.get(idx2);

      for (let k = 0; k < variates; k++) {
        //TODO: make more efficient
        const col = activities.byIndex(k);
        const diff = col.get(idx2) - col.get(idx1);
        diffs[k][pairIdx] = diff;
        if (diff > 0) {
          if (diff > maxActs[k])
            maxActs[k] = diff;
          activityPairsIdxs[k].setBit(pairIdx, true, false);
        }

        mean[k] += diff;
      }

      molNumFrom[pairIdx] = idx1;
      molNumTo[pairIdx] = idx2;
      pairNum[pairIdx] = pairIdx;
      pairsFromSmiles[pairIdx] = mmpRules.smilesFrags[mmpRules.rules[i].smilesRule1];
      pairsToSmiles[pairIdx] = mmpRules.smilesFrags[mmpRules.rules[i].smilesRule2];
      ruleNum[pairIdx] = i;

      pairIdx++;
    }

    for (let k = 0; k < variates; k++) {
      mean[k] /= occasions[i];
      meanDiff[k][i] = mean[k];
    }
  }

  //creating fragments grid
  const fromCol = DG.Column.fromList('string', MMP_COLNAME_FROM, fromFrag);
  const toCol = DG.Column.fromList('string', MMP_COLNAME_TO, toFrag);
  const occasionsCol = DG.Column.fromInt32Array(MMP_COLNAME_PAIRS, occasions);
  fromCol.semType = DG.SEMTYPE.MOLECULE;
  toCol.semType = DG.SEMTYPE.MOLECULE;
  occasionsCol.semType = DG.TYPE.INT;

  const allPairsCols = [fromCol, toCol, occasionsCol];
  const activityMeanNames = Array<string>(variates);
  for (let i = 0; i < variates; i++) {
    const name = MMP_COLNAME_MEANDIFF + ' ' + activities.byIndex(i).name;
    activityMeanNames[i] = name;
    allPairsCols.push(DG.Column.fromFloat32Array(name, meanDiff[i]));
  }

  //TODO: make compatible with trellis
  const colorsPal = new Int32Array(fromCol.length);
  for (let i = 0; i < fromCol.length; i++)
    colorsPal[i] = i < fromCol.length ? 0 : 1;

  const colorCol = DG.Column.fromInt32Array(MMP_COLOR, colorsPal);
  allPairsCols.push(colorCol);

  const dfAllPairs = DG.DataFrame.fromColumns(allPairsCols);
  const allPairsGrid = dfAllPairs.plot.grid();

  //creating cases grid
  const pairsFromCol = DG.Column.fromStrings(MMP_COLNAME_FROM, molFrom);
  const pairsToCol = DG.Column.fromStrings(MMP_COLNAME_TO, molTo);
  const structureDiffFromCol = DG.Column.fromType('object', MMP_STRUCT_DIFF_FROM_NAME, molFrom.length);
  const structureDiffToCol = DG.Column.fromType('object', MMP_STRUCT_DIFF_TO_NAME, molFrom.length);
  const pairNumberCol = DG.Column.fromInt32Array(MMP_COL_PAIRNUM, pairNum);
  const pairNumberFromCol = DG.Column.fromInt32Array(MMP_COL_PAIRNUM_FROM, molNumFrom);
  const pairNumberToCol = DG.Column.fromInt32Array(MMP_COL_PAIRNUM_TO, molNumTo);

  const pairsFromSmilesCol = DG.Column.fromStrings('~smi1', pairsFromSmiles);
  const pairsToSmilesCol = DG.Column.fromStrings('~smi2', pairsToSmiles);
  const ruleNumCol = DG.Column.fromInt32Array('~ruleNum', ruleNum);

  pairsFromSmilesCol.semType = DG.SEMTYPE.MOLECULE;
  pairsToSmilesCol.semType = DG.SEMTYPE.MOLECULE;
  pairsFromCol.semType = DG.SEMTYPE.MOLECULE;
  pairsToCol.semType = DG.SEMTYPE.MOLECULE;
  pairsFromCol.temp[SUBSTRUCT_COL] = MMP_STRUCT_DIFF_FROM_NAME;
  pairsToCol.temp[SUBSTRUCT_COL] = MMP_STRUCT_DIFF_TO_NAME;

  const allTransformationsCols = [pairsFromCol, pairsToCol,
    structureDiffFromCol, structureDiffToCol,
    pairNumberCol, pairNumberFromCol, pairNumberToCol,
    pairsFromSmilesCol, pairsToSmilesCol, ruleNumCol];

  for (let i = 0; i < variates; i++) {
    const name = MMP_COLNAME_DIFF + ' ' + activities.byIndex(i).name;
    allTransformationsCols.push(DG.Column.fromFloat32Array(name, diffs[i]));
  }

  const pairedTransformations = DG.DataFrame.fromColumns(allTransformationsCols);
  const casesGrid = pairedTransformations.plot.grid();

  //creating lines for rendering
  let allCount = 0;
  for (let i = 0; i < variates; i++)
    allCount += activityPairsIdxs[i].trueCount();

  const pointsFrom = new Uint32Array(allCount);
  const pointsTo = new Uint32Array(allCount);
  const linesIdxs = new Uint32Array(allCount);
  const colors = new Array<string>(allCount);
  const linesActivityCorrespondance = new Uint32Array(allCount);

  let pairsCounter = 0;
  for (let j = 0; j < variates; j++) {
    for (let i = -1; (i = activityPairsIdxs[j].findNext(i)) !== -1;) {
      pointsFrom[pairsCounter] = molNumFrom[i];
      pointsTo[pairsCounter] = molNumTo[i];
      linesIdxs[pairsCounter] = i;
      colors[pairsCounter] = palette[j].replace('rgb(', '').replace(')', '');//j == 0 ? '60,177,115' : '255,0,0';
      linesActivityCorrespondance[pairsCounter] = j;
      pairsCounter++;
    }
  }

  const lines: ILineSeries = {
    from: pointsFrom,
    to: pointsTo,
    drawArrows: true,
    opacity: 0.5,
    colors: colors,
    arrowSize: 10,
  };

  return {maxActs, diffs, activityMeanNames, linesIdxs, allPairsGrid, casesGrid, lines, linesActivityCorrespondance};
}
