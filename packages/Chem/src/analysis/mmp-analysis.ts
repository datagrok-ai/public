import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {drawMoleculeToCanvas, getRdKitModule, getRdKitService, getUncommonAtomsAndBonds}
  from '../utils/chem-common-rdkit';

import {DimReductionMethods} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {BitArrayMetrics, BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {chemSpace} from './chem-space';
import $ from 'cash-dom';
import {ISubstruct} from '../rendering/rdkit-cell-renderer';
import {SUBSTRUCT_COL} from '../constants';
import {ILineSeries, MouseOverLineEvent, ScatterPlotCurrentLineStyle, ScatterPlotLinesRenderer}
  from '@datagrok-libraries/utils/src/render-lines-on-sp';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {debounceTime} from 'rxjs/operators';
import {getSigFigs} from '../utils/chem-common';
import {convertMolNotation} from '../package';
import {pickTextColorBasedOnBgColor} from '../utils/ui-utils';

const MMP_COLNAME_FROM = 'From';
const MMP_COLNAME_TO = 'To';
const MMP_COLNAME_PAIRS = 'Pairs';
const MMP_COL_PAIRNUM = '~PairNum';
const MMP_COL_PAIRNUM_FROM = '~PairNumFrom';
const MMP_COL_PAIRNUM_TO = '~PairNumTo';
const MMP_COLNAME_MEANDIFF = 'Mean Difference';
const MMP_COLNAME_DIFF = 'Difference';
const MMP_COLNAME_CHEMSPACE_X = '~X';
const MMP_COLNAME_CHEMSPACE_Y = '~Y';
const MMP_VIEW_NAME = 'MMP Analysis';
const MMP_TAB_TRANSFORMATIONS = 'Transformations';
const MMP_TAB_FRAGMENTS = 'Fragments';
const MMP_TAB_CLIFFS = 'Cliffs';
const MMP_STRUCT_DIFF_FROM_NAME = '~structDiffFrom';
const MMP_STRUCT_DIFF_TO_NAME = '~structDiffTo';

type MmpRules = {
  rules: {
    smilesRule1: number,
    smilesRule2: number,
    pairs: {firstStructure: number, secondStructure: number}[]
  } [],
  smilesFrags: string[]
};

let currentTab = '';
let lastSelectedPair: number | null = null;
let propColNames: string[] = [];
let propLabelsDiv: HTMLDivElement = ui.div();

async function getMmpFrags(molecules: DG.Column): Promise<[string, string][][]> {
  const service = await getRdKitService();
  const res = await service.mmpGetFragments(molecules.toList());
  return res;
}

//returns mmp rules and number of cases
function getMmpRules(frags: [string, string][][]): [MmpRules, number] {
  const mmpRules: MmpRules = {rules: [], smilesFrags: []};
  const dim = frags.length;
  let ruleCounter = 0;
  let allCasesCounter = 0;

  for (let i = 0; i < dim; i++) {
    const dim1 = frags[i].length;
    for (let j = i + 1; j < dim; j++) {
      const dim2 = frags[j].length;
      let core = '';
      let r1 = ''; // molecule minus core for first molecule in pair
      let r2 = ''; // molecule minus core for second molecule in pair

      //here we get the best possible fragment pair
      //TODO: do not process molecular pairs with low similarity
      for (let p1 = 0; p1 < dim1; p1++) {
        for (let p2 = 0; p2 < dim2; p2++) {
          if (frags[i][p1][0] == frags[j][p2][0]) {
            const newCore = frags[i][p1][0];
            if (newCore.length > core.length) {
              core = newCore;
              r1 = frags[i][p1][1];
              r2 = frags[j][p2][1];
            }
          }
        }
      }

      if (core === '' || r1.length / core.length > 0.4 || r2.length / core.length > 0.4)
        continue;

      let ruleSmiles1 = mmpRules.smilesFrags.indexOf(r1);
      let ruleSmiles2 = mmpRules.smilesFrags.indexOf(r2);
      let ruleIndexStraight = -1;
      let ruleIndexInverse = -1;

      for (let ind = 0; ind < mmpRules.rules.length; ind++) {
        if (mmpRules.rules[ind].smilesRule1 == ruleSmiles1 && mmpRules.rules[ind].smilesRule2 == ruleSmiles2)
          ruleIndexStraight = ind;
        if (mmpRules.rules[ind].smilesRule1 == ruleSmiles2 && mmpRules.rules[ind].smilesRule2 == ruleSmiles1)
          ruleIndexInverse = ind;
      }

      if (ruleSmiles1 == -1) {
        mmpRules.smilesFrags.push(r1);
        ruleSmiles1 = mmpRules.smilesFrags.length -1;
      }
      if (ruleSmiles2 == -1) {
        mmpRules.smilesFrags.push(r2);
        ruleSmiles2 = mmpRules.smilesFrags.length -1;
      }

      const indxFirst = ruleSmiles1 < ruleSmiles2;
      if (ruleIndexStraight == -1) {
        mmpRules.rules.push({
          smilesRule1: indxFirst ? ruleSmiles1: ruleSmiles2,
          smilesRule2: indxFirst ? ruleSmiles2: ruleSmiles1,
          pairs: [],
        });
        mmpRules.rules[ruleCounter].pairs.push({firstStructure: indxFirst ? i : j, secondStructure: indxFirst ? j : i});
        ruleCounter++;
        mmpRules.rules.push({
          smilesRule1: indxFirst ? ruleSmiles2: ruleSmiles1,
          smilesRule2: indxFirst ? ruleSmiles1: ruleSmiles2,
          pairs: [],
        });
        mmpRules.rules[ruleCounter].pairs.push({firstStructure: indxFirst ? j : i, secondStructure: indxFirst ? i : j});
        ruleCounter++;
        allCasesCounter += 2;
      } else {
        mmpRules.rules[ruleIndexStraight].pairs
          .push({firstStructure: i, secondStructure: j});
        mmpRules.rules[ruleIndexInverse].pairs
          .push({firstStructure: j, secondStructure: i});
        allCasesCounter += 2;
      }
    }
  }

  return [mmpRules, allCasesCounter];
}

function getMmpActivityPairsAndTransforms(molecules: DG.Column, activities: DG.ColumnList, mmpRules: MmpRules,
  allCasesNumber: number): {
  maxAct: number,
  diffs: Array<Float32Array>,
  activityMeanNames: Array<string>,
  linesIdxs: Uint32Array,
  allPairsGrid: DG.Grid,
  casesGrid: DG.Grid,
  lines: ILineSeries
} {
  const variates = activities.length;

  let maxAct = 0;
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

  const activityPairsIdxs = new BitArray(allCasesNumber);

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
        //TODO: make more universal mask
        if (k == 0 && diff > 0) {
          if (diff > maxAct)
            maxAct = diff;
          activityPairsIdxs.setBit(pairIdx, true, false);
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
  const pointsFrom = new Uint32Array(activityPairsIdxs.trueCount());
  const pointsTo = new Uint32Array(activityPairsIdxs.trueCount());
  const linesIdxs = new Uint32Array(activityPairsIdxs.trueCount());

  let pairsCounter = 0;
  for (let i = -1; (i = activityPairsIdxs.findNext(i)) !== -1;) {
    pointsFrom[pairsCounter] = molNumFrom[i];
    pointsTo[pairsCounter] = molNumTo[i];
    linesIdxs[pairsCounter] = i;
    pairsCounter++;
  }

  const lines: ILineSeries = {
    from: pointsFrom,
    to: pointsTo,
    drawArrows: true,
    opacity: 0.5,
    color: '60,177,115',
  };

  return {maxAct, diffs, activityMeanNames, linesIdxs, allPairsGrid, casesGrid, lines};
}

function getMmpTrellisPlot(allPairsGrid: DG.Grid, activityMeanNames: Array<string>): DG.Viewer {
  const tp = DG.Viewer.fromType(DG.VIEWER.TRELLIS_PLOT, allPairsGrid.table, {
    xColumnNames: [allPairsGrid.table.columns.byIndex(0).name],
    yColumnNames: [allPairsGrid.table.columns.byIndex(1).name],
    viewerType: 'Summary',
    innerViewerLook: {
      columnNames: activityMeanNames,
      aggregations: [DG.STATS.MED],
      visualizationType: 'bars',
      //colorColumnName: 'age',
      //colorAggrType: 'stdev',
      //invertColorScheme: true,
    },
  });

  return tp;
}

function getMmpScatterPlot(table: DG.DataFrame, activities: DG.Column, maxAct: number) :
[sp: DG.Viewer, sliderInput: DG.InputBase, sliderInputValueDiv: HTMLDivElement] {
  const colX = DG.Column.float('~X', table.rowCount);
  const colY = DG.Column.float('~Y', table.rowCount);
  table.columns.add(colX);
  table.columns.add(colY);
  const sp = DG.Viewer.scatterPlot(table, {
    x: '~X',
    y: '~Y',
    zoomAndFilter: 'no action',
    color: activities.name,
    showXSelector: false,
    showYSelector: false,
  });

  const sliderInput = ui.sliderInput('Cutoff', 0, 0, maxAct);
  const sliderInputValueDiv = ui.divText(sliderInput.stringValue, 'ui-input-description');
  sliderInput.addOptions(sliderInputValueDiv);

  return [sp, sliderInput, sliderInputValueDiv];
}

function drawMolPair(molecules: string[], substr: (ISubstruct | null)[], div: HTMLDivElement, tooltip?: boolean) {
  ui.empty(div);
  const hosts = ui.divH([]);
  if (!tooltip)
    hosts.classList.add('chem-mmpa-context-pane-mol-div');
  const canvasWidth = 150;
  const canvasHeight = 75;
  for (let i = 0; i < 2; i++) {
    const imageHost = ui.canvas(canvasWidth, canvasHeight);
    drawMoleculeToCanvas(0, 0, canvasWidth, canvasHeight, imageHost, molecules[i], '', undefined, substr[i]);
    hosts.append(imageHost);
  }
  div.append(hosts);
};

function moleculesPairInfo(line: number, linesIdxs: Uint32Array, pairsDf: DG.DataFrame, diffs: Array<Float32Array>,
  parentTable: DG.DataFrame, rdkitModule: RDModule, tooltip?: boolean): HTMLDivElement {
  const div = ui.divV([], {style: {width: '100%'}});
  const moleculesDiv = ui.divH([]);
  div.append(moleculesDiv);
  const pairIdx = linesIdxs[line];
  const subsrtFrom = pairsDf.get(MMP_STRUCT_DIFF_FROM_NAME, pairIdx);
  const subsrtTo = pairsDf.get(MMP_STRUCT_DIFF_TO_NAME, pairIdx);
  const moleculeFrom = pairsDf.get(MMP_COLNAME_FROM, pairIdx);
  const moleculeTo = pairsDf.get(MMP_COLNAME_TO, pairIdx);
  if (!tooltip) {
    const fromIdx = pairsDf.get(MMP_COL_PAIRNUM_FROM, pairIdx);
    const toIdx = pairsDf.get(MMP_COL_PAIRNUM_TO, pairIdx);
    const props = getMoleculesPropertiesDiv([fromIdx, toIdx], propColNames, parentTable, propLabelsDiv);
    div.append(props);
  } else {
    //TODO: refine
    const diff = ui.tableFromMap({'Diff': getSigFigs(diffs[0][pairIdx], 4)});
    diff.style.maxWidth = '150px';
    div.append(diff);
  }
  if (subsrtFrom || subsrtTo)
    drawMolPair([moleculeFrom, moleculeTo], [subsrtFrom, subsrtTo], moleculesDiv, tooltip);
  else {
    moleculesDiv.append(ui.divText(`Loading...`));
    getInverseSubstructuresAndAlign([moleculeFrom], [moleculeTo], rdkitModule).then((res) => {
      const {inverse1, inverse2, fromAligned, toAligned} = res;
      pairsDf.set(MMP_STRUCT_DIFF_FROM_NAME, pairIdx, inverse1[0]);
      pairsDf.set(MMP_STRUCT_DIFF_TO_NAME, pairIdx, inverse2[0]);
      pairsDf.set(MMP_COLNAME_FROM, pairIdx, fromAligned[0]);
      pairsDf.set(MMP_COLNAME_TO, pairIdx, toAligned[0]);
      drawMolPair([fromAligned[0], toAligned[0]], [inverse1[0], inverse2[0]], moleculesDiv, tooltip);
    });
  }
  return div;
};

function getMoleculesPropertiesDiv(idxs: number[], propertiesColumnsNames: string[], parentTable: DG.DataFrame,
  labelsDiv: HTMLDivElement): HTMLDivElement {
  const grid = grok.shell.tv.grid;
  const properties = [
    ui.divV([], 'chem-mmpa-context-pane-properties-div'),
    ui.divV([], 'chem-mmpa-context-pane-properties-div')];
  for (const col of propertiesColumnsNames) {
    const colorCoding = parentTable.col(col)!.tags[DG.TAGS.COLOR_CODING_TYPE];
    let colorCoded = false;
    if (colorCoding && colorCoding !== DG.COLOR_CODING_TYPE.OFF)
      colorCoded = true;
    for (let i = 0; i < 2; i++) {
      const value = ui.divText(`${parentTable.get(col, idxs[i])}`, 'chem-mmpa-prop-value');
      if (colorCoded) {
        const color = DG.Color.toHtml(grid.cell(col, idxs[i]).color);
        if (grid.col(col)?.isTextColorCoded)
          value.style.color = color;
        else {
          value.style.backgroundColor = color,
          value.style.color = pickTextColorBasedOnBgColor(color, '#FFFFFF', '#000000');
        }
      }
      properties[i].append(value);
    }
  }
  return ui.divH([
    labelsDiv,
    properties[0],
    properties[1],
  ], {style: {width: '100%'}});
}

function getPropsColNamesDiv() {
  const labelsDiv = ui.divV([], {style: {gap: '5px', width: '100px'}});
  for (const col of propColNames) {
    const label = ui.divText(`${col}`, 'chem-similarity-prop-label');
    ui.tooltip.bind(label, col);
    labelsDiv.append(label);
  }
  return labelsDiv;
}

function runMmpChemSpace(table: DG.DataFrame, molecules: DG.Column, sp: DG.Viewer, lines: ILineSeries,
  linesIdxs: Uint32Array, pairsDf: DG.DataFrame, diffs: Array<Float32Array>,
  rdkitModule: RDModule): ScatterPlotLinesRenderer {
  const chemSpaceParams = {
    seqCol: molecules,
    methodName: DimReductionMethods.UMAP,
    similarityMetric: BitArrayMetricsNames.Tanimoto as BitArrayMetrics,
    embedAxesNames: [MMP_COLNAME_CHEMSPACE_X, MMP_COLNAME_CHEMSPACE_Y],
    options: {},
  };

  //@ts-ignore
  const spEditor = new ScatterPlotLinesRenderer(sp as DG.ScatterPlotViewer,
    MMP_COLNAME_CHEMSPACE_X, MMP_COLNAME_CHEMSPACE_Y,
    lines, ScatterPlotCurrentLineStyle.bold);


  spEditor.lineClicked.subscribe((event: MouseOverLineEvent) => {
    spEditor.currentLineId = event.id;
    if (event.id !== -1) {
      grok.shell.o = moleculesPairInfo(event.id, linesIdxs, pairsDf, diffs, table, rdkitModule);
      lastSelectedPair = event.id;
    }
  });

  spEditor.lineHover.pipe(debounceTime(500)).subscribe((event: MouseOverLineEvent) => {
    ui.tooltip.show(
      moleculesPairInfo(event.id, linesIdxs, pairsDf, diffs, table, rdkitModule, true), event.x, event.y);
  });

  const progressBarSpace = DG.TaskBarProgressIndicator.create(`Running Chemical space...`);
  chemSpace(chemSpaceParams).then((res) => {
    const embeddings = res.coordinates;
    for (const col of embeddings)
      table.columns.replace(col.name, col);
    progressBarSpace.close();
  });

  return spEditor;
}

async function getMmpMcs(molecules1: string[], molecules2: string[]): Promise<string[]> {
  const service = await getRdKitService();
  const molecules: [string, string][] = new Array<[string, string]>(molecules1.length);
  for (let i = 0; i <molecules1.length; i++)
    molecules[i] = [molecules1[i], molecules2[i]];

  return await service.mmpGetMcs(molecules);
}

async function getInverseSubstructuresAndAlign(from: string[], to: string[], module: RDModule):
  Promise<{
    inverse1: (ISubstruct | null)[],
    inverse2: (ISubstruct | null)[],
    fromAligned: string[],
    toAligned: string[]}> {
  const fromAligned = new Array<string>(from.length);
  const toAligned = new Array<string>(from.length);
  const res1 = new Array<(ISubstruct | null)>(from.length);
  const res2 = new Array<(ISubstruct | null)>(from.length);

  const mcs = await getMmpMcs(from, to);

  for (let i = 0; i < from.length; i++) {
    //aligning molecules
    const mcsMol = module.get_qmol(mcs[i]);
    mcsMol.set_new_coords();
    //revise
    const mol1 = module.get_mol(from[i]);
    const mol2 = module.get_mol(to[i]);

    mol1.generate_aligned_coords(mcsMol, JSON.stringify({
      useCoordGen: true,
      allowRGroups: true,
      acceptFailure: false,
      alignOnly: true,
    }));

    mol2.generate_aligned_coords(mcsMol, JSON.stringify({
      useCoordGen: true,
      allowRGroups: true,
      acceptFailure: false,
      alignOnly: true,
    }));

    fromAligned[i] = mol1.get_molblock();
    toAligned[i] = mol2.get_molblock();

    //highlight fragment
    // const matches = mol1.get_substruct_matches(mcsMol);
    // for (let i = 0; i < matches.length; i++) {
    // }

    mol1?.delete();
    mol2?.delete();

    const substruct1 = getUncommonAtomsAndBonds(from[i], mcsMol, module, '#bc131f');
    const substruct2 = getUncommonAtomsAndBonds(to[i], mcsMol, module, '#49bead');

    res1[i] = substruct1;
    res2[i] = substruct2;

    mcsMol?.delete();
  }
  return {inverse1: res1, inverse2: res2, fromAligned: fromAligned, toAligned: toAligned};
}

function drawMoleculeLabels(sp: DG.ScatterPlotViewer, table: DG.DataFrame, molCol: DG.Column) {
  sp.onAfterDrawScene.subscribe(() => {
    const maxPointsOnScreen = 20;
    const rowCount = table.rowCount;
    const pointsOnScreen = new Array<DG.Point | null>(maxPointsOnScreen).fill(null);
    const pointsOnScreenIdxs = new Uint32Array(maxPointsOnScreen);
    const {maxX, maxY, minX, minY} = sp.viewBox;
    const maxYWithAxisBox = maxY - sp.xAxisBox.height;
    let counter = 0;
    for (let i = 0; i < rowCount; i++) {
      //@ts-ignore
      const point = sp.pointToScreen(i);
      //check whether point is within a viewport (on the screen)
      if (point.x > minX && point.x < maxX && point.y > minY && point.y < maxYWithAxisBox && table.filter.get(i)) {
        if (counter == maxPointsOnScreen)
          return; // returning if there are more than maximum allowed points in total on the screen
        else {
          pointsOnScreen[counter] = point;
          pointsOnScreenIdxs[counter] = i;
          counter++;
        }
      }
    }

    const minDistance = 50;
    const N = counter * (counter + 1) / 2;
    const pointsToDrawIdxs = new BitArray(counter);
    const distancesMatrix = new Float32Array(N);

    //for each point check that distance to all other points is not less than minimum allowed - in this case draw label
    for (let i = 0; i < counter; i++) {
      let lessThenMinDist = false;
      for (let j = 0; j < counter; j++) {
        if (i === j)
          continue;
        else {
          const dist = j < i ? distancesMatrix[(j * counter - (j - 1) * j / 2) + (i - j) - 1] :
            pointsOnScreen[i]!.distanceTo(pointsOnScreen[j]!);
          if (j > i)
            distancesMatrix[(i * counter - (i - 1) * i / 2) + (j - i) - 1] = dist;
          if (dist < minDistance && pointsOnScreen[i]!.x < pointsOnScreen[j]!.x &&
              pointsOnScreen[i]!.y < pointsOnScreen[j]!.y)
            lessThenMinDist = true;
        }
      }
      if (!lessThenMinDist)
        pointsToDrawIdxs.setBit(i, true, false);
    }

    const ctxMain = sp.getInfo()['canvas'].getContext('2d') as CanvasRenderingContext2D;
    ctxMain.beginPath();
    ctxMain.strokeStyle = `rgba(0, 0, 0, 0.1)`;

    for (let i = -1; (i = pointsToDrawIdxs.findNext(i)) !== -1;) {
      const imageHost = ui.canvas(100, 50);
      let molecule = molCol.get(pointsOnScreenIdxs[i]);
      if (molCol.tags[DG.TAGS.UNITS] === DG.chem.Notation.Smiles) {
        //convert to molFile to draw in coordinates similar to dataframe cell
        molecule = convertMolNotation(molecule, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
      }
      drawMoleculeToCanvas(0, 0, 100, 50, imageHost, molecule);
      ctxMain.drawImage(imageHost, pointsOnScreen[i]!.x, pointsOnScreen[i]!.y, 100, 50);
      ctxMain.strokeRect(pointsOnScreen[i]!.x, pointsOnScreen[i]!.y, 100, 50);
    }
    ctxMain.closePath();
  });
}

export class MmpAnalysis {
  parentTable: DG.DataFrame;
  parentCol: DG.Column;
  mmpRules: MmpRules;
  mmpView: DG.View;
  enableFilters: boolean = true;
  //transformations tab objects
  allPairsGrid: DG.Grid;
  transformationsMask: DG.BitSet;
  casesGrid: DG.Grid;
  pairsMask: DG.BitSet;
  //cliffs tab objects
  diffs: Array<Float32Array>;
  linesIdxs: Uint32Array;
  cutoffMask: DG.BitSet;
  linesMask: BitArray;
  linesRenderer: ScatterPlotLinesRenderer | null = null;
  lines: ILineSeries;
  //rdkit
  rdkitModule: RDModule;

  constructor(table: DG.DataFrame, molecules: DG.Column, rules: MmpRules, diffs: Array<Float32Array>,
    linesIdxs: Uint32Array, allPairsGrid: DG.Grid, casesGrid: DG.Grid, tp: DG.Viewer, sp: DG.Viewer,
    sliderInput: DG.InputBase, sliderInputValueDiv: HTMLDivElement,
    linesEditor: ScatterPlotLinesRenderer, lines: ILineSeries, rdkitModule: RDModule) {
    this.rdkitModule = rdkitModule;

    this.parentTable = table;
    this.parentCol = molecules;
    this.mmpRules = rules;

    this.diffs = diffs;
    this.allPairsGrid = allPairsGrid;
    this.casesGrid = casesGrid;
    this.lines = lines;
    this.linesIdxs = linesIdxs;

    //transformations tab
    this.transformationsMask = DG.BitSet.create(this.allPairsGrid.dataFrame.rowCount);
    this.transformationsMask.setAll(true);
    this.parentTable.onCurrentRowChanged.pipe(debounceTime(1000)).subscribe(() => {
      if (this.parentTable.currentRowIdx !== -1) {
        this.refilterAllPairs(true);
        this.refreshPair(this.rdkitModule);
      }
    });

    this.refilterAllPairs(true);

    this.allPairsGrid.table.onCurrentRowChanged.subscribe(() => {
      if (this.allPairsGrid.table.currentRowIdx !== -1)
        this.refreshPair(this.rdkitModule);
    });

    this.mmpView = DG.View.create();
    this.mmpView.name = MMP_VIEW_NAME;
    this.mmpView.box = true;

    this.pairsMask = DG.BitSet.create(this.casesGrid.dataFrame.rowCount);
    this.pairsMask.setAll(false);
    this.refreshPair(this.rdkitModule);

    //Cliffs tab
    sliderInput.onChanged(() => {
      sliderInputValueDiv.innerText = sliderInput.value === 0 ? '0' : getSigFigs(sliderInput.value, 4).toString();
      this.refilterCliffs(sliderInput.value, true);
    });

    ui.tooltip.bind(sliderInput.captionLabel, 'Select the cutoff by activity difference');
    ui.tooltip.bind(sliderInput.input, 'Activity value cutoff');
    sp.root.style.width = '100%';
    const cliffs = ui.divV([sliderInput.root, ui.box(sp.root, {style: {maxHeight: '100px', paddingRight: '6px'}})]);

    this.cutoffMask = DG.BitSet.create(this.parentTable.rowCount);
    this.cutoffMask.setAll(true);

    this.linesMask = new BitArray(this.casesGrid.dataFrame.rowCount);
    this.linesMask.setAll(true, false);

    this.linesRenderer = linesEditor;

    //tabs
    const tabs = ui.tabControl(null, false);

    tabs.addPane(MMP_TAB_TRANSFORMATIONS, () => {
      const fragmentsHeader = ui.h1('Fragments');
      $(fragmentsHeader).css('padding', '4px');
      $(fragmentsHeader).css('margin', '0px');
      this.allPairsGrid.root.prepend(fragmentsHeader);
      const fragments = ui.splitV([ui.box(fragmentsHeader, {style: {maxHeight: '30px'}}), this.allPairsGrid.root]);

      const casesHeader = ui.h1('Pairs');
      $(casesHeader).css('padding', '4px');
      $(casesHeader).css('margin', '0px');

      this.casesGrid.root.prepend(casesHeader);
      const cases = ui.splitV([ui.box(casesHeader, {style: {maxHeight: '30px'}}), ui.box(this.casesGrid.root)]);

      return ui.splitV([fragments, cases], {}, true);
    });
    tabs.addPane(MMP_TAB_FRAGMENTS, () => {
      return tp.root;
    });
    tabs.addPane(MMP_TAB_CLIFFS, () => {
      return cliffs;
    });

    tabs.onTabChanged.subscribe(() => {
      currentTab = tabs.currentPane.name;
      if (tabs.currentPane.name == MMP_TAB_TRANSFORMATIONS) {
        this.enableFilters = true;
        this.refilterAllPairs(false);
      } else if (tabs.currentPane.name == MMP_TAB_FRAGMENTS) {
        this.refreshFilterAllPairs();
        this.enableFilters = false;
      } else if (tabs.currentPane.name == MMP_TAB_CLIFFS) {
        this.refilterCliffs(sliderInput.value, false);
        if (lastSelectedPair) {
          grok.shell.o = moleculesPairInfo(lastSelectedPair, this.linesIdxs, this.casesGrid.dataFrame,
            this.diffs, this.parentTable, this.rdkitModule);
        }
      }
    });

    const decript1 = 'Shows all fragmental substitutions for a given molecule';
    const decript2 = 'Analysis of fragments versus explored value';
    const decript3 = 'Cliffs analysis';

    tabs.getPane(MMP_TAB_TRANSFORMATIONS).header.onmouseover =
      (ev): void => ui.tooltip.show(decript1, ev.clientX, ev.clientY + 5);
    tabs.getPane(MMP_TAB_FRAGMENTS).header.onmouseover =
      (ev): void => ui.tooltip.show(decript2, ev.clientX, ev.clientY + 5);
    tabs.getPane(MMP_TAB_CLIFFS).header.onmouseover =
      (ev): void => ui.tooltip.show(decript3, ev.clientX, ev.clientY + 5);

    this.mmpView.append(tabs);
  }

  static async init(table: DG.DataFrame, molecules: DG.Column, activities: DG.ColumnList) {
    //rdkit module
    const module = getRdKitModule();

    //initial calculations
    const frags = await getMmpFrags(molecules);
    const [mmpRules, allCasesNumber] = getMmpRules(frags);

    //Transformations tab
    const {maxAct, diffs, activityMeanNames, linesIdxs, allPairsGrid, casesGrid, lines} =
      await getMmpActivityPairsAndTransforms(molecules, activities, mmpRules, allCasesNumber);

    //Fragments tab
    const tp = getMmpTrellisPlot(allPairsGrid, activityMeanNames);

    //Cliffs tab
    const [sp, sliderInput, sliderInputValueDiv] = getMmpScatterPlot(table, activities.byIndex(0), maxAct);
    drawMoleculeLabels(sp as DG.ScatterPlotViewer, table, molecules);

    //running internal chemspace
    const linesEditor = runMmpChemSpace(table, molecules, sp, lines, linesIdxs, casesGrid.dataFrame, diffs, module);

    propColNames = table.columns.names()
      .filter((name) => name !== molecules.name &&
        name !== MMP_COLNAME_CHEMSPACE_X &&
        name !== MMP_COLNAME_CHEMSPACE_Y);

    propLabelsDiv = getPropsColNamesDiv();

    DG.debounce((table.onMetadataChanged), 50)
      .subscribe(async (_: any) => {
        if (lastSelectedPair && currentTab === MMP_TAB_CLIFFS)
          grok.shell.o = moleculesPairInfo(lastSelectedPair, linesIdxs, casesGrid.dataFrame, diffs, table, module);
      });

    return new MmpAnalysis(table, molecules, mmpRules, diffs, linesIdxs, allPairsGrid, casesGrid,
      tp, sp, sliderInput, sliderInputValueDiv, linesEditor, lines, module);
  }

  async refreshPair(rdkitModule: RDModule) {
    const progressBarPairs = DG.TaskBarProgressIndicator.create(`Refreshing pairs...`);
    const idxParent = this.parentTable.currentRowIdx;
    let idxPairs = -1;

    const idx = this.allPairsGrid.table.currentRowIdx;

    const ruleSmi1 = this.allPairsGrid.table.getCol('From').get(idx);
    const ruleSmi2 = this.allPairsGrid.table.getCol('To').get(idx);

    //search for numbers in smiles rules
    const ruleSmiNum1 = this.mmpRules.smilesFrags.indexOf(ruleSmi1);
    const ruleSmiNum2 = this.mmpRules.smilesFrags.indexOf(ruleSmi2);

    this.pairsMask.setAll(false);
    let counter = 0;
    const cases: number[] = [];
    const diffFromSubstrCol = this.casesGrid.dataFrame.getCol(MMP_STRUCT_DIFF_FROM_NAME);
    const diffToSubstrCol = this.casesGrid.dataFrame.getCol(MMP_STRUCT_DIFF_TO_NAME);
    const diffFrom = this.casesGrid.dataFrame.getCol(MMP_COLNAME_FROM);
    const diffTo = this.casesGrid.dataFrame.getCol(MMP_COLNAME_TO);

    //search specific rule

    for (let i = 0; i < this.mmpRules.rules.length; i++) {
      const first = this.mmpRules.rules[i].smilesRule1; //.pairs[j].firstStructure;
      const second = this.mmpRules.rules[i].smilesRule2; //.pairs[j].secondStructure;
      for (let j = 0; j < this.mmpRules.rules[i].pairs.length; j++) {
        if (ruleSmiNum1 == first && ruleSmiNum2 == second) {
          this.pairsMask.set(counter, true, false);
          if (diffFromSubstrCol.get(counter) === null)
            cases.push(counter);
          if (this.mmpRules.rules[i].pairs[j].firstStructure == idxParent)
            idxPairs = counter;
        }
        counter++;
      }
    }

    //recover uncalculated highlights
    const pairsFrom = Array<string>(cases.length);
    const pairsTo = Array<string>(cases.length);
    for (let i = 0; i < cases.length; i++) {
      pairsFrom[i] = diffFrom.get(cases[i]);
      pairsTo[i] = diffTo.get(cases[i]);
    }
    const {inverse1, inverse2, fromAligned, toAligned} =
      await getInverseSubstructuresAndAlign(pairsFrom, pairsTo, rdkitModule);
    for (let i = 0; i < cases.length; i++) {
      diffFrom.set(cases[i], fromAligned[i]);
      diffTo.set(cases[i], toAligned[i]);
      diffFromSubstrCol.set(cases[i], inverse1[i]);
      diffToSubstrCol.set(cases[i], inverse2[i]);
    }

    this.casesGrid.dataFrame.filter.copyFrom(this.pairsMask);

    this.casesGrid.setOptions({
      pinnedRowColumnNames: [],
      pinnedRowValues: [],
    });
    if (idxPairs >= 0) {
      this.casesGrid.setOptions({
        pinnedRowValues: [idxPairs.toString()],
        pinnedRowColumnNames: [MMP_COL_PAIRNUM],
      });
    }

    progressBarPairs.close();
  }

  refreshFilterAllPairs() {
    const consistsBitSet: DG.BitSet = DG.BitSet.create(this.allPairsGrid.dataFrame.rowCount);
    consistsBitSet.setAll(true);
    this.allPairsGrid.dataFrame.filter.copyFrom(consistsBitSet);
  }

  refilterAllPairs(rowChanged: boolean) {
    let idxTrue = -1;
    if (rowChanged) {
      const idx = this.parentTable.currentRowIdx;
      this.transformationsMask.setAll(false);

      for (let i = 0; i < this.mmpRules.rules.length; i++) {
        for (let j = 0; j < this.mmpRules.rules[i].pairs.length; j++) {
          const fs = this.mmpRules.rules[i].pairs[j].firstStructure;
          if (idx == fs) {
            if (idxTrue == -1)
              idxTrue = i;
            this.transformationsMask.set(i, true, false);
            break;
          }
        }
      }
    }

    if (this.enableFilters) {
      this.allPairsGrid.dataFrame.filter.copyFrom(this.transformationsMask);
      this.allPairsGrid.invalidate();
      if (rowChanged)
        this.allPairsGrid.table.currentRowIdx = idxTrue;
    }
  }

  refilterCliffs(cutoff: number, refilter: boolean) {
    if (refilter) {
      this.linesMask.setAll(false, false);
      this.cutoffMask.setAll(false);
      for (let i = 0; i < this.lines.from.length; i++) {
        //TODO: refine
        if (this.diffs[0][this.linesIdxs[i]] >= cutoff) {
          this.cutoffMask.set(this.lines.from[i], true, false);
          this.cutoffMask.set(this.lines.to[i], true, false);
          this.linesMask.setBit(i, true, false);
        }
      }
    }
    this.parentTable.filter.copyFrom(this.cutoffMask);
    this.linesRenderer!.linesVisibility = this.linesMask;
  }
}
