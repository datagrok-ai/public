import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule, getRdKitService, getUncommonAtomsAndBonds} from '../utils/chem-common-rdkit';

import {DimReductionMethods} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {BitArrayMetrics, BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {chemSpace} from './chem-space';
import $ from 'cash-dom';
import {ISubstruct} from '../rendering/rdkit-cell-renderer';
import {SUBSTRUCT_COL} from '../constants';
import { getMCS } from '../utils/most-common-subs';

const MMP_COLNAME_FROM = 'From';
const MMP_COLNAME_TO = 'To';
const MMP_COLNAME_PAIRS = 'Pairs';
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

interface ILine {
  id: number;
  mols: number[];
  a: number[]; // [x, y]
  b: number[]; // [x, y]
}

interface IRenderedLines {
  lines: ILine[];
}

function renderLines(sp: DG.ScatterPlotViewer, xAxis: string, yAxis: string, linesRes: IRenderedLines): ILine [] {
  const lines = linesRes.lines;
  const canvas = (sp.getInfo() as {[index: string] : any})['canvas'];
  const ctx = canvas.getContext('2d') as CanvasRenderingContext2D;
  const x = sp.dataFrame!.columns.byName(xAxis);
  const y = sp.dataFrame!.columns.byName(yAxis);
  for (let i = 0; i < lines.length; i++) {
    const pointFrom = sp.worldToScreen(x.get(lines[i].mols[0]), y.get(lines[i].mols[0]));
    const pointTo = sp.worldToScreen(x.get(lines[i].mols[1]), y.get(lines[i].mols[1]));
    lines[i].a = [pointFrom.x, pointFrom.y];
    lines[i].b = [pointTo.x, pointTo.y];
    const line = new Path2D();
    line.moveTo(lines[i].a[0], lines[i].a[1]);
    const color = '0,128,0';
    const opacity = 0.5;
    ctx.strokeStyle = `rgba(${color},${opacity})`;
    ctx.lineWidth = 2;
    line.lineTo(lines[i].b[0], lines[i].b[1]);
    ctx.stroke(line);
  }
  return lines;
}

type MmpRules = {
  rules: {
    smilesRule1: number,
    smilesRule2: number,
    pairs: {firstStructure: number, secondStructure: number}[]
  } [],
  smilesFrags: string[]
};

async function getMmpFrags(molecules: DG.Column): Promise<[string, string][][]> {
  const service = await getRdKitService();
  const res = await service.mmpGetFragments(molecules.toList());
  return res;
}

function getMmpRules(frags: [string, string][][]): MmpRules {
  const mmpRules: MmpRules = {rules: [], smilesFrags: []};
  const dim = frags.length;
  let ruleCounter = 0;

  for (let i = 0; i < dim; i++) {
    const dim1 = frags[i].length;
    for (let j = i + 1; j < dim; j++) {
      const dim2 = frags[j].length;
      let core = '';
      let r1 = '';
      let r2 = '';

      //here we get the best possibple fragment pair
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

      if (core === '')
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
      } else {
        mmpRules.rules[ruleIndexStraight].pairs
          .push({firstStructure: indxFirst ? i : j, secondStructure: indxFirst ? j : i});
        mmpRules.rules[ruleIndexInverse].pairs
          .push({firstStructure: indxFirst ? j : i, secondStructure: indxFirst ? i : j});
      }
    }
  }

  return mmpRules;
}

function getMmpActivityPairs(activities: DG.Column, mmpRules: MmpRules):
  [maxAct: number,
    activityPairs:{first: number, second: number, diff: number}[],
    allPairsGrid: DG.Grid
  ] {
  let maxAct = 0;
  const from = new Array<string>(mmpRules.rules.length);
  const to = new Array<string>(mmpRules.rules.length);
  const occasions = new Int32Array(mmpRules.rules.length);
  const meanDiff = new Float32Array(mmpRules.rules.length);

  //set activity differences
  const activityPairs: {first: number, second: number, diff: number}[] = [];
  for (let i = 0; i < mmpRules.rules.length; i++) {
    from[i] = mmpRules.smilesFrags[mmpRules.rules[i].smilesRule1];
    to[i] = mmpRules.smilesFrags[mmpRules.rules[i].smilesRule2];
    occasions[i] = mmpRules.rules[i].pairs.length;

    let mean: number = 0;
    for (let j = 0; j < occasions[i]; j++) {
      const idx1 = mmpRules.rules[i].pairs[j].firstStructure;
      const idx2 = mmpRules.rules[i].pairs[j].secondStructure;
      const val1 = activities.get(idx1);
      const val2 = activities.get(idx2);
      const diff = val2 - val1;
      if (diff > 0) {
        if (diff > maxAct)
          maxAct = diff;

        activityPairs.push({first: idx1, second: idx2, diff: diff});
      }

      mean += diff;
    }
    mean/= occasions[i];
    meanDiff[i] = mean;
  }

  const fromCol = DG.Column.fromList('string', MMP_COLNAME_FROM, from);
  const toCol = DG.Column.fromList('string', MMP_COLNAME_TO, to);
  const occasionsCol = DG.Column.fromInt32Array(MMP_COLNAME_PAIRS, occasions);
  const meanDiffCol = DG.Column.fromFloat32Array(MMP_COLNAME_MEANDIFF, meanDiff);

  fromCol.semType = DG.SEMTYPE.MOLECULE;
  toCol.semType = DG.SEMTYPE.MOLECULE;
  occasionsCol.semType = DG.TYPE.INT;

  const dfAllPairs = DG.DataFrame.fromColumns([fromCol, toCol, occasionsCol, meanDiffCol]);
  const allPairsGrid = dfAllPairs.plot.grid();

  return [maxAct, activityPairs, allPairsGrid];
}

function getMmpPairedTransforms(molecules: DG.Column, activities: DG.Column, mmpRules: MmpRules): DG.Grid {
  const pairsFrom: string[] = [];
  const pairsTo: string[] = [];
  const singleDiff: number[] = [];

  for (let i = 0; i < mmpRules.rules.length; i++) {
    for (let j = 0; j < mmpRules.rules[i].pairs.length; j++) {
      const idx1 = mmpRules.rules[i].pairs[j].firstStructure;
      const idx2 = mmpRules.rules[i].pairs[j].secondStructure;
      pairsFrom.push(molecules.get(idx1));
      pairsTo.push(molecules.get(idx2));
      singleDiff.push(activities.get(idx2) - activities.get(idx1));
    }
  }

  const pairsFromCol = DG.Column.fromStrings(MMP_COLNAME_FROM, pairsFrom);
  const pairsToCol = DG.Column.fromStrings(MMP_COLNAME_TO, pairsTo);
  const singleDiffCol = DG.Column.fromList('double', MMP_COLNAME_DIFF, singleDiff);
  const structureDiffFromCol = DG.Column.fromType('object', MMP_STRUCT_DIFF_FROM_NAME, pairsFrom.length);
  const structureDiffToCol = DG.Column.fromType('object', MMP_STRUCT_DIFF_TO_NAME, pairsFrom.length);

  pairsFromCol.semType = DG.SEMTYPE.MOLECULE;
  pairsToCol.semType = DG.SEMTYPE.MOLECULE;
  pairsFromCol.temp[SUBSTRUCT_COL] = MMP_STRUCT_DIFF_FROM_NAME;
  pairsToCol.temp[SUBSTRUCT_COL] = MMP_STRUCT_DIFF_TO_NAME;
  const pairedTransformations = DG.DataFrame.
    fromColumns([pairsFromCol, pairsToCol, singleDiffCol, structureDiffFromCol, structureDiffToCol]);

  return pairedTransformations.plot.grid();
}

function getMmpTrellisPlot(allPairsGrid: DG.Grid): DG.Viewer {
  return DG.Viewer.fromType(DG.VIEWER.TRELLIS_PLOT, allPairsGrid.table, {
    xColumnNames: [allPairsGrid.table.columns.byIndex(0).name],
    yColumnNames: [allPairsGrid.table.columns.byIndex(1).name],
  });
}

function getMmpScatterPlot(table: DG.DataFrame, activities: DG.Column, maxAct: number) :
[sp: DG.Viewer, sliderInput: DG.InputBase] {
  const colX = DG.Column.float('~X', table.rowCount);
  const colY = DG.Column.float('~Y', table.rowCount);
  table.columns.add(colX);
  table.columns.add(colY);
  const sp = DG.Viewer.scatterPlot(table, {
    x: '~X',
    y: '~Y',
    color: activities.name,
  });

  const property =
  {
    'name': 'Cutoff',
    'type': DG.TYPE.FLOAT,
    'showSlider': true,
    'min': 0,
    'max': maxAct,
    'nullable': false,
  };

  const slider = DG.Property.fromOptions(property);
  const initialCutOff = {Cutoff: 0};
  const sliderInput = ui.input.forProperty(slider, initialCutOff);

  return [sp, sliderInput];
}

function mmpRunChemSpace(table: DG.DataFrame, molecules: DG.Column, mmpRules: MmpRules, sp: DG.Viewer): void {
  const chemSpaceParams = {
    seqCol: molecules,
    methodName: DimReductionMethods.UMAP,
    similarityMetric: BitArrayMetricsNames.Tanimoto as BitArrayMetrics,
    embedAxesNames: [MMP_COLNAME_CHEMSPACE_X, MMP_COLNAME_CHEMSPACE_Y],
    options: {},
  };

  const progressBarSpace = DG.TaskBarProgressIndicator.create(`Running Chemical space...`);
  chemSpace(chemSpaceParams).then((res) => {
    const embeddings = res.coordinates;
    for (const col of embeddings)
      table.columns.replace(col.name, col);

    const lines: {id: number, mols: number[], a: number[], b: number[]}[] = [];
    for (let i = 0; i < mmpRules.rules.length; i++) {
      for (let j = 0; j < mmpRules.rules[i].pairs.length; j++) {
        const fs = mmpRules.rules[i].pairs[j].firstStructure;
        const ss = mmpRules.rules[i].pairs[j].secondStructure;
        lines.push({id: i, mols: [fs, ss], a: [], b: []});
      }
    }

    const linesRes = {lines: lines};

    sp!.onEvent('d4-before-draw-scene')
      .subscribe((_: any) => {
        renderLines(sp as DG.ScatterPlotViewer, MMP_COLNAME_CHEMSPACE_X, MMP_COLNAME_CHEMSPACE_Y, linesRes);
      });

    progressBarSpace.close();
  });
}

async function getMmpMcs(molecules1: string[], molecules2: string[]): Promise<string[]> {
  const service = await getRdKitService();
  const molecules: [string, string][] = new Array<[string, string]>(molecules1.length);
  for (let i = 0; i <molecules1.length; i++)
    molecules[i] = [molecules1[i], molecules2[i]];

  return await service.mmpGetMcs(molecules);
}

async function getInverseSubstructures(from: string[], to: string[]):
  Promise<[(ISubstruct | null)[], (ISubstruct | null)[]]> {
  //TODO: the length is known
  const res1: (ISubstruct | null)[] = [];
  const res2: (ISubstruct | null)[] = [];

  const mcs = await getMmpMcs(from, to);

  //const res: string[] = [];
  for (let i = 0; i < from.length; i++) {
    const module = getRdKitModule();
    const mcsMol = module.get_qmol(mcs[i]);

    const coll = DG.Column.fromStrings('ss', [from[i], to[i]]);
    const mcs2 = getMCS(coll, true, true);

    const substruct1 = getUncommonAtomsAndBonds(from[i], mcsMol, module, '#bc131f');
    const substruct2 = getUncommonAtomsAndBonds(to[i], mcsMol, module, '#49bead');

    // const mol1 = module.get_mol(from[i]);
    // const mol2 = module.get_mol(to[i]);

    // mol2.generate_aligned_coords(mol1, JSON.stringify({
    //   useCoordGen: true,
    //   allowRGroups: true,
    //   acceptFailure: false,
    //   alignOnly: true,
    // }));

    // drawMoleculeToCanvas(0, 0, 200, 100, this.canvasMol1, mol1.get_molblock(), '',
    //   {normalizeDepiction: true, straightenDepiction: true}, substruct1);
    // drawMoleculeToCanvas(0, 0, 200, 100, this.canvasMol2, mol2.get_molblock(), '',
    //   {normalizeDepiction: true, straightenDepiction: true}, substruct2);
    res1.push(substruct1);
    res2.push(substruct2);
    //res3.push(mcs[i]);
    //res.push(mcs === null ? '' : mcs);
    mcsMol?.delete();
    // mol1.delete;
    // mol2.delete;
  }
  return [res1, res2];//res;
}

export class MmpAnalysis {
  parentTable: DG.DataFrame;
  parentCol: DG.Column;
  mmpRules: MmpRules;
  mmpView: DG.View;
  enableFilters: boolean = false;
  //transformations tab objects
  allPairsGrid: DG.Grid;
  transformationsMask: DG.BitSet;
  casesGrid: DG.Grid;
  pairsMask: DG.BitSet;
  //cliffs tab objects
  activityPairs: {first: number, second: number, diff: number}[];
  cutoffMask: DG.BitSet;

  constructor(table: DG.DataFrame, molecules: DG.Column, rules: MmpRules,
    activityPairs: {first: number, second: number, diff: number}[], allPairsGrid: DG.Grid, casesGrid: DG.Grid,
    tp: DG.Viewer, sp: DG.Viewer, sliderInput: DG.InputBase) {
    this.parentTable = table;
    this.parentCol = molecules;
    this.mmpRules = rules;

    this.activityPairs = activityPairs;
    this.allPairsGrid = allPairsGrid;
    this.casesGrid = casesGrid;

    //transformations tab
    this.transformationsMask = DG.BitSet.create(this.allPairsGrid.dataFrame.rowCount);
    this.transformationsMask.setAll(true);
    this.parentTable.onCurrentRowChanged.subscribe(() => {
      this.refilterAllPairs(true);
      this.refreshPair();
    });

    this.refilterAllPairs(true);

    this.allPairsGrid.table.currentRowIdx = 0;
    this.allPairsGrid.table.onCurrentRowChanged.subscribe(() => {
      this.refreshPair();
    });

    this.mmpView = DG.View.create();
    this.mmpView.name = MMP_VIEW_NAME;
    this.mmpView.box = true;

    this.pairsMask = DG.BitSet.create(this.casesGrid.dataFrame.rowCount);
    this.pairsMask.setAll(false);
    this.refreshPair();

    //Cliffs tab
    sliderInput.onChanged(() => {
      this.refilterCliffs(sliderInput.value, true);
    });

    sp.root.prepend(sliderInput.root);
    const cliffs = ui.splitV([ui.box(sliderInput.root, {style: {maxHeight: '100px'}}), sp.root]);

    this.cutoffMask = DG.BitSet.create(this.parentTable.rowCount);
    this.cutoffMask.setAll(true);

    //tabs
    const tabs = ui.tabControl(null, false);

    tabs.addPane(MMP_TAB_TRANSFORMATIONS, () => {
      const fragmentsHeader = ui.h1('Fragments');
      $(fragmentsHeader).css('padding', '6px');
      $(fragmentsHeader).css('margin', '0px');
      this.allPairsGrid.root.prepend(fragmentsHeader);
      const fragments = ui.splitV([ui.box(fragmentsHeader, {style: {maxHeight: '30px'}}), this.allPairsGrid.root]);

      const casesHeader = ui.h1('Pairs');
      $(casesHeader).css('padding', '6px');
      $(casesHeader).css('margin', '0px');

      this.casesGrid.root.prepend(casesHeader);
      const cases = ui.splitV([ui.box(casesHeader, {style: {maxHeight: '30px'}}), ui.box(this.casesGrid.root)]);

      return ui.splitV([fragments, cases]);
    });
    tabs.addPane(MMP_TAB_FRAGMENTS, () => {
      return tp.root;
    });
    tabs.addPane(MMP_TAB_CLIFFS, () => {
      return cliffs;
    });

    tabs.onTabChanged.subscribe(() => {
      if (tabs.currentPane.name == MMP_TAB_TRANSFORMATIONS) {
        this.enableFilters = true;
        this.refilterAllPairs(false);
      } else if (tabs.currentPane.name == MMP_TAB_FRAGMENTS) {
        this.refreshFilterAllPairs();
        this.enableFilters = false;
      } else if (tabs.currentPane.name == MMP_TAB_CLIFFS)
        this.refilterCliffs(sliderInput.value, false);
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

  static async init(table: DG.DataFrame, molecules: DG.Column, activities:DG.Column) {
    //initial calculations
    const frags = await getMmpFrags(molecules);
    const mmpRules = getMmpRules(frags);

    //Transformations tab
    const [maxAct, activityPairs, allPairsGrid] = getMmpActivityPairs(activities, mmpRules);
    const casesGrid = getMmpPairedTransforms(molecules, activities, mmpRules);

    //Fragments tab
    const tp = getMmpTrellisPlot(allPairsGrid);

    //Cliffs tab
    const [sp, sliderInput] = getMmpScatterPlot(table, activities, maxAct);

    //running internal chemspace
    mmpRunChemSpace(table, molecules, mmpRules, sp);

    return new MmpAnalysis(table, molecules, mmpRules, activityPairs, allPairsGrid, casesGrid, tp, sp, sliderInput);
  }

  async refreshPair() {
    //TODO: pin columns
    const idxParent = this.parentTable.currentRowIdx;
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
    const substructures = await getInverseSubstructures(pairsFrom, pairsTo);
    for (let i = 0; i < cases.length; i++) {
      diffFromSubstrCol.set(cases[i], substructures[0][i]);
      diffToSubstrCol.set(cases[i], substructures[1][i]);
    }

    this.casesGrid.dataFrame.filter.copyFrom(this.pairsMask);
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
      if (cutoff == 0) {
        this.cutoffMask.setAll(true);
        this.parentTable.filter.copyFrom(this.cutoffMask);
      } else {
        this.cutoffMask.setAll(false);
        for (let i = 0; i < this.activityPairs.length; i++) {
          if (this.activityPairs[i].diff >= cutoff) {
            this.cutoffMask.set(this.activityPairs[i].first, true, false);
            this.cutoffMask.set(this.activityPairs[i].second, true, false);
          }
        }
        this.parentTable.filter.copyFrom(this.cutoffMask);
      }
    } else
      this.parentTable.filter.copyFrom(this.cutoffMask);
  }
}
