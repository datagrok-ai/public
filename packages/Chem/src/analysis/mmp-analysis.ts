import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {drawMoleculeToCanvas, getRdKitModule, getRdKitService, getUncommonAtomsAndBonds} from '../utils/chem-common-rdkit';
import {getMCS} from '../utils/most-common-subs';

import {DimReductionMethods} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {BitArrayMetrics, BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {chemSpace} from './chem-space';
import $ from 'cash-dom';
import {ISubstruct} from '../rendering/rdkit-cell-renderer';
import {HIGHLIGHT_BY_SCAFFOLD_COL, SCAFFOLD_COL} from '../constants';

const MMP_VIEW_NAME = 'MMP Analysis';
const MMP_TAB_TRANSFORMATIONS = 'Transformations';
const MMP_TAB_FRAGMENTS = 'Fragments';
const MMP_TAB_CLIFFS = 'Cliffs';
const MMP_STRUCT_DIFF_NAME = '~structDiff';

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
  const res = await service.getFragments(molecules.toList());
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

function getCommonSubstructures(from: string[], to: string[]): string[] {
  //TODO: the length is known
  // const res1: (ISubstruct | null)[] = [];
  // const res2: (ISubstruct | null)[] = [];
  const res: string[] = [];
  for (let i = 0; i < from.length; i++) {
    const mcsCol = DG.Column.fromStrings('mcs', [from[i], to[i]]);
    const module = getRdKitModule();
    const mcs = getMCS(mcsCol, true, true);
    const mcsMol = module.get_qmol(mcs!);
    const substruct1 = getUncommonAtomsAndBonds(from[i], mcsMol, module);
    const substruct2 = getUncommonAtomsAndBonds(to[i], mcsMol, module);

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
    // res1.push(substruct1);
    // res2.push(substruct2);
    res.push(mcs === null ? '' : mcs);
    mcsMol.delete;
    // mol1.delete;
    // mol2.delete;
  }
  return res;
}

export class MmpAnalysis {
  parentTable: DG.DataFrame = {} as DG.DataFrame;
  parentCol: DG.Column = {} as DG.Column;
  mmpRules: MmpRules = {rules: [], smilesFrags: []};
  mmpView: DG.View = {} as DG.View;
  enableFilters: boolean = true;
  //transformations tab objects
  allPairsGrid: DG.Grid = {} as DG.Grid;
  casesGrid: DG.Grid = {} as DG.Grid;
  canvasMol1: HTMLCanvasElement = {} as HTMLCanvasElement;
  canvasMol2: HTMLCanvasElement = {} as HTMLCanvasElement;
  transformationsMask: DG.BitSet = {} as DG.BitSet;
  pairsMask: DG.BitSet = {} as DG.BitSet;
  pairedTransformations: DG.DataFrame = {} as DG.DataFrame;
  //cliffs tab objects
  cliffsCutoff = 0;
  activityPairs: {first: number, second: number, diff: number}[] = [];
  cutoffMask: DG.BitSet = {} as DG.BitSet;

  async init(table: DG.DataFrame, molecules: DG.Column, activities:DG.Column) {
    //initial calculations
    this.parentTable = table;
    this.parentCol = molecules;
    const frags = await getMmpFrags(molecules);
    this.mmpRules = getMmpRules(frags);
    this.canvasMol1 = ui.canvas();
    this.canvasMol2 = ui.canvas();

    //Transformations tab - calculation of transformations
    const from = new Array<string>(this.mmpRules.rules.length);
    const to = new Array<string>(this.mmpRules.rules.length);
    const occasions = new Array<number>(this.mmpRules.rules.length);
    const meanDiff = new Array<number>(this.mmpRules.rules.length);

    let maxAct = 0;

    //set activity differences
    for (let i = 0; i < this.mmpRules.rules.length; i++) {
      from[i] = this.mmpRules.smilesFrags[this.mmpRules.rules[i].smilesRule1];
      to[i] = this.mmpRules.smilesFrags[this.mmpRules.rules[i].smilesRule2];
      occasions[i] = this.mmpRules.rules[i].pairs.length;

      let mean: number = 0;
      for (let j = 0; j < occasions[i]; j++) {
        const idx1 = this.mmpRules.rules[i].pairs[j].firstStructure;
        const idx2 = this.mmpRules.rules[i].pairs[j].secondStructure;
        const val1 = activities.get(idx1);
        const val2 = activities.get(idx2);
        const diff = val2 - val1;
        if (diff > 0) {
          if (diff > maxAct)
            maxAct = diff;

          this.activityPairs.push({first: idx1, second: idx2, diff: diff});
        }

        mean += diff;
      }
      mean/= occasions[i];
      meanDiff[i] = mean;
    }

    const fromCol = DG.Column.fromList('string', 'From', from);
    const toCol = DG.Column.fromList('string', 'To', to);
    const occasionsCol = DG.Column.fromList('double', 'Pairs', occasions);
    const meanDiffCol = DG.Column.fromList('double', 'MeanDiff', meanDiff);

    fromCol.semType = DG.SEMTYPE.MOLECULE;
    toCol.semType = DG.SEMTYPE.MOLECULE;
    occasionsCol.semType = DG.TYPE.INT;

    const dfAllPairs = DG.DataFrame.fromColumns([fromCol, toCol, occasionsCol, meanDiffCol]);
    this.allPairsGrid = dfAllPairs.plot.grid();
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

    //Transformation tab
    const pairsFrom: string[] = [];
    const pairsTo: string[] = [];
    const singeDiff: number[] = [];

    for (let i = 0; i < this.mmpRules.rules.length; i++) {
      for (let j = 0; j < this.mmpRules.rules[i].pairs.length; j++) {
        const idx1 = this.mmpRules.rules[i].pairs[j].firstStructure;
        const idx2 = this.mmpRules.rules[i].pairs[j].secondStructure;
        pairsFrom.push(this.parentCol.get(idx1));
        pairsTo.push(this.parentCol.get(idx2));
        singeDiff.push(activities.get(idx2) - activities.get(idx1));
      }
    }

    //const structDiffs = getCommonSubstructures(pairsFrom, pairsTo);

    const pairsFromCol = DG.Column.fromStrings('From', pairsFrom);
    const pairsToCol = DG.Column.fromStrings('To', pairsTo);
    const singleDiffCol = DG.Column.fromList('double', 'Difference', singeDiff);
    const structureDiffCol = DG.Column.string(MMP_STRUCT_DIFF_NAME, pairsFrom.length);
    //DG.Column.fromStrings('~structDiff', structDiffs);
    pairsFromCol.semType = DG.SEMTYPE.MOLECULE;
    pairsToCol.semType = DG.SEMTYPE.MOLECULE;
    pairsFromCol.temp[SCAFFOLD_COL] = MMP_STRUCT_DIFF_NAME;
    pairsToCol.temp[SCAFFOLD_COL] = MMP_STRUCT_DIFF_NAME;
    pairsFromCol.temp[HIGHLIGHT_BY_SCAFFOLD_COL] = 'true';
    pairsToCol.temp[HIGHLIGHT_BY_SCAFFOLD_COL] = 'true';
    this.pairedTransformations = DG.DataFrame.fromColumns([pairsFromCol, pairsToCol, singleDiffCol, structureDiffCol]);
    this.casesGrid = this.pairedTransformations.plot.grid();
    this.pairsMask = DG.BitSet.create(this.pairedTransformations.rowCount);
    this.pairsMask.setAll(false);
    this.refreshPair();

    //Fragments tab
    const tp = DG.Viewer.fromType(DG.VIEWER.TRELLIS_PLOT, this.allPairsGrid.table, {
      xColumnNames: [this.allPairsGrid.table.columns.byIndex(0).name],
      yColumnNames: [this.allPairsGrid.table.columns.byIndex(1).name],
    });

    //Cliffs tab
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

    const chemSpaceParams = {
      seqCol: molecules,
      methodName: DimReductionMethods.UMAP,
      similarityMetric: BitArrayMetricsNames.Tanimoto as BitArrayMetrics,
      embedAxesNames: ['~X', '~Y'],
      options: {},
    };

    const progressBarSpace = DG.TaskBarProgressIndicator.create(`Running Chemical space...`);
    chemSpace(chemSpaceParams).then((res) => {
      const embeddings = res.coordinates;
      for (const col of embeddings)
        table.columns.replace(col.name, col);

      const lines: {id: number, mols: number[], a: number[], b: number[]}[] = [];
      for (let i = 0; i < this.mmpRules.rules.length; i++) {
        for (let j = 0; j < this.mmpRules.rules[i].pairs.length; j++) {
          const fs = this.mmpRules.rules[i].pairs[j].firstStructure;
          const ss = this.mmpRules.rules[i].pairs[j].secondStructure;
          lines.push({id: i, mols: [fs, ss], a: [], b: []});
        }
      }

      const linesRes = {lines: lines};

      sp!.onEvent('d4-before-draw-scene')
        .subscribe((_: any) => {
          renderLines(sp as DG.ScatterPlotViewer, '~X', '~Y', linesRes);
        });

      progressBarSpace.close();
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

  refreshPair() {
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
    const diffCol = this.pairedTransformations.getCol(MMP_STRUCT_DIFF_NAME);
    const diffFrom = this.pairedTransformations.getCol('From');
    const diffTo = this.pairedTransformations.getCol('To');

    //search specific rule
    for (let i = 0; i < this.mmpRules.rules.length; i++) {
      const first = this.mmpRules.rules[i].smilesRule1; //.pairs[j].firstStructure;
      const second = this.mmpRules.rules[i].smilesRule2; //.pairs[j].secondStructure;
      for (let j = 0; j < this.mmpRules.rules[i].pairs.length; j++) {
        if (ruleSmiNum1 == first && ruleSmiNum2 == second) {
          this.pairsMask.set(counter, true, false);
          if (diffCol.get(counter) === '')
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
    const substructures = getCommonSubstructures(pairsFrom, pairsTo);
    for (let i = 0; i < cases.length; i++)
      diffCol.set(cases[i], substructures[i]);

    this.pairedTransformations.filter.copyFrom(this.pairsMask);
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
