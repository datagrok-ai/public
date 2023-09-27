import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {drawMoleculeToCanvas, getRdKitModule, getUncommonAtomsAndBonds} from '../utils/chem-common-rdkit';
import {getMCS} from '../utils/most-common-subs';

type MmpRules = {
  rules: {
    smilesRule1: number,
    smilesRule2: number,
    pairs: {firstStructure: number, secondStructure: number}[]
  } [],
  smilesFrags: string[]
};

function getMmpFrags(molecules: DG.Column):[string, string][][] {
  const module = getRdKitModule();
  const frags: [string, string][][] = new Array<[string, string][]>(molecules.length);
  for (let i = 0; i < molecules.length; i++) {
    let mol;
    try {
      mol = module.get_mol(molecules.get(i));
      if (mol) {
        const res = module.get_mmp(mol, 1, 1, 20);
        const fSplit = res.split(';');
        const ffSplit = fSplit[1].split(',');
        ffSplit.pop();
        frags[i] = new Array<[string, string]>(ffSplit.length);

        for (let j = 0; j < ffSplit.length; j++) {
          const fffSplit = ffSplit[j].split('.');
          const firstIsFirst = fffSplit[0].length >= fffSplit[1].length;
          frags[i][j] = [firstIsFirst ? fffSplit[0] : fffSplit[1], firstIsFirst ? fffSplit[1] : fffSplit[0]];
        }
      } else
        frags[i] = new Array<[string, string]>(0);
    } catch (e: any) {
      frags[i] = new Array<[string, string]>(0);
    } finally {
      mol?.delete();
    }
  }

  return frags;
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

export class MmpAnalysis {
  parentTable: DG.DataFrame;
  parentCol: DG.Column;
  mmpRules: MmpRules;
  allPairsGrid: DG.Grid;
  mmpView: DG.View;
  canvasMol1: HTMLCanvasElement;
  canvasMol2: HTMLCanvasElement;

  constructor(table: DG.DataFrame, molecules: DG.Column, activities:DG.Column) {
    this.parentTable = table;
    this.parentCol = molecules;
    const frags = getMmpFrags(molecules);
    this.mmpRules = getMmpRules(frags);
    this.canvasMol1 = ui.canvas();
    this.canvasMol2 = ui.canvas();

    const from = new Array<string>(this.mmpRules.rules.length);
    const to = new Array<string>(this.mmpRules.rules.length);
    const occasions = new Array<number>(this.mmpRules.rules.length);
    const meanDiff = new Array<number>(this.mmpRules.rules.length);

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
        mean += val2 - val1;
      }
      mean/= occasions[i];
      meanDiff[i] = mean;
    }

    const fromCol = DG.Column.fromList('string', 'From', from);
    const toCol = DG.Column.fromList('string', 'To', to);
    const occasionsCol = DG.Column.fromList('double', 'Occasions', occasions);
    const meanDiffCol = DG.Column.fromList('double', 'MeanDiff', meanDiff);

    fromCol.semType = DG.SEMTYPE.MOLECULE;
    toCol.semType = DG.SEMTYPE.MOLECULE;

    const dfAllPairs = DG.DataFrame.fromColumns([fromCol, toCol, occasionsCol, meanDiffCol]);
    this.allPairsGrid = dfAllPairs.plot.grid();

    this.parentTable.onCurrentRowChanged.subscribe(() => {
      this.refilterAllPairs();
      this.refreshPair();
    });

    this.refilterAllPairs();

    this.allPairsGrid.table.currentRowIdx = 0;
    this.allPairsGrid.table.onCurrentRowChanged.subscribe(() => {
      this.refreshPair();
    });

    this.refreshPair();

    this.mmpView = DG.View.create();
    this.mmpView.name = 'MMP Analysis';
    this.mmpView.box = true;
    const tp = DG.Viewer.fromType(DG.VIEWER.TRELLIS_PLOT, this.allPairsGrid.table, {
      xColumnNames: [this.allPairsGrid.table.columns.byIndex(0).name],
      yColumnNames: [this.allPairsGrid.table.columns.byIndex(1).name],
    });
    //const tp = DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, grid.table);
    const tabs = ui.tabControl(null, false);

    tabs.addPane('Transformations', () => {
      //this.refilterAllPairs();

      const molCells = ui.divH([this.canvasMol1, this.canvasMol2]);
      //box.style.maxHeight = '0px';

      return ui.splitV([molCells, this.allPairsGrid.root]);
    });
    tabs.addPane('Fragments', () => {
      this.refreshFilterAllPairs();
      return tp.root;
    });
    this.mmpView.append(tabs);
  }

  refreshPair() {
    const idxParent = this.parentTable.currentRowIdx;
    const idx = this.allPairsGrid.table.currentRowIdx;

    //TODO: add additional structure for algorithm simplification

    const ruleSmi1 = this.allPairsGrid.table.getCol('From').get(idx);
    const ruleSmi2 = this.allPairsGrid.table.getCol('To').get(idx);

    //search for numbers in smiles rules
    const ruleSmiNum1 = this.mmpRules.smilesFrags.indexOf(ruleSmi1);
    const ruleSmiNum2 = this.mmpRules.smilesFrags.indexOf(ruleSmi2);

    let idxSecond = -1;

    //search specific rule
    for (let i = 0; i < this.mmpRules.rules.length; i++) {
      const first = this.mmpRules.rules[i].smilesRule1; //.pairs[j].firstStructure;
      const second = this.mmpRules.rules[i].smilesRule2; //.pairs[j].secondStructure;
      if (ruleSmiNum1 == first && ruleSmiNum2 == second) {
        for (let j = 0; j < this.mmpRules.rules[i].pairs.length; j++) {
          const idxFirst = this.mmpRules.rules[i].pairs[j].firstStructure;

          if (idxFirst == idxParent) {
            idxSecond = this.mmpRules.rules[i].pairs[j].secondStructure;
            break;
          }
        }
      }
      if (idxSecond >= 0)
        break;
    }

    const molecule1 = this.parentCol.get(idxParent);
    const molecule2 = this.parentCol.get(idxSecond);
    const mcsCol = DG.Column.fromStrings('mcs', [molecule1, molecule2]);
    const module = getRdKitModule();
    const mcs = getMCS(mcsCol, true, true);
    const mcsMol = module.get_qmol(mcs!);
    const substruct1 = getUncommonAtomsAndBonds(molecule1, mcsMol, module);
    const substruct2 = getUncommonAtomsAndBonds(molecule2, mcsMol, module);
    // getUncommonAtomsAndBonds();


    drawMoleculeToCanvas(0, 0, 200, 100, this.canvasMol1, molecule1, '',
      {normalizeDepiction: true, straightenDepiction: true}, substruct1);
    drawMoleculeToCanvas(0, 0, 200, 100, this.canvasMol2, molecule2, '',
      {normalizeDepiction: true, straightenDepiction: true}, substruct2);
  }

  refreshFilterAllPairs() {
    const consistsBitSet: DG.BitSet = DG.BitSet.create(this.allPairsGrid.dataFrame.rowCount);
    consistsBitSet.setAll(true);
    this.allPairsGrid.dataFrame.filter.copyFrom(consistsBitSet);
  }

  refilterAllPairs() {
    const idx = this.parentTable.currentRowIdx;
    const consistsBitSet: DG.BitSet = DG.BitSet.create(this.allPairsGrid.dataFrame.rowCount);

    for (let i = 0; i < this.mmpRules.rules.length; i++) {
      for (let j = 0; j < this.mmpRules.rules[i].pairs.length; j++) {
        const fs = this.mmpRules.rules[i].pairs[j].firstStructure;
        if (idx == fs) {
          consistsBitSet.set(i, true, false);
          break;
        }
      }
    }

    this.allPairsGrid.dataFrame.filter.copyFrom(consistsBitSet);
    this.allPairsGrid.invalidate();
  }
}
