import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ChemTemps} from '@datagrok-libraries/chem-meta/src/consts';
import {MMP_NAMES, columnsDescriptions} from './mmp-constants';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {MMPA} from '../mmp-analysis/mmpa';
import {getRdKitModule} from '../../../utils/chem-common-rdkit';
import {createColWithDescription} from './mmp-generations';
import {debounceTime} from 'rxjs/operators';
import {getInverseSubstructuresAndAlign} from './mmp-mol-rendering';
import {MmpInput} from './mmp-viewer';

export class MmpPairedGrids {
  parentTable: DG.DataFrame;
  mmpa: MMPA;
  //fp for fragment pairs
  fpGrid: DG.Grid; //transformations tab specific to molecules and fragments tab the whole
  fpCatsFrom: string [];
  fpCatsTo: string [];
  fpMaskByMolecule: DG.BitSet; //mask for a single molecule
  fpMaskFilter: DG.BitSet; //mask for a number of cases filter
  //mmp for matched molecula  pairs
  mmpGrid: DG.Grid;
  mmpMaskByFragment: DG.BitSet;

  enableFilters: boolean = true;
  rdkit: RDModule;

  constructor(mmpInput: MmpInput, mmpa: MMPA, activityMeanNames: string[] ) {
    this.mmpa = mmpa;
    this.parentTable = mmpInput.table;
    this.fpGrid = getFragmetsPairsGrid(activityMeanNames, mmpa);
    this.fpCatsFrom = this.fpGrid.dataFrame.columns.byName(MMP_NAMES.FROM).categories;
    this.fpCatsTo = this.fpGrid.dataFrame.columns.byName(MMP_NAMES.TO).categories;

    this.mmpGrid = getMatchedPairsGrid(mmpa);
    this.fpMaskByMolecule = DG.BitSet.create(this.fpGrid.dataFrame.rowCount);
    this.fpMaskFilter = DG.BitSet.create(this.fpGrid.dataFrame.rowCount);
    this.mmpMaskByFragment = DG.BitSet.create(this.mmpGrid.dataFrame.rowCount);

    this.rdkit = getRdKitModule();
  }

  setupGrids(): void {
    this.fpMaskFilter!.setAll(true);
    this.fpMaskByMolecule!.setAll(true);
    this.parentTable.onCurrentRowChanged.pipe(debounceTime(1000)).subscribe(() => {
      if (this.parentTable!.currentRowIdx !== -1) {
        this.refilterFragmentPairsByMolecule(true);
        this.refreshMatchedPair(this.rdkit);
      }
    });

    this.refilterFragmentPairsByMolecule(true);

    this.fpGrid.table.onCurrentRowChanged.subscribe(() => {
      if (this.fpGrid!.table.currentRowIdx !== -1)
        this.refreshMatchedPair(this.rdkit);
    });

    this.mmpGrid.table.onCurrentRowChanged.subscribe(() => {
      if (this.mmpGrid!.table.currentRowIdx !== -1)
        this.pinMatchedPair(this.mmpGrid!.table.currentRowIdx);
    });

    this.mmpMaskByFragment!.setAll(false);
    this.refreshMatchedPair(this.rdkit);
  }

  refreshMaskFragmentPairsFilter(): void {
    this.fpGrid.dataFrame.filter.copyFrom(this.fpMaskFilter!);
  };

  refilterFragmentPairsByCases(value: number) : void {
    this.fpMaskFilter.setAll(false);

    const casesCol = this.fpGrid.dataFrame.columns.byName(MMP_NAMES.PAIRS);

    for (let i = 0; i < casesCol.length; i++) {
      if (casesCol.get(i) >= value)
        this.fpMaskFilter.set(i, true, false);
    }

    this.fpGrid.dataFrame.filter.copyFrom(this.fpMaskFilter!);
  }

  /**
  * Gets visible all fragment pairs according to selected molecule in the table.
  * @param {boolean} rowChanged flag if row was changed
  */
  refilterFragmentPairsByMolecule(rowChanged: boolean) : void {
    let idxTrue = -1;
    if (rowChanged) {
      const idx = this.parentTable!.currentRowIdx;
      this.fpMaskByMolecule!.setAll(false);

      for (let i = 0; i < this.mmpa!.rules.rules.length; i++) {
        for (let j = 0; j < this.mmpa!.rules.rules[i].pairs.length; j++) {
          const fs = this.mmpa!.rules.rules[i].pairs[j].firstStructure;
          if (idx == fs) {
            if (idxTrue == -1)
              idxTrue = i;
            this.fpMaskByMolecule!.set(i, true, false);
            break;
          }
        }
      }
    }

    if (this.enableFilters) {
      this.fpGrid!.dataFrame.filter.copyFrom(this.fpMaskByMolecule!);
      this.fpGrid!.invalidate();
      if (rowChanged)
        this.fpGrid!.table.currentRowIdx = idxTrue;
    }
  }

  /**
  * Prepares all the entities to show for selected pair.
  * @param {RDModule} rdkitModule RDkit module instance
  */
  async refreshMatchedPair(rdkitModule: RDModule) : Promise<void> {
    const progressBarPairs = DG.TaskBarProgressIndicator.create(`Refreshing pairs...`);

    this.mmpMaskByFragment.setAll(false);
    const diffFromSubstrCol = this.mmpGrid.dataFrame.getCol(MMP_NAMES.STRUCT_DIFF_FROM_NAME);
    const diffToSubstrCol = this.mmpGrid.dataFrame.getCol(MMP_NAMES.STRUCT_DIFF_TO_NAME);
    const diffFrom = this.mmpGrid.dataFrame.getCol(MMP_NAMES.FROM);
    const diffTo = this.mmpGrid.dataFrame.getCol(MMP_NAMES.TO);

    const [idxPairs, cases] = this.findSpecificRule(diffFromSubstrCol);
    this.recoverHighlights(cases, diffFrom, diffTo, diffFromSubstrCol, diffToSubstrCol, rdkitModule);

    this.mmpGrid.dataFrame.filter.copyFrom(this.mmpMaskByFragment);

    this.unPinMatchedPair();

    if (idxPairs >= 0) {
      this.mmpGrid.setOptions({
        pinnedRowValues: [idxPairs.toString()],
        pinnedRowColumnNames: [MMP_NAMES.PAIRNUM],
      });
    } else {
      this.mmpGrid.setOptions({
        pinnedRowValues: [],
        pinnedRowColumnNames: [],
      });
    }

    progressBarPairs.close();
  }

  /**
  * pin the pair from "Pairs" dataframe to parent table
  * @param {number} idx index of pair
  */
  pinMatchedPair(idx: number): void {
    const columns = this.mmpGrid?.dataFrame.columns;
    const idxFrom: number = columns!.byName(MMP_NAMES.PAIRNUM_FROM).get(idx);
    const idxToTo: number = columns!.byName(MMP_NAMES.PAIRNUM_TO).get(idx);
    const molFrom = this.parentTable!.columns.byName(this.mmpa.initData.molName).get(idxFrom);
    const molTo = this.parentTable!.columns.byName(this.mmpa.initData.molName).get(idxToTo);
    const grid = grok.shell.tv.grid ?? ((grok.shell.view('Browse')! as DG.BrowseView)!.preview! as DG.TableView).grid;
    const indexesPairs = this.mmpMaskByFragment?.getSelectedIndexes();
    const indexesAllFrom = indexesPairs?.map((ip) => columns!.byName(MMP_NAMES.PAIRNUM_FROM).get(ip));
    const indexesAllTo = indexesPairs?.map((ip) => columns!.byName(MMP_NAMES.PAIRNUM_TO).get(ip));
    indexesAllFrom?.forEach((i) => grid.dataFrame.selection.set(i, true));
    indexesAllTo?.forEach((i) => grid.dataFrame.selection.set(i, true));

    grid.setOptions({
      pinnedRowValues: [molFrom, molTo],
      pinnedRowColumnNames: [this.mmpa.initData.molName, this.mmpa.initData.molName],
    });
  }

  unPinMatchedPair(): void {
    const grid = grok.shell.tv.grid ?? ((grok.shell.view('Browse')! as DG.BrowseView)!.preview! as DG.TableView).grid;
    if (grid) {
      grid.setOptions({
        pinnedRowValues: [],
        pinnedRowColumnNames: [],
      });
    }
  }

  findSpecificRule(diffFromSubstrCol: DG.Column): [idxPairs: number, cases: number[]] {
    const idxParent = this.parentTable!.currentRowIdx;
    let idxPairs = -1;
    const cases: number[] = [];
    const idx = this.fpGrid!.table.currentRowIdx;
    if (idx !== -1) {
      const ruleSmi1 = this.fpGrid.table.getCol(MMP_NAMES.FROM).get(idx);
      const ruleSmi2 = this.fpGrid.table.getCol(MMP_NAMES.TO).get(idx);
      const ruleSmiNum1 = this.mmpa.rules.smilesFrags.indexOf(ruleSmi1);
      const ruleSmiNum2 = this.mmpa.rules.smilesFrags.indexOf(ruleSmi2);

      let counter = 0;

      for (let i = 0; i < this.mmpa.rules.rules.length; i++) {
        const first = this.mmpa.rules.rules[i].smilesRule1;
        const second = this.mmpa.rules.rules[i].smilesRule2;
        for (let j = 0; j < this.mmpa.rules.rules[i].pairs.length; j++) {
          if (ruleSmiNum1 == first && ruleSmiNum2 == second) {
            this.mmpMaskByFragment.set(counter, true, false);
            if (diffFromSubstrCol.get(counter) === null)
              cases.push(counter);
            if (this.mmpa.rules.rules[i].pairs[j].firstStructure == idxParent)
              idxPairs = counter;
          }
          counter++;
        }
      }
    }
    return [idxPairs, cases];
  }

  refilterMatchedPairsByFragments(cats: [string, string]): void {
    this.mmpGrid.setOptions({
      pinnedRowValues: [],
      pinnedRowColumnNames: [],
    });
    this.mmpMaskByFragment.setAll(true); //change to false

    const idx = -1;
    for (let i = 0; i < this.fpCatsFrom.length; i++) {
      if (cats[0] == this.fpCatsFrom[i] && cats[1] == this.fpCatsTo[i]) {
        idx == i;
        break;
      }
    }


    if (idx !== -1) {
      const ruleSmi1 = this.fpGrid!.table.getCol(MMP_NAMES.FROM).get(idx);
      const ruleSmi2 = this.fpGrid!.table.getCol(MMP_NAMES.TO).get(idx);
      const ruleSmiNum1 = this.mmpa!.rules.smilesFrags.indexOf(ruleSmi1);
      const ruleSmiNum2 = this.mmpa!.rules.smilesFrags.indexOf(ruleSmi2);

      let counter = 0;

      for (let i = 0; i < this.mmpa!.rules.rules.length; i++) {
        const first = this.mmpa!.rules.rules[i].smilesRule1;
        const second = this.mmpa!.rules.rules[i].smilesRule2;
        for (let j = 0; j < this.mmpa!.rules.rules[i].pairs.length; j++) {
          if (ruleSmiNum1 == first && ruleSmiNum2 == second)
            this.mmpMaskByFragment.set(counter, true, false);
          counter++;
        }
      }
    }

    this.mmpGrid.dataFrame.filter.copyFrom(this.mmpMaskByFragment);
  }

  async recoverHighlights(
    cases: number[], diffFrom : DG.Column, diffTo: DG.Column,
    diffFromSubstrCol: DG.Column, diffToSubstrCol: DG.Column, rdkitModule: RDModule) : Promise<void> {
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
  }
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
  pairsFromCol.temp[ChemTemps.SUBSTRUCT_COL] = MMP_NAMES.STRUCT_DIFF_FROM_NAME;
  pairsToCol.temp[ChemTemps.SUBSTRUCT_COL] = MMP_NAMES.STRUCT_DIFF_TO_NAME;

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
