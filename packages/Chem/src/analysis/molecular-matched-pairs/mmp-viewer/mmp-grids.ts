import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {FRAGMENTS_GRID_HEADER_TOOLTIPS, MMP_NAMES, PAIRS_GRID_HEADER_TOOLTIPS,
  columnsDescriptions} from './mmp-constants';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {MMPA} from '../mmp-analysis/mmpa';
import {getRdKitModule} from '../../../utils/chem-common-rdkit';
import {createColWithDescription} from './mmp-generations';
import {MmpInput} from './mmp-viewer';
import {Subject, Subscription} from 'rxjs';
import {resizeGridColsSize} from '../../../utils/ui-utils';
import {MMPSubstructProvider} from './mmp-substruct-provider';
import {addSubstructProvider} from '@datagrok-libraries/chem-meta/src/types';
import {fillPairInfo} from './mmp-cliffs';

export class MmpPairedGrids {
  parentTable: DG.DataFrame;
  molColName: string;
  mmpa: MMPA;
  //fp for fragment pairs
  fpGrid: DG.Grid; //transformations tab specific to molecules and fragments tab the whole
  fpGridMessage: HTMLElement = ui.div('', 'chem-mmp-grid-info-message');
  fpCatsFrom: string [];
  fpCatsTo: string [];
  fpMaskByMolecule: DG.BitSet; //mask for a single molecule
  //mmp for matched molecula  pairs
  mmpGridTrans: DG.Grid;
  mmpGridTransMessage: HTMLElement = ui.div('', 'chem-mmp-grid-info-message');
  mmpMaskTrans: DG.BitSet;
  mmpMaskTransSelection: DG.BitSet;
  parentTableFragSelection: DG.BitSet;
  parentTableTransSelection: DG.BitSet;

  mmpGridFrag: DG.Grid;
  mmpGridFragMessage: HTMLElement = ui.div('', 'chem-mmp-grid-info-message');
  mmpMaskFrag: DG.BitSet;

  pairsGridCliffsTab: DG.Grid;
  pairsGridCliffsTabMessage: HTMLElement = ui.div('', 'chem-mmp-grid-info-message');
  pairsMaskCliffsTab: DG.BitSet;

  enableFilters: boolean = true;
  rdkit: RDModule;

  showErrorEvent: Subject<boolean> = new Subject();
  showEmptyPairsWarningEvent: Subject<boolean> = new Subject();
  filters: DG.FilterGroup | null = null;
  lastPairIdx: number | null = null;
  lastFragmentIdx: number | null = null;
  currentTab = MMP_NAMES.TAB_TRANSFORMATIONS;
  fragsShowAllMode = true;

  constructor(subs: Subscription[], mmpInput: MmpInput, mmpa: MMPA, activityMeanNames: string[] ) {
    this.rdkit = getRdKitModule();
    this.mmpa = mmpa;
    this.parentTable = mmpInput.table;
    this.molColName = mmpInput.molecules.name;
    this.fpGrid = getFragmetsPairsGrid(activityMeanNames, mmpa);
    this.fpCatsFrom = this.fpGrid.dataFrame.columns.byName(MMP_NAMES.FROM).categories;
    this.fpCatsTo = this.fpGrid.dataFrame.columns.byName(MMP_NAMES.TO).categories;

    this.mmpGridTrans = getMatchedPairsGrid(mmpa, this.rdkit);
    this.fpMaskByMolecule = DG.BitSet.create(this.fpGrid.dataFrame.rowCount);
    this.mmpMaskTrans = DG.BitSet.create(this.mmpGridTrans.dataFrame.rowCount);
    this.mmpMaskTransSelection = DG.BitSet.create(this.mmpGridTrans.dataFrame.rowCount);
    this.parentTableFragSelection = DG.BitSet.create(this.parentTable.rowCount);
    this.parentTableTransSelection = DG.BitSet.create(this.parentTable.rowCount);

    this.mmpGridFrag = this.mmpGridTrans.dataFrame.plot.grid();
    this.mmpMaskFrag = DG.BitSet.create(this.mmpGridFrag.dataFrame.rowCount);

    this.pairsGridCliffsTab = this.mmpGridTrans.dataFrame.plot.grid();
    this.pairsMaskCliffsTab = DG.BitSet.create(this.mmpGridFrag.dataFrame.rowCount).setAll(true);


    subs.push(DG.debounce(this.parentTable.onCurrentRowChanged, 100).subscribe(() => {
      if (this.parentTable!.currentRowIdx !== -1 && !this.fragsShowAllMode) {
        withTimedProgressIndicator('Refreshing pairs...', 200, async () => {
          this.refilterFragmentPairsByMolecule(true);
          await this.refreshMatchedPair();
        })
      }
    }));

    this.initialResizeGridColsSize(this.fpGrid);
    this.initialResizeGridColsSize(this.mmpGridTrans);
    this.initialResizeGridColsSize(this.mmpGridFrag);
    this.initialResizeGridColsSize(this.pairsGridCliffsTab);
    this.updateInfoMessage(this.fpGrid, this.fpGridMessage, 'fragment pair');
    this.updateInfoMessage(this.mmpGridTrans, this.mmpGridTransMessage, 'molecule pair');
    this.updateInfoMessage(this.mmpGridFrag, this.mmpGridFragMessage, 'molecule pair');
    this.updateInfoMessage(this.pairsGridCliffsTab, this.pairsGridCliffsTabMessage, 'molecule pair');
    this.fpGrid.sort([activityMeanNames[0]], [false]);
  }

  initialResizeGridColsSize(grid: DG.Grid) {
    const gridSub = grid.onAfterDrawContent.subscribe(() => {
      resizeGridColsSize(grid, [MMP_NAMES.FROM, MMP_NAMES.TO], 150, 70);
      gridSub.unsubscribe();
    });
  }

  disableFiltersGroup() {
    if (this.filters)
      for (let f of this.filters.filters)
        this.filters.setEnabled(f, false);
  }

  enableFiltersGroup() {
    if (this.filters)
      for (let f of this.filters.filters)
        this.filters.setEnabled(f, true);
  }

  updateInfoMessage(grid: DG.Grid, div: HTMLElement, compName: string) {
    DG.debounce(grid.dataFrame.onFilterChanged, 100).subscribe(() => {
      ui.empty(div);
      const num = grid.dataFrame.filter.trueCount;
      div.append(ui.divText(`Showing ${num} ${compName}${num === 1 ? '' : 's'} of ${grid.dataFrame.rowCount}`));
    });
    const gridSub = grid.onAfterDrawContent.subscribe(() => {
      resizeGridColsSize(grid, [MMP_NAMES.FROM, MMP_NAMES.TO], 150, 70);
      gridSub.unsubscribe();
    });
  }

  setupGrids(): void {
    this.fpMaskByMolecule!.setAll(true);

    this.fpGrid.table.onCurrentRowChanged.subscribe(() => {
      if (this.fpGrid!.table.currentRowIdx !== -1)
        withTimedProgressIndicator('Refreshing pairs...', 200, () => this.refreshMatchedPair());
    });

    this.fpGrid.table.onSelectionChanged.subscribe(() => {
      this.selectPairsWithSubstitutionInParentTable(false);
    });

    this.fpGrid.table.onFilterChanged.subscribe((e: any) => {
      //show error in case grid is empty
      this.showErrorEvent.next(this.fpGrid.dataFrame.filter.trueCount === 0);
    });

    //when swithcing from Fragments to Substitutions tab the following steps are triggered:
    //disable of filters from filters panel -> collaborative filtering -> onRowsFiltering (in case we are aleady on Substitution tab, we set the required mask) -> 
    // -> in onTabChanged we call requestFilter() for fpGrid.dataFrame to trigger filtration
    this.fpGrid.dataFrame.onRowsFiltering.subscribe((_) => {
      if (this.currentTab === MMP_NAMES.TAB_TRANSFORMATIONS)
        this.fpGrid!.dataFrame.filter.copyFrom(this.fpMaskByMolecule!);
    });

    this.mmpGridTrans.table.onCurrentRowChanged.subscribe(() => {
      if (this.mmpGridTrans.table.currentRowIdx !== -1) {
        this.pinMatchedPair(this.mmpGridTrans.table.currentRowIdx, this.mmpGridTrans);
        setTimeout(() => {
          grok.shell.windows.showContextPanel = true;
          grok.shell.o = fillPairInfo(this.mmpa!, this.mmpGridTrans.table.currentRowIdx,
            this.mmpGridTrans.table, this.parentTable!, this.rdkit, this.molColName);
        }, 500);
      }
    });

    this.mmpGridTrans.table.onSelectionChanged.subscribe(() => {
      this.selectPairsWithSubstitutionInParentTable(true);
    });

    this.mmpGridTrans.table.onFilterChanged.subscribe(() => {
      this.showEmptyPairsWarningEvent.next(this.mmpGridTrans.dataFrame.filter.trueCount === 0);
    });

    this.createCustomGridTooltips(this.fpGrid, FRAGMENTS_GRID_HEADER_TOOLTIPS, true);
    this.createCustomGridTooltips(this.mmpGridTrans, PAIRS_GRID_HEADER_TOOLTIPS);

    this.mmpMaskTrans.setAll(true);
  }

  createCustomGridTooltips(grid: DG.Grid, tooltips: {[key: string]: string}, mean?: boolean) {
    grid.onCellTooltip(function(cell: DG.GridCell, x: number, y: number) {
      if (cell.isColHeader && cell.tableColumn) {
        let tooltip = '';
        if (tooltips[cell.tableColumn.name])
          tooltip = tooltips[cell.tableColumn.name];
        else
          tooltip = cell.tableColumn.name.replace('\u0394', mean ? 'Mean difference in' : 'Difference in');
        ui.tooltip.show(ui.divText(tooltip), x, y);
        return true;
      } else
        return false;
    });
  };

  refreshMaskFragmentPairsFilter(): void {
    this.enableFilters = false;
    this.mmpMaskFrag.setAll(true);
    this.unPinMatchedPair();
    this.mmpGridFrag.dataFrame.filter.setAll(true);
  };

  /**
  * Gets visible all fragment pairs according to selected molecule in the table.
  * @param {boolean} rowChanged flag if row was changed
  */
  refilterFragmentPairsByMolecule(rowChanged: boolean) : void {
    let idxTrue = -1;
    if (rowChanged) {
      const idx = this.parentTable!.currentRowIdx;
      this.fpMaskByMolecule!.setAll(false);
      if (idx !== -1) {
        for (let i = 0; i < this.mmpa!.rules.rules.length; i++) {
          for (let j = 0; j < this.mmpa!.rules.rules[i].pairs.length; j++) {
            const fs = this.mmpa!.rules.rules[i].pairs[j].fs;
            if (idx == fs) {
              if (idxTrue == -1)
                idxTrue = i;
              this.fpMaskByMolecule!.set(i, true, false);
              break;
            }
          }
        }
      }
    }

    if (this.enableFilters) {
      this.fpGrid!.dataFrame.filter.copyFrom(this.fpMaskByMolecule!);
      this.fpGrid!.invalidate();
      if (rowChanged) {
        this.fpGrid!.table.currentRowIdx = idxTrue;
        this.lastFragmentIdx = idxTrue;
      }
    }
  }

  /**
  * Prepares all the entities to show for selected pair.
  */
  async refreshMatchedPair() : Promise<void> {
    this.mmpMaskTrans.setAll(false);

    const [idxPairs] = this.findSpecificRule(this.fragsShowAllMode);

    if (idxPairs >= 0) {
      if (this.lastPairIdx !== null)
        this.mmpGridTrans.dataFrame.col(MMP_NAMES.PAIR_SORT)?.set(this.lastPairIdx, false);
      this.mmpGridTrans.dataFrame.col(MMP_NAMES.PAIR_SORT)?.set(idxPairs, true);
      this.lastPairIdx = idxPairs;
      this.mmpGridTrans.sort([MMP_NAMES.PAIR_SORT], [false]);
    }

    this.mmpGridTrans.dataFrame.filter.copyFrom(this.mmpMaskTrans);
  }

  /**
  * pin the pair from "Pairs" dataframe to parent table
  * @param {number} idx index of pair
  */
  pinMatchedPair(idx: number, grid: DG.Grid): void {
    const columns = grid.dataFrame.columns;
    const idxFrom: number = columns!.byName(MMP_NAMES.PAIRNUM_FROM).get(idx);
    const idxToTo: number = columns!.byName(MMP_NAMES.PAIRNUM_TO).get(idx);
    const molFrom = this.parentTable!.columns.byName(this.mmpa.initData.molName).get(idxFrom);
    const molTo = this.parentTable!.columns.byName(this.mmpa.initData.molName).get(idxToTo);
    const gridt = grok.shell.tv.grid;
    gridt.setOptions({
      pinnedRowValues: [molFrom, molTo],
      pinnedRowColumnNames: [this.mmpa.initData.molName, this.mmpa.initData.molName],
    });
  }


  selectPairsWithSubstitutionInParentTable(pairsGrid: boolean) {
    if (!pairsGrid) {
      this.mmpMaskTransSelection.setAll(false);
      //path from selections in fpGrid to selected rows in mmpGridTrans
      const selectedRows = this.fpGrid.table.selection.getSelectedIndexes();
      this.findSpecificRule(this.fragsShowAllMode, selectedRows);
    }
    const mask = pairsGrid ? this.parentTableTransSelection : this.parentTableFragSelection;

    //reset previous selection from mask
    for (let i = -1; (i = mask.findNext(i, true)) != -1;)
      this.parentTable.selection.set(i, false);

    const columns = this.mmpGridTrans.dataFrame.columns;
    const gridt = grok.shell.tv.grid;
    const newSelection = DG.BitSet.create(this.parentTable.rowCount);
    const indexesPairs = pairsGrid ? this.mmpGridTrans.dataFrame.selection.getSelectedIndexes() :
      this.mmpMaskTransSelection.getSelectedIndexes();
    const indexesAllFrom = indexesPairs?.map((ip) => columns!.byName(MMP_NAMES.PAIRNUM_FROM).get(ip));
    const indexesAllTo = indexesPairs?.map((ip) => columns!.byName(MMP_NAMES.PAIRNUM_TO).get(ip));
    indexesAllFrom?.forEach((i) => newSelection.set(i, true));
    indexesAllTo?.forEach((i) => newSelection.set(i, true));
    mask.copyFrom(newSelection);
    gridt.dataFrame.selection.or(newSelection);
  }

  unPinMatchedPair(): void {
    const grid = grok.shell.tv?.grid;
    if (grid) {
      grid.setOptions({
        pinnedRowValues: [],
        pinnedRowColumnNames: [],
      });
    }
  }

  findSpecificRule(showAllFragsMode: boolean,
    selectionIdxs?: Int32Array): [idxPairs: number] {
    const idxParent = this.parentTable!.currentRowIdx;
    let idxPairs = -1;
    const idxs = selectionIdxs ?? [this.fpGrid!.table.currentRowIdx];
    const maskToModify = selectionIdxs ? this.mmpMaskTransSelection : this.mmpMaskTrans;
    if (idxParent !== -1 || showAllFragsMode) {
      for (const idx of idxs) {
        if (idx !== -1) {
          const ruleSmiNum1 = this.mmpa.rules.rules[idx].sr1;
          const ruleSmiNum2 = this.mmpa.rules.rules[idx].sr2;

          let counter = 0;

          for (let i = 0; i < this.mmpa.rules.rules.length; i++) {
            const first = this.mmpa.rules.rules[i].sr1;
            const second = this.mmpa.rules.rules[i].sr2;
            for (let j = 0; j < this.mmpa.rules.rules[i].pairs.length; j++) {
              if (ruleSmiNum1 == first && ruleSmiNum2 == second) {
                maskToModify.set(counter, true, false);
                if (!selectionIdxs) {
                  if (!showAllFragsMode && this.mmpa.rules.rules[i].pairs[j].fs == idxParent)
                    idxPairs = counter;
                }
              }
              counter++;
            }
          }
        }
      }
    }
    return [idxPairs];
  }

  refilterMatchedPairsByFragments(cats: [number, number]): void {
    this.mmpMaskFrag.setAll(false);

    let idx = -1;
    const colFrom = this.fpGrid.dataFrame.columns.byName(MMP_NAMES.FROM);
    const colTo = this.fpGrid.dataFrame.columns.byName(MMP_NAMES.TO);
    const fFrom = this.fpCatsFrom[cats[0]];
    const fTo = this.fpCatsTo[cats[1]];

    for (let i = 0; i < this.fpGrid.dataFrame.rowCount; i++) {
      if (colFrom.get(i) == fFrom && colTo.get(i) == fTo) {
        idx = i;
        break;
      }
    }

    if (idx !== -1) {
      const ruleSmiNum1 = this.mmpa.rules.rules[idx].sr1;
      const ruleSmiNum2 = this.mmpa.rules.rules[idx].sr2;

      let counter = 0;

      for (let i = 0; i < this.mmpa.rules.rules.length; i++) {
        const first = this.mmpa.rules.rules[i].sr1;
        const second = this.mmpa.rules.rules[i].sr2;
        for (let j = 0; j < this.mmpa.rules.rules[i].pairs.length; j++) {
          if (ruleSmiNum1 == first && ruleSmiNum2 == second)
            this.mmpMaskFrag.set(counter, true, false);
          counter++;
        }
      }
    }

    this.mmpGridFrag.dataFrame.filter.copyFrom(this.mmpMaskFrag);
  }
}

function getFragmetsPairsGrid(activityMeanNames: string[], mmpa: MMPA) : DG.Grid {
  const fromCol = createColWithDescription('string', MMP_NAMES.FROM, mmpa.rulesBased.fromFrag);
  const toCol = createColWithDescription('string', MMP_NAMES.TO, mmpa.rulesBased.toFrag);
  const occasionsCol = DG.Column.fromInt32Array(MMP_NAMES.PAIRS_COUNT, mmpa.rulesBased.occasions);
  occasionsCol.setTag('description', columnsDescriptions[MMP_NAMES.PAIRS_COUNT]);
  fromCol.semType = DG.SEMTYPE.MOLECULE;
  toCol.semType = DG.SEMTYPE.MOLECULE;
  occasionsCol.semType = DG.TYPE.INT;

  const fpCols = [fromCol, toCol];

  for (let i = 0; i < activityMeanNames.length; i++)
    fpCols.push(DG.Column.fromFloat32Array(activityMeanNames[i], mmpa.rulesBased.meanDiffs[i]));
  fpCols.push(occasionsCol);

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

function getMatchedPairsGrid(mmpa: MMPA, rdkit: RDModule) : DG.Grid {
  const pairsFromCol = createColWithDescription('string', MMP_NAMES.FROM, mmpa.allCasesBased.molFrom);
  const pairsToCol = createColWithDescription('string', MMP_NAMES.TO, mmpa.allCasesBased.molTo);
  const pairNumberCol = DG.Column.fromInt32Array(MMP_NAMES.PAIRNUM, mmpa.allCasesBased.pairNum);
  const pairNumberSortCol = DG.Column.bool(MMP_NAMES.PAIR_SORT, mmpa.allCasesBased.pairNum.length);
  const pairNumberFromCol = DG.Column.fromInt32Array(MMP_NAMES.PAIRNUM_FROM, mmpa.allCasesBased.molNumFrom);
  const pairNumberToCol = DG.Column.fromInt32Array(MMP_NAMES.PAIRNUM_TO, mmpa.allCasesBased.molNumTo);

  const pairsFromSmilesCol = DG.Column.fromStrings(MMP_NAMES.SMI1, mmpa.allCasesBased.pairsFromSmiles);
  const pairsToSmilesCol = DG.Column.fromStrings(MMP_NAMES.SMI2, mmpa.allCasesBased.pairsToSmiles);
  const ruleNumCol = DG.Column.fromInt32Array(MMP_NAMES.RULENUM, mmpa.allCasesBased.ruleNum);
  const coreNumCol = DG.Column.fromInt32Array(MMP_NAMES.CORE_NUM, mmpa.allCasesBased.coreNums);

  pairsFromSmilesCol.semType = DG.SEMTYPE.MOLECULE;
  pairsToSmilesCol.semType = DG.SEMTYPE.MOLECULE;
  pairsFromCol.semType = DG.SEMTYPE.MOLECULE;
  pairsToCol.semType = DG.SEMTYPE.MOLECULE;

  const allTransformationsCols = [pairsFromCol, pairsToCol, pairNumberCol, pairNumberFromCol, pairNumberToCol,
    pairsFromSmilesCol, pairsToSmilesCol, ruleNumCol, pairNumberSortCol, coreNumCol];

  for (let i = 0; i < mmpa.initData.activitiesCount; i++) {
    const name = MMP_NAMES.DIFF + ' ' + mmpa.initData.activitiesNames[i];
    allTransformationsCols.push(DG.Column.fromFloat32Array(name, mmpa.allCasesBased.diffs[i]));
  }

  const pairedTransformations = DG.DataFrame.fromColumns(allTransformationsCols);

  const substructProviderFrom = new MMPSubstructProvider(mmpa.frags.idToName,
    pairedTransformations, MMP_NAMES.FROM, rdkit, '#bc131f');
  addSubstructProvider(pairsFromCol.temp, substructProviderFrom);

  const substructProviderTo = new MMPSubstructProvider(mmpa.frags.idToName,
    pairedTransformations, MMP_NAMES.TO, rdkit, '#49bead');
  addSubstructProvider(pairsToCol.temp, substructProviderTo);

  return pairedTransformations.plot.grid();
}

/** Creates progress indicator only when inner function is taking more than n ms.
 */
async function withTimedProgressIndicator(pgName: string, time: number, func: () => Promise<void>) {
  let pg: DG.TaskBarProgressIndicator | null = null;
  const timer = setTimeout(() => {
    pg = DG.TaskBarProgressIndicator.create(pgName);
  })
  try {
    await func();
  } catch (e) {
    console.error(e);
  }
  clearTimeout(timer);
  if (pg)
    (pg as unknown as DG.TaskBarProgressIndicator).close();
}