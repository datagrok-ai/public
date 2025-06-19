import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {FRAGMENTS_GRID_HEADER_DESCRIPTIONS, MMP_NAMES, PAIRS_GRID_HEADER_DESCRIPTIONS} from './mmp-constants';
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
import {MmmpaActivityRenderer} from '../mmpa-activity-cell-renderer';
import {PaletteCodes} from './palette';

const MEAN_ACTIVITY_COL_WIDTH = 100;

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
  parentFragmentsFilter: DG.BitSet;
  parentPairsFilter: DG.BitSet;

  mmpGridFrag: DG.Grid;
  mmpGridFragMessage: HTMLElement = ui.div('', 'chem-mmp-grid-info-message');
  mmpMaskFrag: DG.BitSet;

  pairsGridCliffsTab: DG.Grid;
  pairsGridCliffsTabMessage: HTMLElement = ui.div('', 'chem-mmp-grid-info-message');
  pairsMaskCliffsTab: DG.BitSet;

  rdkit: RDModule;

  showErrorEvent: Subject<boolean> = new Subject();
  showEmptyPairsWarningEvent: Subject<boolean> = new Subject();
  filters: DG.FilterGroup | null = null;
  lastPairIdx: number | null = null;
  lastFragmentIdx: number | null = null;
  currentTab = MMP_NAMES.TAB_TRANSFORMATIONS;
  fragsShowAllMode = true;
  currentRowChanged = false;
  followCurrentRowInFragmentsGrid = false;
  activityColors: PaletteCodes | null = null;
  contextPanelDiv = ui.divV([], 'mmpa-substitutions-context-panel');
  moleculePairsContextPanelDiv = ui.div();
  pcPlotContextPanelDiv = ui.div('', {style: {height: '100%'}});

  constructor(subs: Subscription[], mmpInput: MmpInput, mmpa: MMPA, activityMeanNames: string[],
    activityColors: PaletteCodes) {
    this.rdkit = getRdKitModule();
    this.mmpa = mmpa;
    this.activityColors = activityColors;
    this.parentTable = mmpInput.table;
    this.molColName = mmpInput.molecules.name;
    this.fpGrid = getFragmetsPairsGrid(activityMeanNames, mmpa, this.activityColors, subs);
    this.mmpGridTrans = getMatchedPairsGrid(mmpa, this.rdkit);
    this.fpCatsFrom = this.fpGrid.dataFrame.columns.byName(MMP_NAMES.FROM).categories;
    this.fpCatsTo = this.fpGrid.dataFrame.columns.byName(MMP_NAMES.TO).categories;
    this.mmpMaskTransSelection = DG.BitSet.create(this.mmpGridTrans.dataFrame.rowCount);
    this.parentTableFragSelection = DG.BitSet.create(this.parentTable.rowCount);
    this.parentTableTransSelection = DG.BitSet.create(this.parentTable.rowCount);

    this.mmpGridFrag = this.mmpGridTrans.dataFrame.plot.grid({allowColumnMenu: true});

    this.pairsGridCliffsTab = this.mmpGridTrans.dataFrame.plot.grid({allowColumnMenu: true});


    subs.push(this.fpGrid.onCellClick.subscribe((gc) => {
      const actIndex = activityMeanNames.indexOf(gc.cell?.column?.name);
      if (actIndex !== -1) {
        ui.empty(this.pcPlotContextPanelDiv);
        this.pcPlotContextPanelDiv.classList.remove('mmpa-hide-context-panel');
        const filter = `[${
          mmpa.rulesBased.moleculePairIndices[gc.cell?.rowIndex].join(',')
        }].contains(\${${MMP_NAMES.PAIRNUM}})`;
        const fromColName = `${mmpa.initData.activitiesNames[actIndex]}_from`;
        const toColName = `${mmpa.initData.activitiesNames[actIndex]}_to`;
        const title = `${activityMeanNames[actIndex]} for substitution #${gc.cell?.rowIndex}`;
        const pcPlot = DG.Viewer.fromType(DG.VIEWER.PC_PLOT, this.mmpGridTrans.dataFrame, {
          columnNames: [fromColName, toColName],
          filter: filter,
          showTitle: false,
        });
        pcPlot.root.classList.add('mmpa-pc-plot');
        pcPlot.root.prepend(ui.divText(title, 'mmpa-pc-plot-title'));
        this.pcPlotContextPanelDiv.append(pcPlot.root);
        //need to pass this true parameter since current row might still be -1
        //because it updates a bit later after onCellClick
        this.showContextPanel(true);
      }
    }));

    //Filter masks

    //fragment pairs for current molecule in parent df on substitutions tab
    this.fpMaskByMolecule = DG.BitSet.create(this.fpGrid.dataFrame.rowCount);

    //molecule pairs for current fragments pair on substitutions tab
    this.mmpMaskTrans = DG.BitSet.create(this.mmpGridTrans.dataFrame.rowCount);

    //selected cell in trellis plot on Fragments tab (molecule pairs df)
    this.mmpMaskFrag = DG.BitSet.create(this.mmpGridFrag.dataFrame.rowCount);

    //parent table filter (fragments df)
    this.parentFragmentsFilter = DG.BitSet.create(this.fpGrid.dataFrame.rowCount).setAll(true);

    //parent table filter (molecule pairs df)
    this.parentPairsFilter = DG.BitSet.create(this.mmpGridTrans.dataFrame.rowCount).setAll(true);

    //scatter plot activity filters (molecule pairs df) on Cliffs tab
    this.pairsMaskCliffsTab = DG.BitSet.create(this.mmpGridFrag.dataFrame.rowCount).setAll(true);

    this.contextPanelDiv.append(this.moleculePairsContextPanelDiv);
    this.contextPanelDiv.append(this.pcPlotContextPanelDiv);

    // Subscribe to parent table filter changes
    subs.push(DG.debounce(this.parentTable.onFilterChanged, 300).subscribe(() => {
      if (this.currentTab === MMP_NAMES.TAB_GENERATION)
        return;
      this.updateParentFilterMasks(this.currentTab === MMP_NAMES.TAB_CLIFFS);
      if (this.currentTab === MMP_NAMES.TAB_TRANSFORMATIONS) {
      /*   in case we are in current molecule mode and current molecule was
        filtered out - reset fragments and molecule pairs grids */
        if (!this.fragsShowAllMode && this.parentTable.currentRowIdx !== -1 &&
          !this.parentTable.filter.get(this.parentTable.currentRowIdx)) {
            this.fpGrid!.dataFrame.filter.setAll(false);
            this.mmpGridTrans.dataFrame.filter.setAll(false);
            return;
        }
      }
      this.fpGrid.dataFrame.rows.requestFilter();
      this.mmpGridTrans.dataFrame.rows.requestFilter();
    }));

    subs.push(DG.debounce(this.parentTable.onCurrentRowChanged, 100).subscribe(() => {
      if (this.currentTab === MMP_NAMES.TAB_TRANSFORMATIONS)
        this.refreshFragmentsAndPairs(true);
      else if (this.parentTable!.currentRowIdx !== -1 && !this.fragsShowAllMode)
        this.currentRowChanged = true;
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

  switchToSubstitutionsTab() {
    /* if current row was changed while we were on tab other than 'Substitution',
     re-calculate pairs (in case we are not in show all mode) */
    if (this.currentRowChanged) {
      this.currentRowChanged = false;
      this.refreshFragmentsAndPairs(true);
    }
    this.fpGrid.dataFrame.rows.requestFilter();
    this.mmpGridTrans.dataFrame.rows.requestFilter();
  }

  switchToFragmentsTab() {
    this.unPinMatchedPair();
    this.mmpGridTrans.dataFrame.rows.requestFilter();
  }

  async refreshFragmentsAndPairs(parentRowChanged?: boolean) {
    withTimedProgressIndicator('Refreshing pairs...', 200, async () => {
      if (parentRowChanged && this.parentTable!.currentRowIdx !== -1 && !this.fragsShowAllMode)
        this.refilterFragmentPairsByMolecule();
      if (this.followCurrentRowInFragmentsGrid)
        await this.refreshMatchedPair();
      else
        this.mmpMaskTrans.setAll(true);
      this.fpGrid.dataFrame.rows.requestFilter();
      this.mmpGridTrans.dataFrame.rows.requestFilter();
    });
  }

  updateFragmentsDfFilters() {
    if (this.currentTab === MMP_NAMES.TAB_TRANSFORMATIONS && this.fragsShowAllMode === false)
      this.fpGrid!.dataFrame.filter.and(this.fpMaskByMolecule!);
    this.fpGrid.dataFrame.filter.and(this.parentFragmentsFilter);
    this.fpGrid.invalidate();
  }

  updateMoleculePairsFilters(): void {
    const mask = this.currentTab === MMP_NAMES.TAB_TRANSFORMATIONS ? this.mmpMaskTrans :
      this.currentTab === MMP_NAMES.TAB_FRAGMENTS ? this.mmpMaskFrag : this.pairsMaskCliffsTab;
    this.mmpGridFrag.dataFrame.filter.and(mask);
    this.mmpGridFrag.dataFrame.filter.and(this.parentPairsFilter);
    this.currentTab === MMP_NAMES.TAB_TRANSFORMATIONS ? this.mmpGridTrans.invalidate() : this.mmpGridFrag.invalidate();
  };


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
      this.refreshFragmentsAndPairs(false);
    });

    this.fpGrid.table.onSelectionChanged.subscribe(() => {
      this.selectPairsWithSubstitutionInParentTable(false);
    });

    DG.debounce(this.fpGrid.table.onFilterChanged, 300).subscribe(() => {
      //show error in case grid is empty
      if (this.currentTab === MMP_NAMES.TAB_TRANSFORMATIONS)
        this.showErrorEvent.next(this.fpGrid.dataFrame.filter.trueCount === 0);
    });

    // when swithcing from Fragments to Substitutions tab the following steps are triggered:
    // disable of filters from filters panel -> collaborative filtering ->
    // -> onRowsFiltering (in case we are aleady on Substitution tab, we set the required mask)
    // in onTabChanged we call requestFilter() for fpGrid.dataFrame to trigger filtration
    this.fpGrid.dataFrame.onRowsFiltering.subscribe(() => {
      this.updateFragmentsDfFilters();
    });

    this.mmpGridTrans.dataFrame.onRowsFiltering.subscribe(() => {
      this.updateMoleculePairsFilters();
    });

    this.mmpGridTrans.table.onCurrentRowChanged.subscribe(() => {
      if (this.mmpGridTrans.table.currentRowIdx !== -1) {
        this.pinMatchedPair(this.mmpGridTrans.table.currentRowIdx, this.mmpGridTrans);
        setTimeout(() => {
          const pairInfo = fillPairInfo(this.mmpa!, this.mmpGridTrans.table.currentRowIdx,
            this.mmpGridTrans.table, this.parentTable!, this.rdkit, this.molColName);
          this.moleculePairsContextPanelDiv.classList.remove('mmpa-hide-context-panel');
          ui.empty(this.moleculePairsContextPanelDiv);
          this.moleculePairsContextPanelDiv.append(pairInfo);
          this.showContextPanel();
        }, 500);
      }
    });

    this.mmpGridTrans.table.onSelectionChanged.subscribe(() => {
      this.selectPairsWithSubstitutionInParentTable(true);
    });

    this.mmpGridTrans.table.onFilterChanged.subscribe(() => {
      this.showEmptyPairsWarningEvent.next(this.mmpGridTrans.dataFrame.filter.trueCount === 0);
    });

    this.createCustomGridTooltips(this.fpGrid, FRAGMENTS_GRID_HEADER_DESCRIPTIONS);
    this.createCustomGridTooltips(this.mmpGridTrans, PAIRS_GRID_HEADER_DESCRIPTIONS);

    this.mmpMaskTrans.setAll(true);
    this.mmpMaskFrag.setAll(true);
  }

  showContextPanel(show?: boolean) {
    const updateContextPanel = () => {
      if ((this.fpGrid.dataFrame.currentRowIdx === -1 && !show) || this.currentTab !== MMP_NAMES.TAB_TRANSFORMATIONS)
        this.pcPlotContextPanelDiv.classList.add('mmpa-hide-context-panel');
      if (this.mmpGridTrans.dataFrame.currentRowIdx === -1)
        this.moleculePairsContextPanelDiv.classList.add('mmpa-hide-context-panel');
      grok.shell.windows.showContextPanel = true;
      grok.shell.o = this.contextPanelDiv;
    };
    if (!grok.shell.windows.showContextPanel)
      updateContextPanel();
    else {
      try {
        if (!(grok.shell.o as HTMLElement).getElementsByClassName('mmpa-substitutions-context-panel').length)
          updateContextPanel();
      } catch (e: any) {
        updateContextPanel();
      }
    }
  }

  createCustomGridTooltips(grid: DG.Grid, tooltips: {[key: string]: string}) {
    grid.onCellTooltip(function(cell: DG.GridCell, x: number, y: number) {
      if (cell.isColHeader && cell.tableColumn && cell.tableColumn.semType === DG.SEMTYPE.MOLECULE) {
        const tooltip = tooltips[cell.tableColumn.name] ?? '';
        ui.tooltip.show(ui.divText(tooltip), x, y);
        return true;
      } else
        return false;
    });
  };

  refilterFragmentPairsByMolecule() : void {
    let idxTrue = -1;
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

  /**
  * Prepares all the entities to show for selected pair of fragments.
  */
  async refreshMatchedPair() : Promise<void> {
    this.mmpMaskTrans.setAll(false);

    //idxPairs  index of the pair with molecule from parent dataframe in it
    const [idxPairs] = this.findSpecificRule(this.fragsShowAllMode);

    if (idxPairs >= 0) {
      if (this.lastPairIdx !== null)
        this.mmpGridTrans.dataFrame.col(MMP_NAMES.PAIR_SORT)?.set(this.lastPairIdx, false);
      this.mmpGridTrans.dataFrame.col(MMP_NAMES.PAIR_SORT)?.set(idxPairs, true);
      this.lastPairIdx = idxPairs;
      this.mmpGridTrans.sort([MMP_NAMES.PAIR_SORT], [false]);
    }
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
  }

  updateParentFilterMasks(isCliffsTab?: boolean): void {
    this.updateFragmentsParentFilter();
    this.updatePairsParentFilter(isCliffsTab);
  }


  updateFragmentsParentFilter() {
    const parentFilter = this.parentTable.filter;

    this.parentFragmentsFilter.setAll(true);
    const fragmentParentIndices = this.mmpa.rulesBased.fragmentParentIndices;

    for (let i = 0; i < fragmentParentIndices.length; i++) {
      const parentIndices = fragmentParentIndices[i];
      let hasVisibleMolecule = false;
      for (const idx of parentIndices) {
        if (parentFilter.get(idx)) {
          hasVisibleMolecule = true;
          break;
        }
      }
      this.parentFragmentsFilter.set(i, hasVisibleMolecule);
    }
  }

  updatePairsParentFilter(cliffs?: boolean) {
    const parentFilter = this.parentTable.filter;

    this.parentPairsFilter.setAll(true);
    const [fromIndices, toIndices] = this.mmpa.rulesBased.moleculePairParentIndices;

    for (let i = 0; i < fromIndices.length; i++) {
      this.parentPairsFilter.set(i,
        cliffs ? parentFilter.get(fromIndices[i]) && parentFilter.get(toIndices[i]) :
          parentFilter.get(fromIndices[i]) || parentFilter.get(toIndices[i]));
    }
  }
}


function getFragmetsPairsGrid(activityMeanNames: string[], mmpa: MMPA, activityColors: PaletteCodes,
  subs: Subscription[]) : DG.Grid {
  const fromCol = createColWithDescription('string', MMP_NAMES.FROM,
    mmpa.rulesBased.fromFrag, FRAGMENTS_GRID_HEADER_DESCRIPTIONS, '', DG.SEMTYPE.MOLECULE);
  const toCol = createColWithDescription('string', MMP_NAMES.TO,
    mmpa.rulesBased.toFrag, FRAGMENTS_GRID_HEADER_DESCRIPTIONS, '', DG.SEMTYPE.MOLECULE);
  const occasionsCol = createColWithDescription('', MMP_NAMES.PAIRS_COUNT,
    mmpa.rulesBased.occasions, FRAGMENTS_GRID_HEADER_DESCRIPTIONS, 'int32', undefined, true);

  const fpCols = [fromCol, toCol];

  for (let i = 0; i < activityMeanNames.length; i++) {
    fpCols.push(createColWithDescription('', activityMeanNames[i],
      mmpa.rulesBased.meanDiffs[i], FRAGMENTS_GRID_HEADER_DESCRIPTIONS, 'float32', undefined, true));
  }
  fpCols.push(occasionsCol);

  const colorsPal = new Int32Array(fromCol.length);
  for (let i = 0; i < fromCol.length; i++)
    colorsPal[i] = i < fromCol.length ? 0 : 1;

  const colorCol = DG.Column.fromInt32Array(MMP_NAMES.COLOR, colorsPal);
  fpCols.push(colorCol);

  const fpDf = DG.DataFrame.fromColumns(fpCols);
  fpDf.name = MMP_NAMES.FRAGMENTS_GRID;
  const fpGrid = fpDf.plot.grid({allowColumnMenu: true});
  fpGrid.col(MMP_NAMES.COLOR)!.visible = false;

  const renderers: WeakMap<DG.Column, MmmpaActivityRenderer> = new WeakMap();
  for (let i = 0; i < activityMeanNames.length; i++) {
    const col = fpGrid.dataFrame.col(activityMeanNames[i]);
    if (col) {
      const color = activityColors!.hex[i];
      renderers.set(col, new MmmpaActivityRenderer(mmpa, i, color));
      fpGrid.col(activityMeanNames[i])!.width = MEAN_ACTIVITY_COL_WIDTH;
    }
  }

  subs.push(fpGrid.onCellRender.subscribe((e) => {
    if (e && e.bounds && e.cell && e.cell.cell?.column && !e.cell.isColHeader &&
      !e.cell.isRowHeader && renderers.has(e.cell.cell?.column)) {
      e.preventDefault();
      renderers.get(e.cell.cell?.column)?.render(e.g, e.bounds.x, e.bounds.y, e.bounds.width, e.bounds.height, e.cell);
    }
  }));

  return fpGrid;
}

function getMatchedPairsGrid(mmpa: MMPA, rdkit: RDModule) : DG.Grid {
  const pairsFromCol = createColWithDescription('string', MMP_NAMES.FROM,
    mmpa.allCasesBased.molFrom, PAIRS_GRID_HEADER_DESCRIPTIONS, '', DG.SEMTYPE.MOLECULE);
  const pairsToCol = createColWithDescription('string', MMP_NAMES.TO, mmpa.allCasesBased.molTo,
    PAIRS_GRID_HEADER_DESCRIPTIONS, '', DG.SEMTYPE.MOLECULE);
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

  const allTransformationsCols = [pairsFromCol, pairsToCol, pairNumberCol, pairNumberFromCol, pairNumberToCol,
    pairsFromSmilesCol, pairsToSmilesCol, ruleNumCol, pairNumberSortCol, coreNumCol];

  for (let i = 0; i < mmpa.initData.activitiesCount; i++) {
    const name = MMP_NAMES.DIFF + ' ' + mmpa.initData.activitiesNames[i];
    const actCol = createColWithDescription('', name, mmpa.allCasesBased.diffs[i],
      PAIRS_GRID_HEADER_DESCRIPTIONS, 'float32', undefined, false);
    allTransformationsCols.push(actCol);
    const actFromName = `${mmpa.initData.activitiesNames[i]}_from`;
    const actFromCol = createColWithDescription('', actFromName, mmpa.allCasesBased.pairedActivities[i][0],
      PAIRS_GRID_HEADER_DESCRIPTIONS, 'float32', undefined, false);
    allTransformationsCols.push(actFromCol);
    const actToName = `${mmpa.initData.activitiesNames[i]}_to`;
    const actToCol = createColWithDescription('', actToName, mmpa.allCasesBased.pairedActivities[i][1],
      PAIRS_GRID_HEADER_DESCRIPTIONS, 'float32', undefined, false);
    allTransformationsCols.push(actToCol);
  }

  const pairedTransformations = DG.DataFrame.fromColumns(allTransformationsCols);

  const substructProviderFrom = new MMPSubstructProvider(mmpa.frags.idToName,
    pairedTransformations, MMP_NAMES.FROM, rdkit, '#bc131f');
  addSubstructProvider(pairsFromCol.temp, substructProviderFrom);

  const substructProviderTo = new MMPSubstructProvider(mmpa.frags.idToName,
    pairedTransformations, MMP_NAMES.TO, rdkit, '#49bead');
  addSubstructProvider(pairsToCol.temp, substructProviderTo);

  pairedTransformations.name = MMP_NAMES.PAIRS_GRID;

  return pairedTransformations.plot.grid({allowColumnMenu: true});
}

/** Creates progress indicator only when inner function is taking more than n ms.
 */
async function withTimedProgressIndicator(pgName: string, time: number, func: () => Promise<void>) {
  let pg: DG.TaskBarProgressIndicator | null = null;
  const timer = setTimeout(() => {
    pg = DG.TaskBarProgressIndicator.create(pgName);
  });
  try {
    await func();
  } catch (e) {
    console.error(e);
  }
  clearTimeout(timer);
  if (pg)
    (pg as unknown as DG.TaskBarProgressIndicator).close();
}
