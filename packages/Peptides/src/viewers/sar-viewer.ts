import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import * as C from '../utils/constants';
import * as CR from '../utils/cell-renderer';
import {PeptidesModel, PositionStats, VIEWER_TYPE} from '../model';
import wu from 'wu';
import {SelectionItem} from '../utils/types';
import {Stats} from '../utils/statistics';
import {_package} from '../package';

export enum SELECTION_MODE {
  MUTATION_CLIFFS = 'Mutation Cliffs',
  INVARIANT_MAP = 'Invariant Map',
}

export enum MONOMER_POSITION_PROPERTIES {
  COLOR_COLUMN_NAME = 'color',
  AGGREGATION = 'aggregation',
  TARGET = 'target',
};

/** Structure-activity relationship viewer */
export class MonomerPosition extends DG.JsViewer {
  _titleHost = ui.divText(SELECTION_MODE.MUTATION_CLIFFS, {id: 'pep-viewer-title'});
  _viewerGrid!: DG.Grid;
  _model!: PeptidesModel;
  color: string;
  aggregation: string;
  target: string;
  keyPressed: boolean = false;
  currentGridCell: DG.GridCell | null = null;

  constructor() {
    super();
    this.target = this.string(MONOMER_POSITION_PROPERTIES.TARGET, null,
      {category: SELECTION_MODE.MUTATION_CLIFFS, choices: []});
    this.color = this.string(MONOMER_POSITION_PROPERTIES.COLOR_COLUMN_NAME, C.COLUMNS_NAMES.ACTIVITY_SCALED,
      {category: SELECTION_MODE.INVARIANT_MAP,
        choices: wu(grok.shell.t.columns.numerical).toArray().map((col) => col.name)});
    this.aggregation = this.string(MONOMER_POSITION_PROPERTIES.AGGREGATION, DG.AGG.AVG,
      {category: SELECTION_MODE.INVARIANT_MAP,
        choices: Object.values(DG.AGG)
          .filter((agg) => ![DG.AGG.KEY, DG.AGG.PIVOT, DG.AGG.SELECTED_ROWS_COUNT].includes(agg))});
  }

  get name(): string {return VIEWER_TYPE.MONOMER_POSITION;}

  get viewerGrid(): DG.Grid {
    if (!this._viewerGrid)
      this.createMonomerPositionGrid();
    return this._viewerGrid;
  }
  set viewerGrid(grid: DG.Grid) {
    this._viewerGrid = grid;
  }

  get mode(): SELECTION_MODE {
    return this.dataFrame.getTag(C.TAGS.MONOMER_POSITION_MODE) as SELECTION_MODE ??
      SELECTION_MODE.MUTATION_CLIFFS;
  }
  set mode(mode: SELECTION_MODE) {
    this.dataFrame.setTag(C.TAGS.MONOMER_POSITION_MODE, mode);
    this.viewerGrid.invalidate();
  }

  get model(): PeptidesModel {
    this._model ??= PeptidesModel.getInstance(this.dataFrame);
    return this._model;
  }

  detach(): void {this.subs.forEach((sub) => sub.unsubscribe());}

  onTableAttached(): void {
    super.onTableAttached();
    this.helpUrl = 'https://datagrok.ai/help/datagrok/solutions/domains/bio/peptides-sar';
    this.render();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    if (property.name === MONOMER_POSITION_PROPERTIES.TARGET)
      this.model.updateMutationCliffs().then(() => this.render(true));

    this.render(true);
  }

  createMonomerPositionDf(): DG.DataFrame {
    const uniqueMonomers = new Set<string>();
    const splitSeqCols = this.model.positionColumns.toArray();
    for (const col of splitSeqCols) {
      const colCat = col.categories;
      for (const cat of colCat) {
        if (cat !== '')
          uniqueMonomers.add(cat);
      }
    }

    const monomerCol = DG.Column.fromStrings(C.COLUMNS_NAMES.MONOMER, Array.from(uniqueMonomers));
    const monomerPositionDf = DG.DataFrame.fromColumns([monomerCol]);
    monomerPositionDf.name = 'SAR';
    for (const col of splitSeqCols)
      monomerPositionDf.columns.addNewBool(col.name);

    return monomerPositionDf;
  }


  createMonomerPositionGrid(): void {
    const monomerPositionDf = this.createMonomerPositionDf();
    this.viewerGrid = monomerPositionDf.plot.grid();
    this.viewerGrid.sort([C.COLUMNS_NAMES.MONOMER]);
    this.viewerGrid.columns.setOrder([C.COLUMNS_NAMES.MONOMER, ...this.model.positionColumns.toArray().map((col) => col.name)]);
    const monomerCol = monomerPositionDf.getCol(C.COLUMNS_NAMES.MONOMER);
    CR.setMonomerRenderer(monomerCol, this.model.alphabet);
    this.viewerGrid.onCellRender.subscribe((args: DG.GridCellRenderArgs) => renderCell(args, this.model,
      this.mode === SELECTION_MODE.INVARIANT_MAP, this.dataFrame.getCol(this.color), this.aggregation as DG.AGG));

    this.viewerGrid.onCellTooltip((gridCell: DG.GridCell, x: number, y: number) => {
      if (!gridCell.isTableCell) {
        this.model.unhighlight();
        return true;
      }
      const monomerPosition = this.getMonomerPosition(gridCell);
      this.model.highlightMonomerPosition(monomerPosition);
      return this.model.showTooltip(monomerPosition, x, y,
        {fromViewer: true, isMutationCliffs: this.mode === SELECTION_MODE.MUTATION_CLIFFS});
    });
    this.viewerGrid.root.addEventListener('mouseleave', (_ev) => this.model.unhighlight());
    DG.debounce(this.viewerGrid.onCurrentCellChanged, 500).subscribe((gridCell: DG.GridCell) => {
      try {
        if (!this.keyPressed)
          return;
        if (this.currentGridCell !== null) {
          const previousMonomerPosition = this.getMonomerPosition(this.currentGridCell);
          if (this.mode === SELECTION_MODE.INVARIANT_MAP)
            this.model.modifyInvariantMapSelection(previousMonomerPosition, {shiftPressed: true, ctrlPressed: true}, false);
          else if (this.model.mutationCliffs?.get(previousMonomerPosition.monomerOrCluster)?.get(previousMonomerPosition.positionOrClusterType)?.size)
            this.model.modifyMutationCliffsSelection(previousMonomerPosition, {shiftPressed: true, ctrlPressed: true}, false);
        }
        const monomerPosition = this.getMonomerPosition(gridCell);
        if (this.mode === SELECTION_MODE.INVARIANT_MAP)
          this.model.modifyInvariantMapSelection(monomerPosition, {shiftPressed: true, ctrlPressed: false}, true);
        else if (this.model.mutationCliffs?.get(monomerPosition.monomerOrCluster)?.get(monomerPosition.positionOrClusterType)?.size)
          this.model.modifyMutationCliffsSelection(monomerPosition, {shiftPressed: true, ctrlPressed: false}, true);

        this.viewerGrid.invalidate();
      } finally {
        this.keyPressed = false;
        this.currentGridCell = gridCell;
      }
    });
    this.viewerGrid.root.addEventListener('keydown', (ev) => {
      this.keyPressed = ev.key.startsWith('Arrow');
      if (this.keyPressed)
        return;
      if (ev.key === 'Escape' || (ev.code === 'KeyA' && ev.ctrlKey && ev.shiftKey)) {
        if (this.mode === SELECTION_MODE.INVARIANT_MAP)
          this.model.initInvariantMapSelection({notify: false});
        else
          this.model.initMutationCliffsSelection({notify: false});
      } else if (ev.code === 'KeyA' && ev.ctrlKey) {
        const positions = Object.keys(this.model.monomerPositionStats).filter((pos) => pos !== 'general');
        for (const position of positions) {
          const monomers = Object.keys(this.model.monomerPositionStats[position]!).filter((monomer) => monomer !== 'general');
          for (const monomer of monomers) {
            const monomerPosition = {monomerOrCluster: monomer, positionOrClusterType: position};
            if (this.mode === SELECTION_MODE.INVARIANT_MAP)
              this.model.modifyInvariantMapSelection(monomerPosition, {shiftPressed: true, ctrlPressed: false}, false);
            else
              this.model.modifyMutationCliffsSelection(monomerPosition, {shiftPressed: true, ctrlPressed: false}, false);
          }
        }
      }
      this.model.fireBitsetChanged();
      this.viewerGrid.invalidate();
    });
    this.viewerGrid.root.addEventListener('click', (ev) => {
      const gridCell = this.viewerGrid.hitTest(ev.offsetX, ev.offsetY);
      if (!gridCell?.isTableCell || gridCell?.tableColumn?.name === C.COLUMNS_NAMES.MONOMER)
        return;

      const monomerPosition = this.getMonomerPosition(gridCell);
      if (this.mode === SELECTION_MODE.INVARIANT_MAP) {
        this.model.modifyInvariantMapSelection(monomerPosition, {shiftPressed: ev.shiftKey, ctrlPressed: ev.ctrlKey});
        if (this.model.isInvariantMapSelectionEmpty)
          monomerPositionDf.currentRowIdx = -1;
      } else if (this.model.mutationCliffs?.get(monomerPosition.monomerOrCluster)?.get(monomerPosition.positionOrClusterType)?.size)
        this.model.modifyMutationCliffsSelection(monomerPosition, {shiftPressed: ev.shiftKey, ctrlPressed: ev.ctrlKey});

      this.viewerGrid.invalidate();

      this.showHelp();
    });

    setViewerGridProps(this.viewerGrid, false);
  }

  showHelp(): void {
    _package.files.readAsText('help/monomer-position.md').then((text) => {
      grok.shell.windows.help.showHelp(ui.markdown(text));
    }).catch((e) => grok.log.error(e));
  }

  getMonomerPosition(gridCell: DG.GridCell): SelectionItem {
    return {monomerOrCluster: gridCell.cell.dataFrame.get(C.COLUMNS_NAMES.MONOMER, gridCell!.tableRowIndex!) as string,
      positionOrClusterType: gridCell!.tableColumn!.name};
  }

  render(refreshOnly = false): void {
    if (!refreshOnly) {
      $(this.root).empty();
      let switchHost = ui.divText(VIEWER_TYPE.MOST_POTENT_RESIDUES, {id: 'pep-viewer-title'});
      if (this.name === VIEWER_TYPE.MONOMER_POSITION) {
        const mutationCliffsMode = ui.boolInput(SELECTION_MODE.MUTATION_CLIFFS, this.mode === SELECTION_MODE.MUTATION_CLIFFS);
        mutationCliffsMode.root.addEventListener('click', () => {
          invariantMapMode.value = false;
          mutationCliffsMode.value = true;
          this.mode = SELECTION_MODE.MUTATION_CLIFFS;
          this.showHelp();
        });
        mutationCliffsMode.setTooltip('Statistically significant changes in activity');
        const invariantMapMode = ui.boolInput(SELECTION_MODE.INVARIANT_MAP, this.mode === SELECTION_MODE.INVARIANT_MAP);
        invariantMapMode.root.addEventListener('click', () => {
          mutationCliffsMode.value = false;
          invariantMapMode.value = true;
          this.mode = SELECTION_MODE.INVARIANT_MAP;
          this.showHelp();
        });
        invariantMapMode.setTooltip('Number of sequences having monomer-position');
        const setDefaultProperties = (input: DG.InputBase): void => {
          $(input.root).find('.ui-input-editor').css('margin', '0px').attr('type', 'radio');
          $(input.root).find('.ui-input-description').css('padding', '0px').css('padding-right', '16px');
          $(input.captionLabel).addClass('ui-label-right');
        };
        setDefaultProperties(mutationCliffsMode);
        setDefaultProperties(invariantMapMode);

        switchHost = ui.divH([mutationCliffsMode.root, invariantMapMode.root], {id: 'pep-viewer-title'});
        $(switchHost).css('width', 'auto').css('align-self', 'center');
      }
      const viewerRoot = this.viewerGrid.root;
      viewerRoot.style.width = 'auto';
      const header = ui.divH([switchHost], {style: {alignSelf: 'center', lineHeight: 'normal'}});
      this.root.appendChild(ui.divV([header, viewerRoot]));
    }
    this.viewerGrid?.invalidate();
  }
}

/** Vertical structure activity relationship viewer */
export class MostPotentResidues extends DG.JsViewer {
  _titleHost = ui.divText(VIEWER_TYPE.MOST_POTENT_RESIDUES, {id: 'pep-viewer-title'});
  _viewerGrid!: DG.Grid;
  _model!: PeptidesModel;
  keyPressed: boolean = false;
  currentGridRowIdx: number | null = null;

  constructor() {
    super();
  }

  get name(): string {return VIEWER_TYPE.MOST_POTENT_RESIDUES;}

  get viewerGrid(): DG.Grid {
    if (!this._viewerGrid)
      this.createMostPotentResiduesGrid();
    return this._viewerGrid;
  }
  set viewerGrid(grid: DG.Grid) {
    this._viewerGrid = grid;
  }

  get model(): PeptidesModel {
    this._model ??= PeptidesModel.getInstance(this.dataFrame);
    return this._model;
  }

  detach(): void {this.subs.forEach((sub) => sub.unsubscribe());}

  onTableAttached(): void {
    super.onTableAttached();
    this.helpUrl = 'https://datagrok.ai/help/datagrok/solutions/domains/bio/peptides-sar';
    this.render();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    this.render();
  }

  createMostPotentResiduesDf(): DG.DataFrame {
    const monomerPositionStatsEntries = Object.entries(this.model.monomerPositionStats) as [string, PositionStats][];
    const posData: number[] = new Array(monomerPositionStatsEntries.length - 1);
    const monomerData: string[] = new Array(posData.length);
    const mdData: number[] = new Array(posData.length);
    const pValData: (number | null)[] = new Array(posData.length);
    const countData: number[] = new Array(posData.length);
    const ratioData: number[] = new Array(posData.length);
    const meanData: number[] = new Array(posData.length);

    let i = 0;
    for (const [position, positionStats] of monomerPositionStatsEntries) {
      const generalPositionStats = positionStats.general;
      if (!generalPositionStats)
        continue;
      if (Object.entries(positionStats).length === 1)
        continue;

      const filteredMonomerStats: [string, Stats][] = [];
      for (const [monomer, monomerStats] of Object.entries(positionStats)) {
        if (monomer === 'general')
          continue;
        if ((monomerStats as Stats).count > 1 && (monomerStats as Stats).pValue === null)
          filteredMonomerStats.push([monomer, monomerStats as Stats]);

        if ((monomerStats as Stats).pValue === generalPositionStats.minPValue)
          filteredMonomerStats.push([monomer, monomerStats as Stats]);
      }

      if (filteredMonomerStats.length === 0)
        continue;

      let maxEntry: [string, Stats] | null = null;
      for (const [monomer, monomerStats] of filteredMonomerStats) {
        if (maxEntry === null || maxEntry[1].meanDifference < monomerStats.meanDifference)
          maxEntry = [monomer, monomerStats];
      }

      if (maxEntry === null)
        continue;

      posData[i] = parseInt(position);
      monomerData[i] = maxEntry![0];
      mdData[i] = maxEntry![1].meanDifference;
      pValData[i] = maxEntry![1].pValue;
      countData[i] = maxEntry![1].count;
      ratioData[i] = maxEntry![1].ratio;
      meanData[i] = maxEntry![1].mean;
      ++i;
    }

    posData.length = i;
    monomerData.length = i;
    mdData.length = i;
    pValData.length = i;
    countData.length = i;
    ratioData.length = i;

    const mprDf = DG.DataFrame.create(i); // Subtract 'general' entry from mp-stats
    const mprDfCols = mprDf.columns;
    mprDfCols.add(DG.Column.fromList(DG.TYPE.INT, C.COLUMNS_NAMES.POSITION, posData));
    mprDfCols.add(DG.Column.fromList(DG.TYPE.STRING, C.COLUMNS_NAMES.MONOMER, monomerData));
    mprDfCols.add(DG.Column.fromList(DG.TYPE.FLOAT, C.COLUMNS_NAMES.MEAN_DIFFERENCE, mdData));
    mprDfCols.add(DG.Column.fromList(DG.TYPE.FLOAT, C.COLUMNS_NAMES.MEAN, meanData));
    mprDfCols.add(DG.Column.fromList(DG.TYPE.FLOAT, C.COLUMNS_NAMES.P_VALUE, pValData));
    mprDfCols.add(DG.Column.fromList(DG.TYPE.INT, C.COLUMNS_NAMES.COUNT, countData));
    mprDfCols.add(DG.Column.fromList(DG.TYPE.FLOAT, C.COLUMNS_NAMES.RATIO, ratioData));

    return mprDf;
  }

  createMostPotentResiduesGrid(): void {
    const mostPotentResiduesDf = this.createMostPotentResiduesDf();
    this.viewerGrid = mostPotentResiduesDf.plot.grid();
    this.viewerGrid.sort([C.COLUMNS_NAMES.POSITION]);
    const pValGridCol = this.viewerGrid.col(C.COLUMNS_NAMES.P_VALUE)!;
    pValGridCol.format = '#.000';
    pValGridCol.name = 'P-value';
    const monomerCol = mostPotentResiduesDf.getCol(C.COLUMNS_NAMES.MONOMER);

    // Setting Monomer column renderer
    CR.setMonomerRenderer(monomerCol, this.model.alphabet);
    this.viewerGrid.onCellRender.subscribe(
      (args: DG.GridCellRenderArgs) => renderCell(args, this.model, false, undefined, undefined));

    this.viewerGrid.onCellTooltip((gridCell: DG.GridCell, x: number, y: number) => {
      if (!gridCell.isTableCell) {
        this.model.unhighlight();
        return true;
      }
      const monomerPosition = this.getMonomerPosition(gridCell);
      this.model.highlightMonomerPosition(monomerPosition);
      if (gridCell.tableColumn?.name === C.COLUMNS_NAMES.MONOMER)
        monomerPosition.positionOrClusterType = C.COLUMNS_NAMES.MONOMER;
      else if (gridCell.tableColumn?.name !== C.COLUMNS_NAMES.MEAN_DIFFERENCE)
        return false;
      return this.model.showTooltip(monomerPosition, x, y, {fromViewer: true, isMutationCliffs: true});
    });
    DG.debounce(this.viewerGrid.onCurrentCellChanged, 500).subscribe((gridCell: DG.GridCell) => {
      try {
        if ((this.keyPressed && mostPotentResiduesDf.currentCol.name !== C.COLUMNS_NAMES.MEAN_DIFFERENCE) || !this.keyPressed)
          return;
        const monomerPosition = this.getMonomerPosition(gridCell);
        if (this.currentGridRowIdx !== null) {
          const previousMonomerPosition = this.getMonomerPosition(this.viewerGrid.cell('Diff', this.currentGridRowIdx));
          this.model.modifyMutationCliffsSelection(previousMonomerPosition, {shiftPressed: true, ctrlPressed: true}, false);
        }
        if (this.model.mutationCliffs?.get(monomerPosition.monomerOrCluster)?.get(monomerPosition.positionOrClusterType)?.size)
          this.model.modifyMutationCliffsSelection(monomerPosition, {shiftPressed: true, ctrlPressed: false});
        this.viewerGrid.invalidate();
      } finally {
        this.keyPressed = false;
        this.currentGridRowIdx = gridCell.gridRow;
      }
    });
    this.viewerGrid.root.addEventListener('keydown', (ev) => {
      this.keyPressed = ev.key.startsWith('Arrow');
      if (this.keyPressed)
        return;
      if (ev.key === 'Escape' || (ev.code === 'KeyA' && ev.ctrlKey && ev.shiftKey)) {
        this.model.initMutationCliffsSelection({notify: false});
      } else if (ev.code === 'KeyA' && ev.ctrlKey) {
        for (let rowIdx = 0; rowIdx < mostPotentResiduesDf.rowCount; ++rowIdx) {
          const monomerPosition = this.getMonomerPosition(this.viewerGrid.cell('Diff', rowIdx));
          this.model.modifyMutationCliffsSelection(monomerPosition, {shiftPressed: true, ctrlPressed: false}, false);
        }
      }
      this.model.fireBitsetChanged();
      this.viewerGrid.invalidate();
    });
    this.viewerGrid.root.addEventListener('mouseleave', (_ev) => this.model.unhighlight());
    this.viewerGrid.root.addEventListener('click', (ev) => {
      const gridCell = this.viewerGrid.hitTest(ev.offsetX, ev.offsetY);
      if (!gridCell?.isTableCell || gridCell!.tableColumn!.name !== C.COLUMNS_NAMES.MEAN_DIFFERENCE)
        return;

      const monomerPosition = this.getMonomerPosition(gridCell);
      if (!this.model.mutationCliffs?.get(monomerPosition.monomerOrCluster)?.get(monomerPosition.positionOrClusterType)?.size)
        return;
      this.model.modifyMutationCliffsSelection(monomerPosition, {shiftPressed: ev.shiftKey, ctrlPressed: ev.ctrlKey});
      this.viewerGrid.invalidate();

      _package.files.readAsText('help/most-potent-residues.md').then((text) => {
        grok.shell.windows.help.showHelp(ui.markdown(text));
      }).catch((e) => grok.log.error(e));
    });
    const mdCol: DG.GridColumn = this.viewerGrid.col(C.COLUMNS_NAMES.MEAN_DIFFERENCE)!;
    mdCol.name = 'Diff';
    setViewerGridProps(this.viewerGrid, true);
  }

  getMonomerPosition(gridCell: DG.GridCell): SelectionItem {
    return {monomerOrCluster: gridCell.cell.dataFrame.get(C.COLUMNS_NAMES.MONOMER, gridCell!.tableRowIndex!),
      positionOrClusterType: `${gridCell.cell.dataFrame.get(C.COLUMNS_NAMES.POSITION, gridCell!.tableRowIndex!)}`};
  }

  render(refreshOnly = false): void {
    if (!refreshOnly) {
      $(this.root).empty();
      const switchHost = ui.divText(VIEWER_TYPE.MOST_POTENT_RESIDUES, {id: 'pep-viewer-title'});
      const viewerRoot = this.viewerGrid.root;
      viewerRoot.style.width = 'auto';
      const header = ui.divH([switchHost], {style: {alignSelf: 'center', lineHeight: 'normal'}});
      this.root.appendChild(ui.divV([header, viewerRoot]));
    }
    this.viewerGrid?.invalidate();
  }
}

function renderCell(args: DG.GridCellRenderArgs, model: PeptidesModel, isInvariantMap?: boolean,
  colorCol?: DG.Column<number>, colorAgg?: DG.AGG, renderNums?: boolean): void {
  const renderColNames = [...model.positionColumns.toArray().map((col) => col.name), C.COLUMNS_NAMES.MEAN_DIFFERENCE];
  const canvasContext = args.g;
  const bound = args.bounds;

  canvasContext.save();
  canvasContext.beginPath();
  canvasContext.rect(bound.x, bound.y, bound.width, bound.height);
  canvasContext.clip();

  // Hide row column
  const cell = args.cell;
  if (cell.isRowHeader && cell.gridColumn.visible) {
    cell.gridColumn.visible = false;
    args.preventDefault();
    canvasContext.restore();
    return;
  }

  const tableColName = cell.tableColumn?.name;
  const tableRowIndex = cell.tableRowIndex!;
  if (!cell.isTableCell || renderColNames.indexOf(tableColName!) === -1) {
    canvasContext.restore();
    return;
  }

  const gridTable = cell.grid.table;
  const currentMonomer: string = gridTable.get(C.COLUMNS_NAMES.MONOMER, tableRowIndex);
  const currentPosition: string = tableColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ? tableColName :
    gridTable.get(C.COLUMNS_NAMES.POSITION, tableRowIndex).toFixed();
  const currentPosStats = model.monomerPositionStats[currentPosition];

  if (!currentPosStats![currentMonomer]) {
    args.preventDefault();
    canvasContext.restore();
    return;
  }

  if (isInvariantMap) {
    const value = currentPosStats![currentMonomer]!.count;
    const positionCol = model.df.getCol(currentPosition);
    const positionColData = positionCol.getRawData();
    const positionColCategories = positionCol.categories;

    const colorColData = colorCol!.getRawData();
    const colorValuesIndexes: number[] = [];
    for (let i = 0; i < positionCol.length; ++i) {
      if (positionColCategories[positionColData[i]] === currentMonomer)
        colorValuesIndexes.push(i);
    }
    const cellColorDataCol = DG.Column.float('color', colorValuesIndexes.length)
      .init((i) => colorColData[colorValuesIndexes[i]]);
    const colorColStats = colorCol!.stats;

    const color = DG.Color.scaleColor(cellColorDataCol.aggregate(colorAgg!), colorColStats.min, colorColStats.max);
    CR.renderInvaraintMapCell(
      canvasContext, currentMonomer, currentPosition, model.invariantMapSelection, value, bound, color);
  } else {
    CR.renderMutationCliffCell(canvasContext, currentMonomer, currentPosition, model.monomerPositionStats, bound,
      model.mutationCliffsSelection, model.mutationCliffs, renderNums);
  }
  args.preventDefault();
  canvasContext.restore();
}

function setViewerGridProps(grid: DG.Grid, isMostPotentResiduesGrid: boolean): void {
  const gridProps = grid.props;
  gridProps.allowEdit = false;
  gridProps.allowRowSelection = false;
  gridProps.allowBlockSelection = false;
  gridProps.allowColSelection = false;
  gridProps.showRowHeader = false;
  gridProps.showCurrentRowIndicator = false;

  gridProps.rowHeight = 20;
  const girdCols = grid.columns;
  const colNum = girdCols.length;
  for (let i = 0; i < colNum; ++i) {
    const col = girdCols.byIndex(i)!;
    const colName = col.name;
    col.width = isMostPotentResiduesGrid && colName !== 'Diff' && colName !== C.COLUMNS_NAMES.MONOMER ? 50 :
      gridProps.rowHeight + 10;
  }
}
