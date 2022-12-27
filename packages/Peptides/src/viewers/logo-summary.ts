import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {PeptidesModel} from '../model';
import * as C from '../utils/constants';
import * as CR from '../utils/cell-renderer';
import {Stats} from '../utils/statistics';
import {PositionHeight, UnitsHandler, TAGS as bioTAGS} from '@datagrok-libraries/bio';

export class LogoSummary extends DG.JsViewer {
  _titleHost = ui.divText('Logo Summary Table', {id: 'pep-viewer-title'});
  model!: PeptidesModel;
  viewerGrid!: DG.Grid;
  initialized: boolean = false;
  webLogoMode: string;
  membersRatioThreshold: number;

  constructor() {
    super();

    this.webLogoMode = this.string('webLogoMode', PositionHeight.full,
      {choices: [PositionHeight.full, PositionHeight.Entropy]});
    this.membersRatioThreshold = this.float('membersRatioThreshold', 0.7, {min: 0.01, max: 1.0});
  }

  onTableAttached(): void {
    super.onTableAttached();

    this.model = PeptidesModel.getInstance(this.dataFrame);
    this.subs.push(this.model.onSettingsChanged.subscribe(() => {
      this.createLogoSummaryGrid();
      this.render();
    }));

    this.createLogoSummaryGrid();
    this.initialized = true;
    this.render();
  }

  detach(): void {this.subs.forEach((sub) => sub.unsubscribe());}

  render(): void {
    if (this.initialized) {
      $(this.root).empty();
      this.viewerGrid.root.style.width = 'auto';
      this.root.appendChild(ui.divV([this._titleHost, this.viewerGrid.root]));
      this.viewerGrid.invalidate();
    }
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    if (property.name == 'membersRatioThreshold')
      this.updateFilter();
    this.render();
  }

  createLogoSummaryGrid(): DG.Grid {
    const clustersColName = this.model.settings.clustersColumnName!;
    let summaryTableBuilder = this.dataFrame.groupBy([clustersColName]);
    for (const [colName, aggregationFunc] of Object.entries(this.model.settings.columns ?? {}))
      summaryTableBuilder = summaryTableBuilder.add(aggregationFunc as any, colName, `${aggregationFunc}(${colName})`);

    const summaryTable = summaryTableBuilder.aggregate();
    const summaryTableCols = summaryTable.columns;
    const webLogoCol: DG.Column<string> = summaryTableCols.addNewString('WebLogo');
    const clustersCol: DG.Column<string> = summaryTable.getCol(clustersColName);
    const clustersColData = clustersCol.getRawData();
    const clustersColCategories = clustersCol.categories;
    const summaryTableLength = clustersColData.length;
    const membersColData = summaryTableCols.addNewInt(C.COLUMNS_NAMES.MEMBERS).getRawData();
    const tempWebLogoDfPlotList: DG.DataFramePlotHelper[] = new Array(summaryTableLength);
    const tempDistributionDfPlotList: DG.DataFramePlotHelper[] = new Array(summaryTableLength);
    const originalClustersCol = this.dataFrame.getCol(clustersColName);
    const originalClustersColData = originalClustersCol.getRawData();
    const originalClustersColCategories = originalClustersCol.categories;
    const originalClustersColLength = originalClustersColData.length;
    const peptideCol: DG.Column<string> = this.dataFrame.getCol(this.model.settings.sequenceColumnName!);
    const peptideColData = peptideCol.getRawData();
    const peptideColCategories = peptideCol.categories;
    const peptideColTags = peptideCol.tags;
    const meanDifferenceColData = summaryTableCols.addNewFloat('Mean difference').getRawData();
    const pValColData = summaryTableCols.addNewFloat('P-Value').getRawData();
    const ratioColData = summaryTableCols.addNewFloat('Ratio').getRawData();
    const distributionCol = summaryTableCols.addNewString('Distribution');
    const activityCol = this.dataFrame.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);

    for (let summaryTableRowIndex = 0; summaryTableRowIndex < summaryTableLength; ++summaryTableRowIndex) {
      const indexes: number[] = [];
      const currentClusterCategoryIndex = clustersColData[summaryTableRowIndex];
      const currentCluster = clustersColCategories[currentClusterCategoryIndex];
      for (let j = 0; j < originalClustersColLength; ++j) {
        if (originalClustersColCategories[originalClustersColData[j]] == currentCluster)
          indexes.push(j);
      }
      const tCol = DG.Column.string('peptides', indexes.length);
      tCol.init((i) => peptideColCategories[peptideColData[indexes[i]]]);

      for (const tag of peptideColTags)
        tCol.setTag(tag[0], tag[1]);

      const uh = new UnitsHandler(tCol);
      tCol.setTag(bioTAGS.alphabetSize, uh.getAlphabetSize().toString());

      //TODO: use bitset instead of splitCol
      const splitCol = DG.Column.bool(C.COLUMNS_NAMES.SPLIT_COL, activityCol.length);
      splitCol.init((splitColIndex) => originalClustersColData[splitColIndex] == currentClusterCategoryIndex);
      const distributionTable = DG.DataFrame.fromColumns([activityCol, splitCol]);

      const stats = this.model.clusterStats[currentClusterCategoryIndex];
      const dfSlice = DG.DataFrame.fromColumns([tCol]);
      tempWebLogoDfPlotList[summaryTableRowIndex] = dfSlice.plot;
      tempDistributionDfPlotList[summaryTableRowIndex] = distributionTable.plot;
      membersColData[summaryTableRowIndex] = stats.count;
      meanDifferenceColData[summaryTableRowIndex] = stats.meanDifference;
      pValColData[summaryTableRowIndex] = stats.pValue;
      ratioColData[summaryTableRowIndex] = stats.ratio;
    }
    webLogoCol.setTag(DG.TAGS.CELL_RENDERER, 'html');
    distributionCol.setTag(DG.TAGS.CELL_RENDERER, 'html');

    this.viewerGrid = summaryTable.plot.grid();
    this.updateFilter();
    const gridClustersCol = this.viewerGrid.col(clustersColName)!;
    gridClustersCol.name = 'Clusters';
    gridClustersCol.visible = true;
    this.viewerGrid.columns.rowHeader!.visible = false;
    this.viewerGrid.props.rowHeight = 55;
    this.viewerGrid.onCellPrepare((cell) => {
      const currentRowIdx = cell.tableRowIndex;
      if (!cell.isTableCell || currentRowIdx == null || currentRowIdx == -1)
        return;

      if (cell.tableColumn?.name == 'WebLogo') {
        tempWebLogoDfPlotList[currentRowIdx]
          .fromType('WebLogo', {maxHeight: cell.grid.props.rowHeight - 5, positionHeight: this.webLogoMode})
          .then((viewer) => cell.element = viewer.root);
      } else if (cell.tableColumn?.name == 'Distribution') {
        const viewerRoot = tempDistributionDfPlotList[currentRowIdx].histogram({
          filteringEnabled: false,
          valueColumnName: C.COLUMNS_NAMES.ACTIVITY_SCALED,
          splitColumnName: C.COLUMNS_NAMES.SPLIT_COL,
          legendVisibility: 'Never',
          showXAxis: true,
          showColumnSelector: false,
          showRangeSlider: false,
          showBinSelector: false,
          backColor: '#fffff',
        }).root;

        viewerRoot.style.width = 'auto';
        const height = (cell.grid.props.rowHeight - 5) / 2 * 3;
        viewerRoot.style.height = `${height}px`;
        cell.element = viewerRoot;
      }
    });
    this.viewerGrid.root.addEventListener('click', (ev) => {
      const cell = this.viewerGrid.hitTest(ev.offsetX, ev.offsetY);
      if (!cell || !cell.isTableCell)
        return;

      const clusterIdx = clustersColData[cell.tableRowIndex!];
      summaryTable.currentRowIdx = -1;
      if (ev.shiftKey)
        this.model.modifyClusterSelection(clusterIdx);
      else
        this.model.initClusterSelection(clusterIdx);
      this.viewerGrid.invalidate();
    });
    this.viewerGrid.onCellRender.subscribe((gridCellArgs) => {
      const gc = gridCellArgs.cell;
      if (gc.tableColumn?.name !== clustersColName || gc.isColHeader)
        return;
      const canvasContext = gridCellArgs.g;
      const bound = gridCellArgs.bounds;
      canvasContext.save();
      canvasContext.beginPath();
      canvasContext.rect(bound.x, bound.y, bound.width, bound.height);
      canvasContext.clip();
      const cellRawData = clustersColData[gc.cell.rowIndex];
      CR.renderLogoSummaryCell(canvasContext, gc.cell.value, cellRawData, this.model.logoSummarySelection, bound);
      gridCellArgs.preventDefault();
      canvasContext.restore();
    });
    this.viewerGrid.onCellTooltip((cell, x, y) => {
      if (!cell.isColHeader && cell.tableColumn?.name === clustersColName)
        this.model.showTooltipCluster(clustersColData[cell.cell.rowIndex], x, y);
      return true;
    });
    const webLogoGridCol = this.viewerGrid.columns.byName('WebLogo')!;
    webLogoGridCol.cellType = 'html';
    webLogoGridCol.width = 350;

    const gridProps = this.viewerGrid.props;
    gridProps.allowEdit = false;
    gridProps.allowRowSelection = false;
    gridProps.allowBlockSelection = false;
    gridProps.allowColSelection = false;

    return this.viewerGrid;
  }

  updateFilter(): void {
    const table = this.viewerGrid.table;
    const memberstCol = table.getCol(C.COLUMNS_NAMES.MEMBERS);
    const membersColData = memberstCol.getRawData();
    const maxCount = memberstCol.stats.max;
    table.filter.init((i) => membersColData[i] > Math.ceil(maxCount * this.membersRatioThreshold));
  }
}
