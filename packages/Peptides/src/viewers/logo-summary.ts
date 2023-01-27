import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {PeptidesModel} from '../model';
import * as C from '../utils/constants';
import * as CR from '../utils/cell-renderer';
import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {PositionHeight} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {getStats, MaskInfo, Stats} from '../utils/statistics';
import wu from 'wu';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

export class LogoSummary extends DG.JsViewer {
  _titleHost = ui.divText('Logo Summary Table', {id: 'pep-viewer-title'});
  model!: PeptidesModel;
  viewerGrid!: DG.Grid;
  initialized: boolean = false;
  webLogoMode: string;
  membersRatioThreshold: number;
  newClusterName: string;
  webLogoDfPlot: DG.DataFramePlotHelper[] = [];
  distributionDfPlot: DG.DataFramePlotHelper[] = [];

  constructor() {
    super();

    this.webLogoMode = this.string('webLogoMode', PositionHeight.full,
      {choices: [PositionHeight.full, PositionHeight.Entropy]});
    this.membersRatioThreshold = this.float('membersRatioThreshold', 0.7, {min: 0, max: 1.0});
    this.newClusterName = this.string('newClusterName', 'New cluster');
  }

  onTableAttached(): void {
    super.onTableAttached();

    this.model = PeptidesModel.getInstance(this.dataFrame);
    this.subs.push(this.model.onSettingsChanged.subscribe(() => {
      this.createLogoSummaryGrid();
      this.render();
    }));
    this.subs.push(this.model.onNewCluster.subscribe(() => this.clusterFromSelection()));
    this.subs.push(this.model.onRemoveCluster.subscribe(() => this.removeCluster()));

    this.createLogoSummaryGrid();
    this.initialized = true;
    this.render();
  }

  detach(): void {this.subs.forEach((sub) => sub.unsubscribe());}

  render(): void {
    if (this.initialized) {
      $(this.root).empty();
      const df = this.viewerGrid.dataFrame;
      if (!df.filter.anyTrue) {
        const emptyDf = ui.divText('No clusters to satisfy the threshold. ' +
          'Please, lower the threshold in viewer proeperties to include clusters');
        this.root.appendChild(ui.divV([this._titleHost, emptyDf]));
        return;
      }
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
    const activityCol = this.dataFrame.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);

    const clustersColName = this.model.settings.clustersColumnName!;
    const originalClustersCol = this.dataFrame.getCol(clustersColName);
    const originalClustersColData = originalClustersCol.getRawData();
    const originalClustersColCategories = originalClustersCol.categories;
    const originalClustersColLength = originalClustersColData.length;

    const customClustersColumnsList = wu(this.model.customClusters).toArray();
    const getAggregatedColName = (aggF: string, colName: string) => `${aggF}(${colName})`;

    let summaryTableBuilder = this.dataFrame.groupBy([clustersColName]);
    const aggregateColumnsEntries = Object.entries(this.model.settings.columns ?? {});
    for (const [colName, aggregationFunc] of aggregateColumnsEntries) {
      summaryTableBuilder = summaryTableBuilder.add(
        aggregationFunc as any, colName, getAggregatedColName(aggregationFunc, colName));
    }

    const tempSummaryTable = summaryTableBuilder.aggregate();
    const tempSummaryTableLength = tempSummaryTable.rowCount;
    const tempClustersCol: DG.Column<string> = tempSummaryTable.getCol(clustersColName);
    const summaryTableLength = tempSummaryTableLength + customClustersColumnsList.length;
    const summaryTable = DG.DataFrame.create(summaryTableLength);
    const summaryTableCols = summaryTable.columns;

    const clustersCol = summaryTableCols.addNewString(clustersColName);
    for (let i = 0; i < summaryTableLength; ++i) {
      clustersCol.set(i, i < tempSummaryTableLength ? tempClustersCol.get(i) :
        customClustersColumnsList[i - tempSummaryTableLength].name);
    }
    const clustersColData = clustersCol.getRawData();
    const clustersColCategories = clustersCol.categories;

    const peptideCol: DG.Column<string> = this.dataFrame.getCol(this.model.settings.sequenceColumnName!);
    const peptideColData = peptideCol.getRawData();
    const peptideColCategories = peptideCol.categories;
    const peptideColTags = peptideCol.tags;

    const membersColData = summaryTableCols.addNewInt(C.LST_COLUMN_NAMES.MEMBERS).getRawData();
    const webLogoCol = summaryTableCols.addNewString(C.LST_COLUMN_NAMES.WEB_LOGO);
    const distributionCol = summaryTableCols.addNewString(C.LST_COLUMN_NAMES.DISTRIBUTION);
    const meanDifferenceColData = summaryTableCols.addNewFloat(C.LST_COLUMN_NAMES.MEAN_DIFFERENCE).getRawData();
    const pValColData = summaryTableCols.addNewFloat(C.LST_COLUMN_NAMES.P_VALUE).getRawData();
    const ratioColData = summaryTableCols.addNewFloat(C.LST_COLUMN_NAMES.RATIO).getRawData();

    for (const [colName, aggregationFunc] of aggregateColumnsEntries) {
      const tempSummaryTableCol = tempSummaryTable.getCol(getAggregatedColName(aggregationFunc, colName));
      const summaryTableCol = summaryTableCols.addNew(tempSummaryTableCol.name, tempSummaryTableCol.type);
      summaryTableCol.init((i) => i < tempSummaryTableLength ? tempSummaryTableCol.get(i) : null);
    }

    this.webLogoDfPlot = new Array(summaryTableLength);
    this.distributionDfPlot = new Array(summaryTableLength);

    for (let summaryTableRowIndex = 0; summaryTableRowIndex < summaryTableLength; ++summaryTableRowIndex) {
      const isOriginalCluster = summaryTableRowIndex < tempSummaryTableLength;
      const currentClusterCategoryIndex = clustersColData[summaryTableRowIndex];
      const currentCluster = clustersColCategories[currentClusterCategoryIndex]; // Cluster name
      const customClusterColData = customClustersColumnsList.find((col) => col.name == currentCluster)?.toList();

      const isValidIndex = isOriginalCluster ?
        (j: number) => originalClustersColCategories[originalClustersColData[j]] == currentCluster :
        (j: number) => customClusterColData![j];

      const stats = this.model.clusterStats[currentClusterCategoryIndex];
      const tCol = DG.Column.string('peptides', stats.count);
      let tColIdx = 0;
      for (let j = 0; j < originalClustersColLength; ++j) {
        if (isValidIndex(j))
          tCol.set(tColIdx++, peptideColCategories[peptideColData[j]]);
      }

      for (const tag of peptideColTags)
        tCol.setTag(tag[0], tag[1]);

      const uh = new UnitsHandler(tCol);
      tCol.setTag(bioTAGS.alphabetSize, uh.getAlphabetSize().toString());

      //TODO: use bitset instead of splitCol
      const splitCol = DG.Column.bool(C.COLUMNS_NAMES.SPLIT_COL, activityCol.length);
      const getSplitColValueAt = isOriginalCluster ?
        (splitColIndex: number) => originalClustersColData[splitColIndex] == currentClusterCategoryIndex :
        (splitColIndex: number) => customClusterColData![splitColIndex];
      splitCol.init((i) => getSplitColValueAt(i));

      const distributionTable = DG.DataFrame.fromColumns([activityCol, splitCol]);
      const dfSlice = DG.DataFrame.fromColumns([tCol]);

      this.webLogoDfPlot[summaryTableRowIndex] = dfSlice.plot;
      this.distributionDfPlot[summaryTableRowIndex] = distributionTable.plot;

      membersColData[summaryTableRowIndex] = stats.count;
      meanDifferenceColData[summaryTableRowIndex] = stats.meanDifference;
      pValColData[summaryTableRowIndex] = stats.pValue;
      ratioColData[summaryTableRowIndex] = stats.ratio;

      //Setting aggregated col values
      if (!isOriginalCluster) {
        for (const [colName, aggregationFunc] of aggregateColumnsEntries) {
          const arrayBuffer = this.dataFrame.getCol(colName).getRawData();
          const clusterMask = DG.BitSet.fromBytes(arrayBuffer.buffer, arrayBuffer.byteLength / 4);
          const subDf = this.dataFrame.clone(clusterMask, [colName]);
          const newColName = getAggregatedColName(aggregationFunc, colName);
          const aggregatedDf = subDf.groupBy()
            .add(aggregationFunc as any, colName, newColName)
            .aggregate();
          const value = aggregatedDf.get(newColName, 0);
          summaryTable.set(newColName, summaryTableRowIndex, value);
        }
      }
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
        this.webLogoDfPlot[currentRowIdx]
          .fromType('WebLogo', {maxHeight: cell.grid.props.rowHeight - 5, positionHeight: this.webLogoMode,
            horizontalAlignment: 'left'})
          .then((viewer) => cell.element = viewer.root);
      } else if (cell.tableColumn?.name == 'Distribution') {
        const viewerRoot = this.distributionDfPlot[currentRowIdx].histogram({
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

      summaryTable.currentRowIdx = -1;
      if (ev.shiftKey)
        this.model.modifyClusterSelection(cell.tableRowIndex!);
      else
        this.model.initClusterSelection(cell.tableRowIndex!);
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
      CR.renderLogoSummaryCell(canvasContext, gc.cell.value, gc.cell.rowIndex, this.model.logoSummarySelection, bound);
      gridCellArgs.preventDefault();
      canvasContext.restore();
    });
    this.viewerGrid.onCellTooltip((cell, x, y) => {
      if (!cell.isColHeader && cell.tableColumn?.name === clustersColName)
        this.model.showTooltipCluster(cell.cell.rowIndex, x, y, cell.cell.value);
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
    const memberstCol = table.getCol(C.LST_COLUMN_NAMES.MEMBERS);
    const membersColData = memberstCol.getRawData();
    const maxCount = memberstCol.stats.max;
    const minMembers = Math.ceil(maxCount * this.membersRatioThreshold);
    table.filter.init((i) => membersColData[i] > minMembers);
  }

  clusterFromSelection(): void {
    const selection = this.dataFrame.selection;
    const viewerDf = this.viewerGrid.dataFrame;
    const viewerDfCols = viewerDf.columns;
    const viewerDfColsLength = viewerDfCols.length;
    const newClusterVals = new Array(viewerDfCols.length);

    const activityScaledCol = this.dataFrame.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    const maskInfo: MaskInfo = {
      mask: selection.getBuffer(),
      trueCount: selection.trueCount,
      falseCount: selection.falseCount,
    };
    const stats = getStats(activityScaledCol.getRawData(), maskInfo);
    const distributionTable = DG.DataFrame.fromColumns([activityScaledCol, this.model.splitCol]);

    const peptideCol: DG.Column<string> = this.dataFrame.getCol(this.model.settings.sequenceColumnName!);
    const peptideColData = peptideCol.getRawData();
    const peptideColCategories = peptideCol.categories;
    const peptideColTags = peptideCol.tags;
    const selectedIndexes = selection.getSelectedIndexes();
    const tCol = DG.Column.string('peptides', selectedIndexes.length);

    for (let i = 0; i < selectedIndexes.length; ++i)
      tCol.set(i, peptideColCategories[peptideColData[selectedIndexes[i]]]);
    for (const tag of peptideColTags)
      tCol.setTag(tag[0], tag[1]);

    const uh = new UnitsHandler(tCol);
    tCol.setTag(bioTAGS.alphabetSize, uh.getAlphabetSize().toString());

    const webLogoTable = DG.DataFrame.fromColumns([tCol]);
    this.webLogoDfPlot.push(webLogoTable.plot);
    this.distributionDfPlot.push(distributionTable.plot);

    const colCategories = viewerDfCols.byName(this.model.settings.clustersColumnName!).categories;
    let newClusterName = this.newClusterName;
    let clusterNum = 1;
    const getString = !isNaN(parseInt(newClusterName)) ? () => `${parseInt(newClusterName) + 1}` :
      newClusterName == '' ? () => `${clusterNum++}` :
        () => `${this.newClusterName} ${clusterNum++}`;
    while (colCategories.includes(newClusterName))
      newClusterName = getString();

    this.getProperty('newClusterName')?.set(this, getString());

    for (let i = 0; i < viewerDfColsLength; ++i) {
      const col = viewerDfCols.byIndex(i);
      newClusterVals[i] = col.name == this.model.settings.clustersColumnName! ? newClusterName :
        col.name == C.LST_COLUMN_NAMES.MEMBERS ? maskInfo.trueCount :
          col.name == C.LST_COLUMN_NAMES.WEB_LOGO ? null :
            col.name == C.LST_COLUMN_NAMES.DISTRIBUTION ? null :
              col.name == C.LST_COLUMN_NAMES.MEAN_DIFFERENCE ? stats.meanDifference:
                col.name == C.LST_COLUMN_NAMES.P_VALUE ? stats.pValue:
                  col.name == C.LST_COLUMN_NAMES.RATIO ? stats.ratio:
        console.warn(`PeptidesLSTWarn: value for column ${col.name} is undefined`)! || null;
    }
    viewerDf.rows.addNew(newClusterVals);

    this.model.clusterStats.push(stats);
    this.model.addNewCluster(newClusterName);
  }

  removeCluster(): void {
    const lss = this.model.logoSummarySelection;

    const viewerDf = this.viewerGrid.dataFrame;
    const viewerClustersCol = viewerDf.getCol(this.model.settings.clustersColumnName!);
    const originalClusterAmount = this.dataFrame.getCol(this.model.settings.clustersColumnName!).categories.length;

    const removeClusterIndexesList = lss.filter((clusterId) => clusterId + 1 > originalClusterAmount);
    if (removeClusterIndexesList.length == 0)
      return grok.shell.info('Nothing removed. Please select a created cluster to remove');

    for (const clusterId of removeClusterIndexesList) {
      lss.splice(lss.indexOf(clusterId), 1);
      this.dataFrame.columns.remove(viewerClustersCol.get(clusterId));
    }

    this.model.logoSummarySelection = lss;
    this.model.clusterStats = this.model.calculateClusterStatistics();
    this.createLogoSummaryGrid();
    this.render();
  }
}
