import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import * as C from '../utils/constants';
import * as CR from '../utils/cell-renderer';
import {HorizontalAlignments, IWebLogoViewer, PositionHeight} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {ClusterTypeStats, getAggregatedColumnValues, getAggregatedValue, getStats, Stats} from '../utils/statistics';
import wu from 'wu';
import {getActivityDistribution, getStatsTableMap} from '../widgets/distribution';
import {getDistributionPanel, getDistributionTable} from '../utils/misc';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import * as type from '../utils/types';
import {SelectionItem} from '../utils/types';
import {_package} from '../package';
import {calculateClusterStatistics} from "../utils/algorithms";

const getAggregatedColName = (aggF: string, colName: string): string => `${aggF}(${colName})`;

export enum CLUSTER_TYPE {
  ORIGINAL = 'original',
  CUSTOM = 'custom',
}

export type ClusterType = `${CLUSTER_TYPE}`;

export const enum LST_PROPERTIES {
  WEB_LOGO_MODE = 'webLogoMode',
  MEMBERS_RATIO_THRESHOLD = 'membersRatioThreshold',
  SEQUENCE_COLUMN_NAME = 'sequence',
  CLUSTERS_COLUMN_NAME = 'clusters',
}

enum LST_CATEGORIES {
  GENERAL = 'General',
  STYLE = 'Style',
}

export class LogoSummaryTable extends DG.JsViewer {
  _titleHost = ui.divText(VIEWER_TYPE.LOGO_SUMMARY_TABLE, {id: 'pep-viewer-title'});
  sequenceColumnName: string;
  clustersColumnName: string;
  initialized: boolean = false;
  webLogoMode: string;
  membersRatioThreshold: number;
  bitsets: DG.BitSet[] = [];
  keyPress: boolean = false;
  currentRowIndex: number | null = null;

  constructor() {
    super();

    this.sequenceColumnName = this.column(LST_PROPERTIES.SEQUENCE_COLUMN_NAME, {category: LST_CATEGORIES.GENERAL, nullable: false});
    this.clustersColumnName = this.column(LST_PROPERTIES.CLUSTERS_COLUMN_NAME, {category: LST_CATEGORIES.GENERAL, nullable: false});

    this.webLogoMode = this.string(LST_PROPERTIES.WEB_LOGO_MODE, PositionHeight.Entropy,
      {choices: [PositionHeight.full, PositionHeight.Entropy], category: LST_CATEGORIES.STYLE});
    this.membersRatioThreshold = this.float(LST_PROPERTIES.MEMBERS_RATIO_THRESHOLD, 0.3,
      {min: 0, max: 1.0, category: LST_CATEGORIES.STYLE});
  }

  _model!: PeptidesModel;

  get model(): PeptidesModel {
    this._model ??= PeptidesModel.getInstance(this.dataFrame);
    return this._model;
  }

  _viewerGrid!: DG.Grid;

  get viewerGrid(): DG.Grid {
    this._viewerGrid ??= this.createLogoSummaryTableGrid();
    this.initialized = true;
    return this._viewerGrid;
  }

  set viewerGrid(grid: DG.Grid) {
    this._viewerGrid = grid;
    this.render();
  }

  _clusterStats?: ClusterTypeStats;

  get clusterStats(): ClusterTypeStats {
    this._clusterStats ??= calculateClusterStatistics(this.dataFrame, this.clustersColumnName, this.customClusters.toArray());
    return this._clusterStats;
  }

  set clusterStats(clusterStats: ClusterTypeStats) {
    this._clusterStats = clusterStats;
    this.viewerGrid = this.createLogoSummaryTableGrid();
  }

  _clusterSelection?: type.Selection;

  get clusterSelection(): type.Selection {
    const tagSelection = this.dataFrame.getTag(C.TAGS.CLUSTER_SELECTION);
    this._clusterSelection ??= tagSelection === null ? this.initClusterSelection({notify: false}) :
      JSON.parse(tagSelection);
    return this._clusterSelection!;
  }

  set clusterSelection(selection: type.Selection) {
    this._clusterSelection = selection;
    this.dataFrame.tags[C.TAGS.CLUSTER_SELECTION] = JSON.stringify(selection);
    this.model.fireBitsetChanged();
    this.model.analysisView.grid.invalidate();
  }

  _logoSummaryTable?: DG.DataFrame;

  get logoSummaryTable(): DG.DataFrame {
    this._logoSummaryTable ??= this.createLogoSummaryTable();
    return this._logoSummaryTable;
  }

  set logoSummaryTable(df: DG.DataFrame) {
    this._logoSummaryTable = df;
  }

  get isClusterSelectionEmpty(): boolean {
    return (this.clusterSelection[CLUSTER_TYPE.ORIGINAL].length + this.clusterSelection[CLUSTER_TYPE.CUSTOM].length) === 0;
  }

  get customClusters(): wu.WuIterable<DG.Column<boolean>> {
    const query: { [key: string]: string } = {};
    query[C.TAGS.CUSTOM_CLUSTER] = '1';
    return wu(this.dataFrame.columns.byTags(query));
  }

  onTableAttached(): void {
    super.onTableAttached();
    this.render();
  }

  detach(): void {this.subs.forEach((sub) => sub.unsubscribe());}

  render(): void {
    $(this.root).empty();
    if (this.clustersColumnName == null || this.sequenceColumnName == null) {
      this.root.appendChild(ui.divText('Please, select a sequence and cluster columns in the viewer properties'));
      return;
    }
    if (!this.initialized) {
      return;
    }

    if (!this.logoSummaryTable.filter.anyTrue) {
      const emptyDf = ui.divText('No clusters to satisfy the threshold. ' +
        'Please, lower the threshold in viewer proeperties to include clusters');
      this.root.appendChild(ui.divV([this._titleHost, emptyDf]));
      return;
    }
    const expand = ui.iconFA('expand-alt', () => {
      const dialog = ui.dialog('Logo Summary Table');
      dialog.add(this.viewerGrid.root);
      dialog.onCancel(() => this.render());
      dialog.showModal(true);
      this.viewerGrid.invalidate();
    }, 'Show Logo Summary Table in full screen');
    $(expand).addClass('pep-help-icon');
    this.viewerGrid.root.style.width = 'auto';
    this.root.appendChild(ui.divV([
      ui.divH([this._titleHost, expand], {style: {alignSelf: 'center', lineHeight: 'normal'}}),
      this.viewerGrid.root,
    ]));
    this.viewerGrid.invalidate();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    if (property.name === 'membersRatioThreshold') {
      this.updateFilter();
    } else if (property.name === 'sequenceColumnName' || property.name === 'clustersColumnName') {
      // this.model.settings = {sequenceColumnName: this.sequenceColName, clustersColumnName: this.clustersColumnName};
      if (this.sequenceColumnName !== null && this.clustersColumnName !== null)
        this.viewerGrid = this.createLogoSummaryTableGrid();
    }
    // this.render();
  }

  initClusterSelection(options: { notify?: boolean } = {}): type.Selection {
    options.notify ??= true;

    const newClusterSelection = {} as type.Selection;
    newClusterSelection[CLUSTER_TYPE.ORIGINAL] = [];
    newClusterSelection[CLUSTER_TYPE.CUSTOM] = [];
    if (options.notify)
      this.clusterSelection = newClusterSelection;
    else
      this._clusterSelection = newClusterSelection;

    return this.clusterSelection;
  }

  createLogoSummaryTable(): DG.DataFrame {
    const clustersColName = this.clustersColumnName;
    const isDfFiltered = this.dataFrame.filter.anyFalse;
    const filteredDf = isDfFiltered ? this.dataFrame.clone(this.dataFrame.filter) : this.dataFrame;
    const filteredDfCols = filteredDf.columns;
    const filteredDfRowCount = filteredDf.rowCount;
    const activityCol = filteredDf.getCol(C.COLUMNS_NAMES.ACTIVITY);
    const activityColData = activityCol.getRawData();

    const filteredDfClustCol = filteredDf.getCol(clustersColName);
    const filteredDfClustColData = filteredDfClustCol.getRawData();
    const filteredDfClustColCat = filteredDfClustCol.categories;

    const query: { [key: string]: string } = {};
    query[C.TAGS.CUSTOM_CLUSTER] = '1';
    const customClustColList: DG.Column<boolean>[] = wu(filteredDfCols.byTags(query)).filter((c) => c.max > 0).toArray();

    const customLST = DG.DataFrame.create(customClustColList.length);
    const customLSTCols = customLST.columns;
    const customLSTClustCol = customLSTCols.addNewString(C.LST_COLUMN_NAMES.CLUSTER);

    const customMembersColData = customLSTCols.addNewInt(C.LST_COLUMN_NAMES.MEMBERS).getRawData();
    const customWebLogoCol = customLSTCols.addNewString(C.LST_COLUMN_NAMES.WEB_LOGO);
    const customDistCol = customLSTCols.addNewString(C.LST_COLUMN_NAMES.DISTRIBUTION);
    const customMDColData = customLSTCols.addNewFloat(C.LST_COLUMN_NAMES.MEAN_DIFFERENCE).getRawData();
    const customPValColData = customLSTCols.addNewFloat(C.LST_COLUMN_NAMES.P_VALUE).getRawData();
    const customRatioColData = customLSTCols.addNewFloat(C.LST_COLUMN_NAMES.RATIO).getRawData();

    let origLSTBuilder = filteredDf.groupBy([clustersColName]);
    const aggColsEntries = Object.entries(this.model.settings.columns ?? {});
    const aggColNames = aggColsEntries.map(([colName, aggFn]) => getAggregatedColName(aggFn, colName));
    const customAggRawCols = new Array(aggColNames.length);
    const colAggEntries = aggColsEntries.map(
      ([colName, aggFn]) => [filteredDf.getCol(colName), aggFn] as [DG.Column<number>, DG.AggregationType]);

    for (let aggIdx = 0; aggIdx < aggColsEntries.length; ++aggIdx) {
      const [colName, aggFn] = aggColsEntries[aggIdx];
      origLSTBuilder = origLSTBuilder.add(aggFn, colName, aggColNames[aggIdx]);
      const customLSTAggCol = customLSTCols.addNewFloat(aggColNames[aggIdx]);
      customAggRawCols[aggIdx] = customLSTAggCol.getRawData();
    }

    // BEGIN: fill LST part with custom clusters
    const customBitsets: DG.BitSet[] = new Array(customClustColList.length);

    for (let rowIdx = 0; rowIdx < customClustColList.length; ++rowIdx) {
      const customClustCol = customClustColList[rowIdx];
      customLSTClustCol.set(rowIdx, customClustCol.name);
      const bitArray = BitArray.fromUint32Array(filteredDfRowCount, customClustCol.getRawData() as Uint32Array);
      if (bitArray.allFalse || bitArray.allTrue)
        continue;
      const bsMask = DG.BitSet.fromBytes(bitArray.buffer.buffer, filteredDfRowCount);

      const stats: Stats = isDfFiltered ? getStats(activityColData, bitArray) :
        this.clusterStats[CLUSTER_TYPE.CUSTOM][customClustCol.name];

      customMembersColData[rowIdx] = stats.count;
      customBitsets[rowIdx] = bsMask;
      customMDColData[rowIdx] = stats.meanDifference;
      customPValColData[rowIdx] = stats.pValue ?? DG.FLOAT_NULL;
      customRatioColData[rowIdx] = stats.ratio;

      for (let aggColIdx = 0; aggColIdx < aggColNames.length; ++aggColIdx) {
        const [col, aggFn] = colAggEntries[aggColIdx];
        customAggRawCols[aggColIdx][rowIdx] = getAggregatedValue(col, aggFn, bsMask);
      }
    }

    customWebLogoCol.setTag(DG.TAGS.CELL_RENDERER, 'html');
    customDistCol.setTag(DG.TAGS.CELL_RENDERER, 'html');

    // END

    // BEGIN: fill LST part with original clusters
    const origLST = origLSTBuilder.aggregate();
    const origLSTLen = origLST.rowCount;
    const origLSTCols = origLST.columns;
    const origLSTClustCol: DG.Column<string> = origLST.getCol(clustersColName);
    origLSTClustCol.name = C.LST_COLUMN_NAMES.CLUSTER;

    const origLSTClustColCat = origLSTClustCol.categories;

    const origMembersColData = origLSTCols.addNewInt(C.LST_COLUMN_NAMES.MEMBERS).getRawData();
    const origWebLogoCol = origLSTCols.addNewString(C.LST_COLUMN_NAMES.WEB_LOGO);
    const origDistCol = origLSTCols.addNewString(C.LST_COLUMN_NAMES.DISTRIBUTION);
    const origMDColData = origLSTCols.addNewFloat(C.LST_COLUMN_NAMES.MEAN_DIFFERENCE).getRawData();
    const origPValColData = origLSTCols.addNewFloat(C.LST_COLUMN_NAMES.P_VALUE).getRawData();
    const origRatioColData = origLSTCols.addNewFloat(C.LST_COLUMN_NAMES.RATIO).getRawData();
    const origBitsets: DG.BitSet[] = new Array(origLSTLen);

    const origClustMasks = Array.from({length: origLSTLen},
      () => new BitArray(filteredDfRowCount, false));

    for (let rowIdx = 0; rowIdx < filteredDfRowCount; ++rowIdx) {
      const filteredClustName = filteredDfClustColCat[filteredDfClustColData[rowIdx]];
      const origClustIdx = origLSTClustColCat.indexOf(filteredClustName);
      origClustMasks[origClustIdx].setTrue(rowIdx);
    }

    for (let rowIdx = 0; rowIdx < origLSTLen; ++rowIdx) {
      const mask = origClustMasks[rowIdx];
      if (mask.allFalse || mask.allTrue)
        continue;
      const bsMask = DG.BitSet.fromBytes(mask.buffer.buffer, filteredDfRowCount);

      const stats = isDfFiltered ? getStats(activityColData, mask) :
        this.clusterStats[CLUSTER_TYPE.ORIGINAL][origLSTClustColCat[rowIdx]];

      origMembersColData[rowIdx] = stats.count;
      origBitsets[rowIdx] = bsMask;
      origMDColData[rowIdx] = stats.meanDifference;
      origPValColData[rowIdx] = stats.pValue ?? DG.FLOAT_NULL;
      origRatioColData[rowIdx] = stats.ratio;
    }

    origWebLogoCol.setTag(DG.TAGS.CELL_RENDERER, 'html');
    origDistCol.setTag(DG.TAGS.CELL_RENDERER, 'html');
    // END

    // combine LSTs and create a grid
    const summaryTable = origLST.append(customLST);
    this.bitsets = origBitsets.concat(customBitsets);
    return summaryTable;
  }

  createLogoSummaryTableGrid(): DG.Grid {
    // const summaryTable = this.createLogoSummaryTable();
    const isDfFiltered = this.dataFrame.filter.anyFalse;
    const filteredDf = isDfFiltered ? this.dataFrame.clone(this.dataFrame.filter) : this.dataFrame;
    const aggColsEntries = Object.entries(this.model.settings.columns ?? {});
    const aggColNames = aggColsEntries.map(([colName, aggFn]) => getAggregatedColName(aggFn, colName));

    const grid = this.logoSummaryTable.plot.grid();
    grid.sort([C.LST_COLUMN_NAMES.MEMBERS], [false]);
    this.updateFilter();
    const gridClustersCol = grid.col(C.LST_COLUMN_NAMES.CLUSTER)!;
    gridClustersCol.visible = true;
    grid.columns.setOrder([C.LST_COLUMN_NAMES.CLUSTER, C.LST_COLUMN_NAMES.MEMBERS,
      C.LST_COLUMN_NAMES.WEB_LOGO, C.LST_COLUMN_NAMES.DISTRIBUTION, C.LST_COLUMN_NAMES.MEAN_DIFFERENCE,
      C.LST_COLUMN_NAMES.P_VALUE, C.LST_COLUMN_NAMES.RATIO, ...aggColNames]);
    grid.columns.rowHeader!.visible = false;
    grid.props.rowHeight = 55;

    const webLogoCache = new DG.LruCache<number, DG.Viewer & IWebLogoViewer>();
    const distCache = new DG.LruCache<number, DG.Viewer<DG.IHistogramLookSettings>>();
    const maxSequenceLen = this.model.positionColumns.toArray().length;
    const webLogoGridCol = grid.columns.byName(C.LST_COLUMN_NAMES.WEB_LOGO)!;
    webLogoGridCol.cellType = 'html';
    webLogoGridCol.width = 350;
    const activityCol = filteredDf.getCol(C.COLUMNS_NAMES.ACTIVITY);
    const pepCol: DG.Column<string> = filteredDf.getCol(this.sequenceColumnName);

    grid.onCellRender.subscribe(async (gridCellArgs) => {
      const gridCell = gridCellArgs.cell;
      const currentRowIdx = gridCell.tableRowIndex;
      if (!gridCell.isTableCell || currentRowIdx === null || currentRowIdx === -1)
        return;

      const canvasContext = gridCellArgs.g;
      const bound = gridCellArgs.bounds;
      canvasContext.save();
      canvasContext.beginPath();
      canvasContext.rect(bound.x, bound.y, bound.width, bound.height);
      canvasContext.clip();

      try {
        const height = Math.max(gridCell.bounds.height - 2, 0);
        const clusterBitSet = this.bitsets[currentRowIdx];

        if (gridCell.tableColumn?.name === C.LST_COLUMN_NAMES.CLUSTER) {
          CR.renderLogoSummaryCell(canvasContext, gridCell.cell.value, this.clusterSelection, bound);
          gridCellArgs.preventDefault();
        } else if (gridCell.tableColumn?.name === C.LST_COLUMN_NAMES.WEB_LOGO) {
          const positionWidth = Math.floor((gridCell.bounds.width - 2 - (4 * (maxSequenceLen - 1))) / maxSequenceLen);

          let viewer = webLogoCache.get(currentRowIdx);
          if (viewer !== undefined) {
            const viewerProps = viewer.getProperties();

            for (const prop of viewerProps) {
              if (prop.name === 'positionHeight' && prop.get(viewer) !== this.webLogoMode)
                prop.set(viewer, this.webLogoMode);
              else if (prop.name === 'positionWidth' && prop.get(viewer) !== positionWidth)
                prop.set(viewer, positionWidth);
              else if (prop.name === 'minHeight' && prop.get(viewer) !== height)
                prop.set(viewer, height);
            }
            const viewerRoot = $(viewer.root).css('height', `${height}px`);//;
            viewerRoot.children().first().css('overflow-y', 'hidden !important');
          } else {
            const webLogoTable = this.createWebLogoDf(pepCol, clusterBitSet);
            viewer = await webLogoTable.plot
              .fromType('WebLogo', {positionHeight: this.webLogoMode, horizontalAlignment: HorizontalAlignments.LEFT,
                maxHeight: 1000, minHeight: height, positionWidth: positionWidth, showPositionLabels: false});
            webLogoCache.set(currentRowIdx, viewer);
          }
          gridCell.element = viewer.root;
          gridCellArgs.preventDefault();
        } else if (gridCell.tableColumn?.name === C.LST_COLUMN_NAMES.DISTRIBUTION) {
          let viewer = distCache.get(currentRowIdx);
          if (viewer === undefined) {
            const distributionDf = this.createDistributionDf(activityCol, clusterBitSet);
            viewer = distributionDf.plot.histogram({
              filteringEnabled: false,
              valueColumnName: C.COLUMNS_NAMES.ACTIVITY,
              splitColumnName: C.COLUMNS_NAMES.SPLIT_COL,
              legendVisibility: 'Never',
              showXAxis: false,
              showColumnSelector: false,
              showRangeSlider: false,
              showBinSelector: false,
              backColor: DG.Color.toHtml(DG.Color.white),
              xAxisHeight: 1,
            });
            viewer.root.style.width = 'auto';
            distCache.set(currentRowIdx, viewer);
          }
          viewer.root.style.height = `${height}px`;
          gridCell.element = viewer.root;
          gridCellArgs.preventDefault();
        }
      } finally {
        canvasContext.restore();
      }
    });
    grid.root.addEventListener('mouseleave', (_ev) => this.model.unhighlight());
    DG.debounce(grid.onCurrentCellChanged, 500).subscribe((gridCell) => {
      if (!gridCell.isTableCell)
        return;

      try {
        if (!this.keyPress || gridCell.tableColumn?.name !== C.LST_COLUMN_NAMES.CLUSTER)
          return;
        if (this.currentRowIndex !== null && this.currentRowIndex !== -1)
          this.modifyClusterSelection(this.getCluster(grid.cell(C.LST_COLUMN_NAMES.CLUSTER, this.currentRowIndex)), {shiftPressed: true, ctrlPressed: true}, false);

        this.modifyClusterSelection(this.getCluster(gridCell), {shiftPressed: true, ctrlPressed: false});
        grid.invalidate();
      } finally {
        this.keyPress = false;
        this.currentRowIndex = gridCell.gridRow;
      }
    });
    grid.root.addEventListener('keydown', (ev) => {
      this.keyPress = ev.key.startsWith('Arrow');
      if (this.keyPress)
        return;
      if (ev.key === 'Escape' || (ev.code === 'KeyA' && ev.shiftKey && ev.ctrlKey))
        this.initClusterSelection({notify: false});
      else if (ev.code === 'KeyA' && ev.ctrlKey) {
        for (let rowIdx = 0; rowIdx < this.logoSummaryTable.rowCount; ++rowIdx) {
          this.modifyClusterSelection(this.getCluster(grid.cell(C.LST_COLUMN_NAMES.CLUSTER, rowIdx)),
            {shiftPressed: true, ctrlPressed: false}, false);
        }
      }
      this.model.fireBitsetChanged();
      grid.invalidate();
    });
    grid.root.addEventListener('click', (ev) => {
      const gridCell = grid.hitTest(ev.offsetX, ev.offsetY);
      if (!gridCell || !gridCell.isTableCell || gridCell.tableColumn?.name !== C.LST_COLUMN_NAMES.CLUSTER)
        return;

      const selection = this.getCluster(gridCell);
      this.modifyClusterSelection(selection, {shiftPressed: ev.shiftKey, ctrlPressed: ev.ctrlKey});
      grid.invalidate();

      _package.files.readAsText('help/logo-summary-table.md').then((text) => {
        grok.shell.windows.help.showHelp(ui.markdown(text));
      }).catch((e) => grok.log.error(e));
    });
    grid.onCellTooltip((gridCell, x, y) => {
      if (!gridCell.isTableCell) {
        this.model.unhighlight();
        return true;
      }

      const cluster = this.getCluster(gridCell);
      this.highlightCluster(cluster);
      if (gridCell.tableColumn?.name === C.LST_COLUMN_NAMES.CLUSTER)
        this.showTooltip(cluster, x, y);
      return true;
    });

    const gridProps = grid.props;
    gridProps.allowEdit = false;
    gridProps.allowRowSelection = false;
    gridProps.allowBlockSelection = false;
    gridProps.allowColSelection = false;
    gridProps.showCurrentRowIndicator = false;

    return grid;
  }

  highlightCluster(cluster: type.SelectionItem): void {
    const bitArray = this.clusterStats[cluster.positionOrClusterType as ClusterType][cluster.monomerOrCluster].mask;
    this.dataFrame.rows.highlight((i) => bitArray.getBit(i));
    this.model.isHighlighting = true;
  }

  getCluster(gridCell: DG.GridCell): SelectionItem {
    const clustName = this.logoSummaryTable.get(C.LST_COLUMN_NAMES.CLUSTER, gridCell.tableRowIndex!);
    const clustColCat = this.dataFrame.getCol(this.clustersColumnName).categories;
    return {positionOrClusterType: clustColCat.includes(clustName) ? CLUSTER_TYPE.ORIGINAL : CLUSTER_TYPE.CUSTOM,
      monomerOrCluster: clustName};
  }

  updateFilter(): void {
    const memberstCol = this.logoSummaryTable.getCol(C.LST_COLUMN_NAMES.MEMBERS);
    const membersColData = memberstCol.getRawData();
    const maxCount = memberstCol.stats.max;
    const minMembers = Math.ceil(maxCount * this.membersRatioThreshold);
    this.logoSummaryTable.filter.init((i) => membersColData[i] > minMembers);
  }

  clusterFromSelection(): void {
    const currentSelection = this.model.getVisibleSelection();
    const viewerDfCols = this.logoSummaryTable.columns;
    const viewerDfColsLength = viewerDfCols.length;
    const newClusterVals = new Array(viewerDfCols.length);
    const activityScaledCol = this.dataFrame.getCol(C.COLUMNS_NAMES.ACTIVITY);
    const bitArray = BitArray.fromString(currentSelection.toBinaryString());
    const stats = getStats(activityScaledCol.getRawData(), bitArray);

    this.bitsets.push(currentSelection.clone());

    const newClusterName = this.dataFrame.columns.getUnusedName('New Cluster');
    const aggregatedValues: {[colName: string]: number} = {};
    const aggColsEntries = Object.entries(this.model.settings.columns ?? {});
    for (const [colName, aggFn] of aggColsEntries) {
      const newColName = getAggregatedColName(aggFn, colName);
      const col = this.dataFrame.getCol(colName);
      aggregatedValues[newColName] = getAggregatedValue(col, aggFn, currentSelection);
    }

    for (let i = 0; i < viewerDfColsLength; ++i) {
      const col = viewerDfCols.byIndex(i);
      newClusterVals[i] = col.name === C.LST_COLUMN_NAMES.CLUSTER ? newClusterName :
        col.name === C.LST_COLUMN_NAMES.MEMBERS ? currentSelection.trueCount :
          col.name === C.LST_COLUMN_NAMES.WEB_LOGO ? null :
            col.name === C.LST_COLUMN_NAMES.DISTRIBUTION ? null :
              col.name === C.LST_COLUMN_NAMES.MEAN_DIFFERENCE ? stats.meanDifference:
                col.name === C.LST_COLUMN_NAMES.P_VALUE ? stats.pValue:
                  col.name === C.LST_COLUMN_NAMES.RATIO ? stats.ratio:
                    col.name in aggregatedValues ? aggregatedValues[col.name] :
                      undefined;
      if (typeof newClusterVals[i] === 'undefined')
        _package.logger.warning(`PeptidesLSTWarn: value for column ${col.name} is undefined`);
    }
    this.logoSummaryTable.rows.addNew(newClusterVals);

    this.clusterStats[CLUSTER_TYPE.CUSTOM][newClusterName] = stats;
    this.addNewCluster(newClusterName);
  }

  removeCluster(): void {
    const lss = this.clusterSelection[CLUSTER_TYPE.CUSTOM];

    // // Names of the clusters to remove
    if (lss.length === 0)
      return grok.shell.warning('No custom clusters selected to be removed');

    const viewerDfRows = this.logoSummaryTable.rows;
    const clustCol = this.logoSummaryTable.getCol(C.LST_COLUMN_NAMES.CLUSTER);
    const clustColCat = clustCol.categories;
    const dfCols = this.dataFrame.columns;

    for (const cluster of lss) {
      lss.splice(lss.indexOf(cluster), 1);
      dfCols.remove(cluster);
      delete this.clusterStats[CLUSTER_TYPE.CUSTOM][cluster];
      const clustIdx = clustColCat.indexOf(cluster);
      viewerDfRows.removeAt(clustIdx);
      this.bitsets.splice(clustIdx, 1);
    }

    clustCol.compact();

    this.clusterSelection[CLUSTER_TYPE.CUSTOM] = lss;
    // this.model.clusterSelection = this.model.clusterSelection;
    this.render();
  }

  addNewCluster(clusterName: string): void {
    const newClusterCol = DG.Column.fromBitSet(clusterName, this.model.getVisibleSelection());
    newClusterCol.setTag(C.TAGS.CUSTOM_CLUSTER, '1');
    newClusterCol.setTag(C.TAGS.ANALYSIS_COL, `${true}`);
    this.dataFrame.columns.add(newClusterCol);
    this.model.analysisView.grid.col(newClusterCol.name)!.visible = false;
  }

  modifyClusterSelection(cluster: type.SelectionItem, options: type.SelectionOptions = {shiftPressed: false, ctrlPressed: false}, notify: boolean = true): void {
    if (notify)
      this.clusterSelection = this.model.modifySelection(this.clusterSelection, cluster, options);
    else
      this._clusterSelection = this.model.modifySelection(this.clusterSelection, cluster, options);
  }

  showTooltip(cluster: SelectionItem, x: number, y: number): HTMLDivElement | null {
    const bs = this.dataFrame.filter;
    const filteredDf = bs.anyFalse ? this.dataFrame.clone(bs) : this.dataFrame;
    const rowCount = filteredDf.rowCount;
    const bitArray = new BitArray(rowCount, false);
    const activityCol = filteredDf.getCol(C.COLUMNS_NAMES.ACTIVITY);
    const activityColData = activityCol.getRawData();

    if (cluster.positionOrClusterType === CLUSTER_TYPE.ORIGINAL) {
      const origClustCol = filteredDf.getCol(this.clustersColumnName);
      const origClustColData = origClustCol.getRawData();
      const origClustColCategories = origClustCol.categories;
      const seekValue = origClustColCategories.indexOf(cluster.monomerOrCluster);

      for (let i = 0; i < rowCount; ++i) {
        if (origClustColData[i] === seekValue)
          bitArray.setTrue(i);
      }
      bitArray.incrementVersion();
    } else {
      const clustCol: DG.Column<boolean> = filteredDf.getCol(cluster.monomerOrCluster);
      bitArray.buffer = clustCol.getRawData() as Uint32Array;
    }

    const stats = bs.anyFalse ? getStats(activityColData, bitArray) :
      this.clusterStats[cluster.positionOrClusterType as ClusterType][cluster.monomerOrCluster];

    if (!stats.count)
      return null;

    const mask = DG.BitSet.fromBytes(bitArray.buffer.buffer, rowCount);
    const distributionTable = this.createDistributionDf(activityCol, mask);
    const hist = getActivityDistribution(distributionTable, true);
    const tableMap = getStatsTableMap(stats);
    const aggregatedColMap = getAggregatedColumnValues(this.dataFrame, this.model.settings.columns!, {filterDf: true, mask: mask});
    const resultMap: {[key: string]: any} = {...tableMap, ...aggregatedColMap};
    const tooltip = getDistributionPanel(hist, resultMap);

    ui.tooltip.show(tooltip, x, y);

    return tooltip;
  }

  createWebLogoDf(pepCol: DG.Column<string>, mask: DG.BitSet): DG.DataFrame {
    const newDf = DG.DataFrame.fromColumns([pepCol]);
    newDf.filter.copyFrom(mask);
    return newDf;
  }

  createDistributionDf(activityCol: DG.Column<number>, splitMask: DG.BitSet): DG.DataFrame {
    // const table = DG.DataFrame.fromColumns([activityCol, DG.Column.fromBitSet(C.COLUMNS_NAMES.SPLIT_COL, splitMask)]);
    return getDistributionTable(activityCol, splitMask);
  }
}
