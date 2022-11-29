import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {PeptidesModel} from '../model';
import * as C from '../utils/constants';
import * as CR from '../utils/cell-renderer';
import * as bio from '@datagrok-libraries/bio';

export class LogoSummary extends DG.JsViewer {
  _titleHost = ui.divText('Logo Summary Table', {id: 'pep-viewer-title'});
  model!: PeptidesModel;
  viewerGrid!: DG.Grid;
  initialized: boolean = false;
  webLogoMode: string;
  importanceThreshold: number;

  constructor() {
    super();

    this.webLogoMode = this.string('webLogoMode', bio.PositionHeight.full,
      {choices: [bio.PositionHeight.full, bio.PositionHeight.Entropy]});
    this.importanceThreshold = this.float('importanceThreshold', 0.7);
  }

  onTableAttached(): void {
    super.onTableAttached();

    this.model = PeptidesModel.getInstance(this.dataFrame);
    this.subs.push(this.model.onSettingsChanged.subscribe(() => {
      this.createLogoSummaryGrid();
      this.render();
    }));
    // this.subs.push(this.model.onLogoSummaryGridChanged.subscribe((grid) => {
    //   this.viewerGrid = grid;
    //   this.render();
    // }));
    // this.model.updateDefault();
    // this.viewerGrid = this.model.logoSummaryGrid;
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

  onPropertyChanged(): void {
    this.createLogoSummaryGrid();
    this.render();
  }

  createLogoSummaryGrid(): DG.Grid {
    let summaryTableBuilder = this.dataFrame.groupBy([C.COLUMNS_NAMES.CLUSTERS]);
    for (const [colName, aggregationFunc] of Object.entries(this.model.settings.columns ?? {}))
      summaryTableBuilder = summaryTableBuilder.add(aggregationFunc as any, colName, `${aggregationFunc}(${colName})`);

    const summaryTable = summaryTableBuilder.aggregate();
    const summaryTableLength = summaryTable.rowCount;
    const clustersCol: DG.Column<number> = summaryTable.getCol(C.COLUMNS_NAMES.CLUSTERS);
    const membersCol: DG.Column<number> = summaryTable.columns.addNewInt('Members');
    const webLogoCol: DG.Column<string> = summaryTable.columns.addNew('WebLogo', DG.COLUMN_TYPE.STRING);
    const tempDfList: DG.DataFrame[] = new Array(summaryTableLength);
    const originalClustersCol = this.dataFrame.getCol(C.COLUMNS_NAMES.CLUSTERS);
    const peptideCol: DG.Column<string> = this.dataFrame.getCol(C.COLUMNS_NAMES.MACROMOLECULE);

    for (let index = 0; index < summaryTableLength; ++index) {
      const indexes: number[] = [];
      for (let j = 0; j < originalClustersCol.length; ++j) {
        if (originalClustersCol.get(j) === clustersCol.get(index))
          indexes.push(j);
      }
      const tCol = DG.Column.string('peptides', indexes.length);
      tCol.init((i) => peptideCol.get(indexes[i]));

      for (const tag of peptideCol.tags)
        tCol.setTag(tag[0], tag[1]);

      const uh = new bio.UnitsHandler(tCol);
      tCol.setTag(bio.TAGS.alphabetSize, uh.getAlphabetSize().toString());

      const dfSlice = DG.DataFrame.fromColumns([tCol]);
      tempDfList[index] = dfSlice;
      webLogoCol.set(index, index.toString());
      membersCol.set(index, dfSlice.rowCount);
      //TODO: user should be able to choose threshold
      const maxCount = this.model.clusterStatsDf.getCol(C.COLUMNS_NAMES.COUNT).stats.max;
      if (dfSlice.rowCount <= Math.ceil(maxCount * this.importanceThreshold))
        summaryTable.filter.set(index, false, false);
    }
    webLogoCol.setTag(DG.TAGS.CELL_RENDERER, 'html');

    this.viewerGrid = summaryTable.plot.grid();
    const gridClustersCol = this.viewerGrid.col(C.COLUMNS_NAMES.CLUSTERS)!;
    gridClustersCol.name = 'Clusters';
    gridClustersCol.visible = true;
    this.viewerGrid.columns.rowHeader!.visible = false;
    this.viewerGrid.props.rowHeight = 55;
    this.viewerGrid.onCellPrepare((cell) => {
      if (cell.isTableCell && cell.tableColumn?.name === 'WebLogo') {
        tempDfList[parseInt(cell.cell.value)].plot
          .fromType('WebLogo', {maxHeight: 50, positionHeight: this.webLogoMode})
          .then((viewer) => cell.element = viewer.root);
      }
    });
    this.viewerGrid.root.addEventListener('click', (ev) => {
      const cell = this.viewerGrid.hitTest(ev.offsetX, ev.offsetY);
      if (!cell || !cell.isTableCell)
        return;

      const cluster = clustersCol.get(cell.tableRowIndex!)!;
      summaryTable.currentRowIdx = -1;
      if (ev.shiftKey)
        this.model.modifyClusterSelection(cluster);
      else
        this.model.initClusterSelection(cluster);
    });
    this.viewerGrid.onCellRender.subscribe((gridCellArgs) => {
      const gc = gridCellArgs.cell;
      if (gc.tableColumn?.name !== C.COLUMNS_NAMES.CLUSTERS || gc.isColHeader)
        return;
      const canvasContext = gridCellArgs.g;
      const bound = gridCellArgs.bounds;
      canvasContext.save();
      canvasContext.beginPath();
      canvasContext.rect(bound.x, bound.y, bound.width, bound.height);
      canvasContext.clip();
      CR.renderLogoSummaryCell(canvasContext, gc.cell.value, this.model.logoSummarySelection, bound);
      gridCellArgs.preventDefault();
      canvasContext.restore();
    });
    this.viewerGrid.onCellTooltip((cell, x, y) => {
      if (!cell.isColHeader && cell.tableColumn?.name === C.COLUMNS_NAMES.CLUSTERS)
        this.model.showTooltipCluster(cell.cell.value, x, y);
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
}
