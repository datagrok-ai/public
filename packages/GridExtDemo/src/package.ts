import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as GridUtils from '@datagrok-libraries/gridext/src/utils/GridUtils';
import {PinnedColumn} from '@datagrok-libraries/gridext/src/pinned/PinnedColumn';
import * as PinnedUtils from '@datagrok-libraries/gridext/src/pinned/PinnedUtils';
import {ClickableTextRenderer} from "@datagrok-libraries/gridext/src/renderer/ClickableTextRenderer";
import {RendererUIManager} from "@datagrok-libraries/gridext/src/renderer/RendererUIManager";
import {DateCellRenderer} from "@datagrok-libraries/gridext/src/renderer/DateCellRenderer";


//tags: app
//name: Demo Date Renderer
export async function demoDateRenderer() {
  const nRowCount = 100;
  const nColCount = 5;
  const dframe: DG.DataFrame = grok.data.demo.randomWalk(nRowCount, nColCount);
  const colDate = dframe.columns.addNewDateTime('New Date');
  for (let nR = 0; nR < nRowCount; ++nR) {
    colDate.set(nR, DG.DateTime.fromDate(new Date()));
  }

  const view: DG.TableView = grok.shell.addTableView(dframe);
  const grid = view.grid;
  const colGrid = grid.columns.byName("New Date");
  if(colGrid !== null) {
    GridUtils.setGridColumnRenderer(colGrid, new DateCellRenderer());
  }
  RendererUIManager.register(grid);
}


//tags: app
//name: Demo Clickable Renderer
export async function demoClickableRenderer() {
  const nRowCount = 100;
  const nColCount = 500;
  const dframe: DG.DataFrame = grok.data.demo.randomWalk(nRowCount, nColCount);

  const view: DG.TableView = grok.shell.addTableView(dframe);
  const grid = view.grid;

  let colGrid : DG.GridColumn | null;
  const renderer = new ClickableTextRenderer();
  for(let nC=3; nC<4; ++nC) {
    colGrid = grid.columns.byIndex(nC);
    if(colGrid === null)
      continue;

    GridUtils.setGridColumnRenderer(colGrid, renderer);
  }

  let bSuccess = RendererUIManager.register(view.grid);
  if(!bSuccess)
    throw new Error("Failed to register Renderer UI manager");
 }

//tags: app
//name: Demo Pinned Columns
export async function demoPinnedColumns() {

  const nRowCount = 100;
  const nColCount = 500;
  const dframe : DG.DataFrame = grok.data.demo.randomWalk(nRowCount, nColCount);
  const view : DG.TableView = grok.shell.addTableView(dframe);
  view.grid.setOptions({
    colHeaderHeight: 30,
    rowHeight: 25
  });

  grok.events.onContextMenu.subscribe((args) => {
    PinnedUtils.handleContextMenu(args, (menu : DG.Menu, colGridOrPinned : DG.GridColumn | PinnedColumn, grid : DG.Grid) => {

      if(colGridOrPinned instanceof PinnedColumn) {
        const colGrid = colGridOrPinned.getGridColumn();
        if(colGrid !== null && !GridUtils.isRowHeader(colGrid)) {
        menu.item("Unpin Column", () => {
          colGridOrPinned.close();
        });
      }
          menu.item("Unpin All Columns", () => {
          PinnedUtils.closeAllPinnedColumns(grid);
        });
      }
      else {
        menu.item("Pin Column", async() => {
          PinnedUtils.addPinnedColumn(colGridOrPinned);
        });
      }
    });
  });
}



