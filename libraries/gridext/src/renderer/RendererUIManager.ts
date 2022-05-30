import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as rxjs from 'rxjs';
import * as GridUtils from '../utils/GridUtils';
import {GridCellRendererEx} from "./GridCellRendererEx";

function getRenderer(colGrid : DG.GridColumn) : GridCellRendererEx | null {
  //const colGrid = cellGrid.gridColumn;
  if (colGrid === null || colGrid === undefined) {
    throw new Error('Grid cell is detached from the Grid column');
  }

  let renderer : GridCellRendererEx | null = GridUtils.getGridColumnRenderer(colGrid);
  if(renderer !== null) {
    return renderer;
  }

  //renderer = cellGrid.renderer;
  //return renderer;
  return null;
}

function convertToCellXY(arXY : Array<number>, cellGrid : DG.GridCell, e : MouseEvent) : void {
  const grid : DG.Grid = cellGrid.grid;
  const colG = cellGrid.gridColumn;
  const sbHorz = grid.horzScroll;
  const nXOnCell = e.offsetX - colG.left + Math.floor(sbHorz.min);

  const nHColHeader = GridUtils.getGridColumnHeaderHeight(grid);
  const nHRow = GridUtils.getGridRowHeight(grid);
  const nYOnCell = e.offsetY - nHColHeader - nHRow * (cellGrid.gridRow - Math.floor(grid.vertScroll.min));
  arXY[0] = nXOnCell;
  arXY[1] = nYOnCell;
}



class RendererUIManagerImpl {
  constructor(grid : DG.Grid) {
    this.m_grid = grid;
    const dart = DG.toDart(grid);
    dart.m_managerUIRenderer = this;

    this.m_handlerCellRender = grid.onCellRender.subscribe((args) => {

      if (args.cell.isTableCell) {
        const renderer = args.cell.gridColumn === null ? null : GridUtils.getGridColumnRenderer(args.cell.gridColumn);
        if(renderer !== null) {
          renderer.render(args.g, args.bounds.x, args.bounds.y, args.bounds.width, args.bounds.height, args.cell, args.cell.style);
          args.preventDefault();
        }
      }
    });

    this.m_handlerColResized = grid.onColumnResized.subscribe((ev) => {

      const colGrid = ev.args.column;
      const renderer = getRenderer(colGrid);

      if (renderer !== null) {
        const grid = colGrid.grid;
        const bAdjusting = ev.args.dragging;
        const nW = colGrid.width;
        renderer.onResizeWidth(colGrid, grid, nW, bAdjusting);
      }
    });

    const arColRowIdxs = new Array<number>(4);
    this.m_handlerRowResized = grid.onRowsResized.subscribe((ev) => {

      const bAdjusting = ev.args.dragging;
      GridUtils.fillVisibleViewportGridCells(arColRowIdxs, grid);

      const nHRow = GridUtils.getGridRowHeight(grid);

      const nColMin = arColRowIdxs[0];
      const nColMax = arColRowIdxs[1];

      let colGrid = null;
      let renderer= null;

      //console.log("Requesting cols " + nColMin + " " + nColMax + " for " + nHRow);
      let nColVisited = 0;
      for(let nCol=nColMin; nCol<=nColMax; ++nCol) {
        if(nCol === 0)
          continue;

        colGrid = grid.columns.byIndex(nCol);
        if(colGrid === null || !colGrid.visible)
          continue;

        renderer = getRenderer(colGrid);
        if(renderer === null)
          continue;

        renderer.onResizeHeight(colGrid, grid, nHRow, bAdjusting);
        ++nColVisited;
      }
      //grok.shell.info("Resizing row height: " + (ev.args.dragging ? "in progress" : "done"));
    });



    const arXY = [-1,-1];
    let nXDown = -1;
    let nYDown = -1;
    this.m_handlerMouseDown = rxjs.fromEvent(grid.overlay, 'mousedown').subscribe((e) => {
        const eMouse = e as MouseEvent;
        const cell = grid.hitTest(eMouse.offsetX, eMouse.offsetY);
        if (cell === null || cell === undefined || cell.dart === undefined || !cell.isTableCell || cell.gridColumn === null) {
          return;
        }
        const renderer = getRenderer(cell.gridColumn);
        if(renderer instanceof GridCellRendererEx) {
          convertToCellXY(arXY, cell, eMouse);
          renderer.onMouseDownEx(cell, eMouse, arXY[0], arXY[1]);
        }
        nXDown = eMouse.offsetX;
        nYDown = eMouse.offsetY;
      });

    this.m_handlerMouseUp = rxjs.fromEvent(grid.overlay, 'mouseup').subscribe((e) => {
        const eMouse = e as MouseEvent;
        const cell = grid.hitTest(eMouse.offsetX, eMouse.offsetY);
        if (cell === null || cell === undefined || cell.dart === undefined || !cell.isTableCell || cell.gridColumn === null) {
          return;
        }

      const renderer = getRenderer(cell.gridColumn);
      if(renderer instanceof GridCellRendererEx) {
        convertToCellXY(arXY, cell, eMouse);
        renderer.onMouseUpEx(cell, eMouse, arXY[0], arXY[1]);
      }
      });

    this.m_handlerClick = rxjs.fromEvent(grid.overlay, 'click').subscribe((e) => {
        const eMouse = e as MouseEvent;
        const cell = grid.hitTest(eMouse.offsetX, eMouse.offsetY);
        if (cell === null || cell === undefined || cell.dart === undefined || !cell.isTableCell || cell.gridColumn === null) {
          return;
        }
        if (nXDown === eMouse.offsetX && nYDown === eMouse.offsetY) {
          const renderer = getRenderer(cell.gridColumn);
          if(renderer instanceof GridCellRendererEx) {
            convertToCellXY(arXY, cell, eMouse);
            renderer.onClickEx(cell, eMouse, arXY[0], arXY[1]);
          }
        }

        nXDown = -1;
        nYDown = -1;
      });
    let cellCurrent : DG.GridCell | null = null;
    this.m_handlerMouseMove = rxjs.fromEvent(grid.overlay, 'mousemove').subscribe((e) => {
        const eMouse = e as MouseEvent;

        let cell = grid.hitTest(eMouse.offsetX, eMouse.offsetY);
        if (cell === null || cell === undefined || cell.dart === undefined || !cell.isTableCell || cell.gridColumn === null) {
          cell = null;
        }
        if (cellCurrent === null && cell != null) {
          const renderer = getRenderer(cell.gridColumn);
          if(renderer instanceof GridCellRendererEx) {
            convertToCellXY(arXY, cell, eMouse);
            renderer.onMouseEnterEx(cell, eMouse, arXY[0], arXY[1]);
          }
        }

        if (cellCurrent !== null && cellCurrent.gridColumn !== null && (cell === null || cellCurrent.gridColumn.name !== cell.gridColumn.name || cellCurrent.gridRow !== cell.gridRow)) {

          const renderer = getRenderer(cellCurrent.gridColumn);
          if(renderer instanceof GridCellRendererEx) {
            convertToCellXY(arXY, cellCurrent, eMouse);
            renderer.onMouseLeaveEx(cellCurrent, eMouse, arXY[0], arXY[1]);
          }
          //onMouseEvent(e, cellCurrent, 'onMouseLeave');
          if (cell !== null && cell.gridColumn !== null) {
            const renderer = getRenderer(cell.gridColumn);
            if(renderer instanceof GridCellRendererEx) {
              convertToCellXY(arXY, cell, eMouse);
              renderer.onMouseEnterEx(cell, eMouse, arXY[0], arXY[1]);
            }
            //onMouseEvent(e, cell, 'onMouseEnter');
          }
        }

        if (cell !== null && cell.gridColumn !== null) {
          const renderer = getRenderer(cell.gridColumn);
          if(renderer instanceof GridCellRendererEx) {
            convertToCellXY(arXY, cell, eMouse);
            renderer.onMouseMoveEx(cell, eMouse, arXY[0], arXY[1]);
          }
        }

        cellCurrent = cell;
      });

    this.m_handlerMouseOut = rxjs.fromEvent(grid.overlay, 'mouseout').subscribe((e) => {
        if (cellCurrent !== null && cellCurrent.gridColumn !== null) {

          const renderer = getRenderer(cellCurrent.gridColumn);
          if(renderer instanceof GridCellRendererEx) {
            const eMouse = e as MouseEvent;
            convertToCellXY(arXY, cellCurrent, eMouse);
            renderer.onMouseLeaveEx(cellCurrent, eMouse, arXY[0], arXY[1]);
          }
        }

        cellCurrent = null;
      });

    const managerThis = this;
    this.m_handlerViewerClosed = grok.events.onViewerClosed.subscribe((args) => {
      const viewer = args.args.viewer;
      if(DG.toDart(viewer) === DG.toDart(this.m_grid)){
        if(this.m_grid === null)
          throw new Error("Grid cannot be null.");

        managerThis.dispose();
      }
    });
  }

  dispose() : void {

    this.m_handlerCellRender.unsubscribe();
    this.m_handlerCellRender = null;

    this.m_handlerColResized.unsubscribe();
    this.m_handlerColResized = null;

    this.m_handlerRowResized.unsubscribe();
    this.m_handlerRowResized = null;

    this.m_handlerMouseDown.unsubscribe();
    this.m_handlerMouseDown = null;

    this.m_handlerMouseUp.unsubscribe();
    this.m_handlerMouseUp = null;

    this.m_handlerClick.unsubscribe();
    this.m_handlerClick = null;

    this.m_handlerMouseMove.unsubscribe();
    this.m_handlerMouseMove = null;

    this.m_handlerMouseOut.unsubscribe();
    this.m_handlerMouseOut = null;

    this.m_handlerViewerClosed.unsubscribe();
    this.m_handlerViewerClosed = null;

    const dart = DG.toDart(this.m_grid);
    dart.m_managerUIRenderer = null;
    this.m_grid = null;
  }

  private m_grid : DG.Grid | null;
  private m_handlerCellRender : any;
  private m_handlerColResized: any;
  private m_handlerRowResized: any;
  private m_handlerMouseDown : any;
  private m_handlerMouseUp : any;
  private m_handlerClick : any;
  private m_handlerMouseMove : any;
  private m_handlerMouseOut : any;
  private m_handlerViewerClosed: any;
}

export class RendererUIManager {
  constructor() {
    throw new Error("Cannot create instances of the RendererUIManager class");
  }

  static isRegistered(grid : DG.Grid) : boolean {
    const dart = DG.toDart(grid);
    const b = dart.m_managerUIRenderer instanceof RendererUIManagerImpl;
    return b;
  }


  static register(grid : DG.Grid) : boolean {
    if(RendererUIManager.isRegistered(grid)) {
      return false;
    }

    const manager = new RendererUIManagerImpl(grid);
    return true;
  }

  static unregister(grid : DG.Grid) : boolean {
    if(!RendererUIManager.isRegistered(grid)) {
      return false;
    }

    const dart = DG.toDart(grid);
    const manager = dart.m_managerUIRenderer;
    manager.dispose();
    return true;
  }
}
