import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {GridCell} from "datagrok-api/dg";

export class TestCellRenderer extends DG.GridCellRenderer {
  get name() { return 'test'; }
  get cellType() { return 'test'; }

  get defaultWidth(): number | null { return 200; }
  get defaultHeight(): number | null { return 100; }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    g.fillStyle = 'black';
    g.fillText(`test`, x + 25, y + 10);
  }

  onKeyDown(gridCell: GridCell, e: KeyboardEvent): void { console.log('keyDown', gridCell, e); }
  onKeyPress(gridCell: GridCell, e: KeyboardEvent): void { console.log('keyPress', gridCell, e); }

  onMouseLeave(gridCell: GridCell, e: MouseEvent): void { console.log('mouseLeave', gridCell, e); }
  onMouseDown(gridCell: GridCell, e: MouseEvent): void { console.log('mouseDown', gridCell, e); }
  onMouseUp(gridCell: GridCell, e: MouseEvent): void { console.log('mouseUp', gridCell, e); }
  onClick(gridCell: GridCell, e: MouseEvent): void { console.log('click', gridCell, e); }
  onDoubleClick(gridCell: GridCell, e: MouseEvent): void { console.log('doubleClick', gridCell, e); }

  // rendering in mouse event coordinates
  onMouseMove(gridCell: GridCell, e: MouseEvent): void {
    console.log('mouseMove', gridCell, e);
    const g = gridCell.grid.canvas.getContext('2d')!;
    g.fillStyle = 'blue';
    g.fillRect(e.offsetX, e.offsetY, 2, 2);
  }

  // rendering in cell coordinates
  onMouseEnter(gridCell: GridCell, e: MouseEvent): void {
    console.log('mouseEnter', gridCell, e);
    const r = gridCell.bounds;
    const g = gridCell.grid.canvas.getContext('2d')!;
    g.fillStyle = 'red';
    g.fillRect(r.x, r.y, r.width, r.height);
  }
}


export class HtmlTestCellRenderer extends DG.GridCellRenderer {
  get name() { return 'htest'; }
  get cellType() { return 'htest'; }

  get defaultWidth(): number | null { return 200; }
  get defaultHeight(): number | null { return 100; }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    gridCell.element = ui.button('foo', () => grok.shell.info('' + gridCell.gridRow));
  }
}