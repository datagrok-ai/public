import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


export class TestCellRenderer extends DG.GridCellRenderer {
  get name() { return 'test'; }

  get cellType() { return 'testUnitsKg'; }

  get defaultWidth(): number | null { return 200; }

  get defaultHeight(): number | null { return 100; }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    g.fillStyle = 'black';
    g.fillText(`test`, x + 25, y + 10);
  }

  onKeyDown(gridCell: DG.GridCell, e: KeyboardEvent): void { console.log('keyDown', gridCell, e); }

  onKeyPress(gridCell: DG.GridCell, e: KeyboardEvent): void { console.log('keyPress', gridCell, e); }

  onMouseLeave(gridCell: DG.GridCell, e: MouseEvent): void { console.log('mouseLeave', gridCell, e); }

  onMouseDown(gridCell: DG.GridCell, e: MouseEvent): void { console.log('mouseDown', gridCell, e); }

  onMouseUp(gridCell: DG.GridCell, e: MouseEvent): void { console.log('mouseUp', gridCell, e); }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void { console.log('click', gridCell, e); }

  onDoubleClick(gridCell: DG.GridCell, e: MouseEvent): void { console.log('doubleClick', gridCell, e); }

  // rendering in mouse event coordinates
  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    console.log('mouseMove', gridCell, e);
    const g = gridCell.grid.canvas.getContext('2d')!;
    g.fillStyle = 'blue';
    g.fillRect(e.offsetX, e.offsetY, 2, 2);
  }

  // rendering in cell coordinates
  onMouseEnter(gridCell: DG.GridCell, e: MouseEvent): void {
    console.log('mouseEnter', gridCell, e);
    const r = gridCell.bounds;
    const g = gridCell.grid.canvas.getContext('2d')!;
    g.fillStyle = 'red';
    g.fillRect(r.x, r.y, r.width, r.height);
  }
}


@grok.decorators.cellRenderer({
  name: 'htestCellRenderer',
  cellType: 'htest',
})
export class HtmlTestCellRenderer extends DG.GridCellRenderer {
  get name() { return 'htest'; }

  get cellType() { return 'testUnitsTon'; }

  get defaultWidth(): number | null { return 200; }

  get defaultHeight(): number | null { return 100; }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    gridCell.element = ui.button('foo', () => grok.shell.info('' + gridCell.gridRow));
  }
}
