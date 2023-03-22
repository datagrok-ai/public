import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {GridCell} from "datagrok-api/dg";

export class CellRenderViewer<TRenderer extends DG.GridCellRenderer = DG.GridCellRenderer> extends DG.JsViewer {
  canvas: HTMLCanvasElement = ui.canvas();
  renderer: TRenderer;
  gridCell?: DG.GridCell;
  gridCellStyle: DG.GridCellStyle = DG.GridCellStyle.create();

  constructor(renderer: TRenderer) {
    super();

    this.renderer = renderer;
    this.root.append(this.canvas);

    this.canvas.style.width = '100%';
    this.canvas.style.height = '100%';
    ui.tools.handleResize(this.canvas, (w, h) => {
      this.canvas.width = w;
      this.canvas.height = h;
      this.render();
    });

    const bind = (f: (gridCell: GridCell, e: MouseEvent) => void) => (e: MouseEvent) => {
      if (this.gridCell)
        f(this.gridCell, e);
    }

    this.canvas.onmousemove = bind(this.renderer.onMouseMove);
    this.canvas.onmouseenter = bind(this.renderer.onMouseEnter);
    this.canvas.onmouseleave = bind(this.renderer.onMouseLeave);
  }

  render() {
    if (!this.gridCell)
      return;

    const g = this.canvas.getContext('2d')!
    g.clearRect(0, 0, this.canvas.width, this.canvas.height);

    this.renderer.renderInternal(g, 0, 0, this.canvas.width, this.canvas.height,
      this.gridCell.dart, this.gridCellStyle);
  }
}
