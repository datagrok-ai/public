import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

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

    ui.tools.handleResize(this.canvas, (w: number, h: number) => {
      this.canvas.width = w;
      this.canvas.height = h;
      this.render();
    });

    this.canvas.onmousemove = (e: MouseEvent) => this.renderer.onMouseMove(this.gridCell!, e);
    this.canvas.onmouseenter = (e: MouseEvent) => this.renderer.onMouseEnter(this.gridCell!, e);
    this.canvas.onmouseleave = (e: MouseEvent) => this.renderer.onMouseLeave(this.gridCell!, e);
  }

  static fromGridCell(gridCell: DG.GridCell, renderer: DG.GridCellRenderer): CellRenderViewer {
    const viewer = new CellRenderViewer(renderer);
    viewer.gridCell = gridCell;
    viewer.render();
    return viewer;
  }

  render(): void {
    if (!this.gridCell)
      return;

    const g = this.canvas.getContext('2d')!;
    g.clearRect(0, 0, this.canvas.width, this.canvas.height);

    this.renderer.renderInternal(g, 0, 0, this.canvas.width, this.canvas.height,
      this.gridCell.dart, this.gridCellStyle);
  }
}
