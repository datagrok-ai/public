import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export class Plate {
  constructor() {
    this.barcode = '123456';
    this.rows = 8;
    this.cols = 12;
  }
}

export class PlateGridCellRenderer extends DG.GridCellRenderer {
  get name() { return 'Plate cell renderer'; }
  get cellType() { return 'demo_plate'; }
  get defaultWidth() { return 50; }
  get defaultHeight() { return 50; }

  render(g, x, y, w, h, gridCell, cellStyle) {
    g.fillStyle = 'black';
    g.textBaseline = 'top';
    g.fillText(`b ${gridCell.cell.value}`, x + 3, y + 3);
    g.beginPath();
    g.lineWidth = "1";
    g.strokeStyle = "gray";
    g.rect(x + 1, y + 1, w - 2, h - 2);
    g.stroke();
  }
}

// Defines the way Datagrok handles entities of the specified type
export class PlateHandler extends DG.ObjectHandler {
  get type() { return 'demo_plate' }

  // Checks whether this is the handler for [x]
  isApplicable(x) { return x instanceof Plate; }

  getCanvasRenderer(x) { return new FruitCanvasRenderer(); }
  getGridCellRenderer(x) { return new PlateGridCellRenderer(); }

  renderIcon(x) { return ui.iconFA('grip-horizontal'); }
  renderMarkup(x) { return ui.span([this.renderIcon(), ' plate ', x.barcode]); }
}