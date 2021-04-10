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

export class BriefPlateRenderer extends DG.GridCellRenderer {
  get name() { return 'Brief plate renderer'; }
  get cellType() { return 'demo_plate_brief'; }

  render(g, x, y, w, h, gridCell, cellStyle) {
    g.fillStyle = 'black';
    g.textBaseline = 'top';
    g.fillText(`b ${gridCell.cell.value}`, x + 3, y + 3);
  }
}

export class FullPlateRenderer extends DG.GridCellRenderer {
  get name() { return 'Full plate renderer'; }
  get cellType() { return 'demo_plate_full'; }
  get defaultWidth() { return 50; }
  get defaultHeight() { return 50; }

  render(g, x, y, w, h, gridCell, cellStyle) {
    g.fillStyle = 'black';
    g.textBaseline = 'top';
    g.fillText(`b ${gridCell.cell.value}`, x + 3, y + 3);
    g.beginPath();
    g.lineWidth = 1;
    g.strokeStyle = "gray";
    g.rect(x + 1, y + 1, w - 2, h - 2);
    g.stroke();
  }
}


// Defines the way Datagrok handles entities of the specified type
export class FullPlateHandler extends DG.ObjectHandler {
  get type() { return 'demo_plate'; }
  get name() { return 'full handler' }

  // Checks whether this is the handler for [x]
  isApplicable(x) { return x.barcode; }
  getGridCellRenderer(x) { return new FullPlateRenderer(); }
  renderIcon(x) { return ui.iconFA('grip-horizontal'); }
  renderMarkup(x) { return ui.span([this.renderIcon(), ' plate ', x.barcode]); }
}


// Defines the way Datagrok handles entities of the specified type
export class BriefPlateHandler extends DG.ObjectHandler {
  get type() { return 'brief handler'; }

  // Checks whether this is the handler for [x]
  isApplicable(x) { return x.barcode; }
  getGridCellRenderer(x) { return new BriefPlateRenderer(); }
  renderIcon(x) { return ui.iconFA('grip-horizontal'); }
  renderMarkup(x) { return ui.span([this.renderIcon(), ' plate ', x.barcode]); }
}