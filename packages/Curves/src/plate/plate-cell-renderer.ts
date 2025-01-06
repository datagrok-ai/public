import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PlateWidget} from "./plate-widget";
import {GridColumn, x} from "datagrok-api/dg";

export class PlateCellHandler extends DG.ObjectHandler {
  get type(): string { return 'Plate'; }

  isApplicable(x: any): boolean {
    return x instanceof DG.GridCell && x.cellType === 'Plate';
  }

  getGridCellRenderer(): DG.GridCellRenderer | null {
    return new PlateGridCellRenderer();
  }
}


@grok.decorators.cellRenderer({cellType: 'Plate'})
export class PlateGridCellRenderer extends DG.GridCellRenderer {
  plate = new PlateWidget();

  get name(): string { return 'Plate'; }
  get cellType(): string { return 'Plate'; }

  getDefaultSize(gridColumn: GridColumn) { return { width: 120, height: 80 }  }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    this.plate.plateData = gridCell.value;
    this.plate.grid.props.colHeaderHeight = h > 120 ? 16 : 0;
    this.plate.grid.props.showRowHeader = w > 100;

    this.plate.grid.render(g, new DG.Rect(x, y, w, h).cutTop(Math.min(h / 20, 5)).cutRight(Math.min(w / 20, 5)));
  }

  onClick(gridCell: DG.GridCell, e: MouseEvent) {
    grok.shell.o = PlateWidget.detailedView(gridCell.value).root;
  }
}
