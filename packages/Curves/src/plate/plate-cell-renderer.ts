import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PlateWidget} from './plate-widget';
import {GridColumn, x} from 'datagrok-api/dg';

export class PlateCellHandler extends DG.ObjectHandler {
  get type(): string { return 'Plate'; }

  isApplicable(x: any): boolean {
    return x instanceof DG.GridCell && x.cellType === 'Plate';
  }

  getGridCellRenderer(): DG.GridCellRenderer | null {
    return new PlateGridCellRenderer();
  }
}


@grok.decorators.cellRenderer({cellType: 'Plate', name: 'PlateGridCellRenderer'})
export class PlateGridCellRenderer extends DG.GridCellRenderer {
  plate = new PlateWidget();

  get name(): string { return 'Plate'; }
  get cellType(): string { return 'Plate'; }

  getDefaultSize(_gridColumn: GridColumn) { return {width: 240, height: 160}; }

  // eslint-disable-next-line max-len
  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, _cellStyle: DG.GridCellStyle) {
    if (gridCell.value == null) {
      g.fillStyle = 'white';
      g.fillRect(x, y, w, h);
      return;
    }

    this.plate.plateData = gridCell.value;
    this.plate.grid.props.colHeaderHeight = h > 120 ? 16 : 0;
    this.plate.grid.props.showRowHeader = w > 100;

    this.plate.grid.render(g, new DG.Rect(x, y, w, h).cutTop(Math.min(h / 20, 5)).cutRight(Math.min(w / 20, 5)));
  }

  onClick(gridCell: DG.GridCell, _e: MouseEvent) {
    grok.shell.o = PlateWidget.detailedView(gridCell.value).root;
  }
}
