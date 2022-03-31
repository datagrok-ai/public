import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {GridCell, Point} from "datagrok-api/src/grid";
import {Paint} from "datagrok-api/src/utils";
import {Color} from "datagrok-api/src/widgets";
import {MARKER_TYPE} from "datagrok-api/src/const";

function getDataColumns(gc: DG.GridColumn): DG.Column[] {
  if (gc.settings == null)
    gc.settings = gc.grid.dataFrame.columns.byNames(Array.from(gc.grid.dataFrame.columns.numerical).map((c: any) => c.name));

  return gc.settings;
}

export class SparklineCellRenderer extends DG.GridCellRenderer {
  get name() { return 'sparkline'; }
  get cellType() { return 'sparkline'; }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: GridCell, cellStyle: DG.GridCellStyle) {

    x += 2; y += 2; w -= 4; h -= 4;

    if (w < 20 || h < 10) return;

    g.strokeStyle = 'lightgrey';
    g.lineWidth = 1;

    let row = gridCell.cell.row.idx;
    let cols = getDataColumns(gridCell.gridColumn);

    function getPos(col: number, row: number): Point {
      return new Point(x + w * (cols.length == 1 ? 0 : col / (cols.length - 1)),
        (y + h) - h * cols[col].scale(row));
    }

    g.beginPath();
    for (let i = 0; i < cols.length; i++) {
      let p = getPos(i, row);

      if (i == 0 || cols[i].isNone(row))
        g.moveTo(p.x, p.y);
      else
        g.lineTo(p.x, p.y);
    }
    g.stroke();

    for (let i = 0; i < cols.length; i++) {
      let p = getPos(i, row);
      Paint.marker(g, MARKER_TYPE.CIRCLE, p.x, p.y, Color.red, 3);
    }
  }

  renderSettings(gridColumn: DG.GridColumn): HTMLElement {
    return ui.columnsInput('Sparkline columns', gridColumn.grid.dataFrame, (names: any) => {
      gridColumn.settings = names;
      gridColumn.grid.invalidate();
    }).root;
  }
}
