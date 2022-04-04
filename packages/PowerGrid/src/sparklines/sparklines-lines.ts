import * as DG from 'datagrok-api/dg';
import {InputBase, Property, TYPE} from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {GridCell, Point} from "datagrok-api/src/grid";
import {Paint} from "datagrok-api/src/utils";
import {Color} from "datagrok-api/src/widgets";
import {MARKER_TYPE} from "datagrok-api/src/const";

interface SparklineSettings {
  normalize: boolean;
  columnNames: string[];
}

function names(columns: Iterable<DG.Column>): string[] {
  return Array.from(columns).map((c: any) => c.name)
}

function getSettings(gc: DG.GridColumn): SparklineSettings {
  return gc.settings ??= {
    normalize: true,
    columnNames: names(gc.grid.dataFrame.columns.numerical)
  }
}

function getDataColumns(gc: DG.GridColumn): DG.Column[] {
  return gc.grid.dataFrame.columns.byNames(getSettings(gc).columnNames);
}

export class SparklineCellRenderer extends DG.GridCellRenderer {
  get name() { return 'sparkline'; }
  get cellType() { return 'sparkline'; }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: GridCell, cellStyle: DG.GridCellStyle) {

    const settings = getSettings(gridCell.gridColumn);
    x += 4; y += 2; w -= 6; h -= 4;

    if (w < 20 || h < 10) return;

    g.strokeStyle = 'lightgrey';
    g.lineWidth = 1;

    let row = gridCell.cell.row.idx;
    let cols = gridCell.grid.dataFrame.columns.byNames(settings.columnNames);
    let gmin = settings.normalize ? 0 : Math.min(...cols.map((c: DG.Column) => c.min));
    let gmax = settings.normalize ? 0 : Math.max(...cols.map((c: DG.Column) => c.max));

    function getPos(col: number, row: number): Point {
      const r: number = settings.normalize ? cols[col].scale(row) : (cols[col].get(row) - gmin) / (gmax - gmin);
      return new Point(x + w * (cols.length == 1 ? 0 : col / (cols.length - 1)), (y + h) - h * r);
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
    gridColumn.settings ??= { normalize: true };
    const settings: SparklineSettings = gridColumn.settings;

    const normalizeInput = InputBase.forProperty(Property.js('normalize', TYPE.BOOL), settings);
    normalizeInput.onChanged(() => gridColumn.grid.invalidate());

    return ui.inputs([
      normalizeInput,
      ui.columnsInput('Sparkline columns', gridColumn.grid.dataFrame, (columns) => {
        settings.columnNames = names(columns);
        gridColumn.grid.invalidate();
        console.log(JSON.stringify(gridColumn.settings));
      }, {
        available: names(gridColumn.grid.dataFrame.columns.numerical),
        checked: settings?.columnNames ?? names(gridColumn.grid.dataFrame.columns.numerical)
      })
    ]);
  }
}
