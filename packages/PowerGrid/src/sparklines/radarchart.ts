import * as DG from 'datagrok-api/dg';
import {InputBase, Property, TYPE} from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {GridCell, GridColumn, Point, Rect} from "datagrok-api/src/grid";
import {Paint} from "datagrok-api/src/utils";
import {Color} from "datagrok-api/src/widgets";
import {MARKER_TYPE} from "datagrok-api/src/const";
import {getSettingsBase, names, SummarySettingsBase} from "./base";

import "../rect-extensions";
import "../canvas-extensions";

class it {
  static range = (n: number) => [...Array(n).keys()]
}


interface RadarchartSettings extends SummarySettingsBase {
  // radius: number;
}

function getSettings(gc: DG.GridColumn): RadarchartSettings {
  return gc.settings ??= {
    ...getSettingsBase(gc),
    // ...{radius: 10,},
  }
}

export class RadarchartCellRender extends DG.GridCellRenderer {
  get name() { return 'radar ts'; }

  get cellType() { return 'radarchart_ts'; }

  // getPreferredCellSize(col: DG.GridColumn) {
  //   return new Size(80,80);
  // }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: GridCell, cellStyle: DG.GridCellStyle
  ) {
    const settings = getSettings(gridCell.gridColumn);

    if (w < 20 || h < 10) return;

    const b = new Rect(x, y, w, h).inflate(-2, -2);

    const row = gridCell.cell.row.idx;
    const cols = gridCell.grid.dataFrame.columns.byNames(settings.columnNames);
    const box = new Rect(x, y, w, h).fitSquare().inflate(-2, -2);

    g.strokeStyle = 'lightgray';

    // axises' point calculator
    const p = (col: number, ratio: number) => new Point(
      box.midX + ratio * box.width * Math.cos(2 * Math.PI * col / (cols.length)) / 2,
      box.midY + ratio * box.width * Math.sin(2 * Math.PI * col / (cols.length)) / 2);

    // Axises
    for (let i = 0; i < cols.length; i++) {
      g.line(box.midX, box.midY, p(i, 1).x, p(i, 1).y, Color.lightGray);
    }

    // Grid
    for (let i = 1; i <= 4; i++) {
      g.setStrokeStyle(Color.toRgb(Color.lightGray))
        .polygon(it.range(cols.length).map((col) => p(col, i / 4)))
        .stroke()
    }

    const path = it.range(cols.length)
      .map((col_i) => p(col_i, !cols[col_i].isNone(row) ? cols[col_i].scale(row) : 0));
    g.setFillStyle("#c0ffc0")
      .polygon(path)
      .fill();
  }

  renderSettings(gc: GridColumn): Element {
    gc.settings ??= getSettings(gc)
    const settings = gc.settings;

    return ui.inputs([
      ui.columnsInput('Radar columns', gc.grid.dataFrame, (columns) => {
        settings.columnNames = names(columns);
        gc.grid.invalidate();
        console.log(JSON.stringify(gc.settings));
      }, {
        available: names(gc.grid.dataFrame.columns.numerical),
        checked: settings?.columnNames ?? names(gc.grid.dataFrame.columns.numerical)
      })
    ]);
  }
}
