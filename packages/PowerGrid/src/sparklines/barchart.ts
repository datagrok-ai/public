import * as DG from 'datagrok-api/dg';
import {InputBase, Property, TYPE} from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {GridCell, GridColumn, Point, Rect} from "datagrok-api/src/grid";
import {Paint} from "datagrok-api/src/utils";
import {Color} from "datagrok-api/src/widgets";
import {MARKER_TYPE} from "datagrok-api/src/const";
import {getSettingsBase, names, SummarySettingsBase} from "./base";
import "../rect-extensions";

interface BarchartSettings extends SummarySettingsBase {
  // normalize: boolean;
}

function getSettings(gc: DG.GridColumn): BarchartSettings {
  return gc.settings ??= {
    ...getSettingsBase(gc),
    // ...{normalize: true},
  }
}

export class BarchartCellRenderer extends DG.GridCellRenderer {
  get name() { return 'bar ts';}

  get cellType() { return 'barchart_ts';}

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

    // g.fillRect(g, bb, Color.lightGray);
    g.fillStyle = '#8080ff';

    for (let col_i = 0; col_i < cols.length; col_i++) {
      if (!cols[col_i].isNone(row)) {
        const bb = b
          .getLeftPart(cols.length, col_i)
          .getBottomScaled(cols[col_i].scale(row))
          .inflateRel(0.9, 1);

        g.fillRect(bb.left, bb.top, bb.width, bb.height);
      }
    }

    try {
      // throw new Error("test");
    } catch (e) {
      g.beginPath();
      g.strokeStyle = 'red';
      g.moveTo(b.left, b.top);
      g.lineTo(b.right, b.bottom);
      g.moveTo(b.right, b.top);
      g.lineTo(b.left, b.bottom);
      g.stroke();
    }
  }

  renderSettings(gc: GridColumn): Element {
    gc.settings ??= getSettings(gc)
    const settings = gc.settings;

    return ui.inputs([
      ui.columnsInput('Barchart columns', gc.grid.dataFrame, (columns) => {
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
