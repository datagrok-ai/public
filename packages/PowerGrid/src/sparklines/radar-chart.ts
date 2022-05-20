import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {getSettingsBase, names, SummarySettingsBase} from './shared';


class it {
  static range = (n: number) => [...Array(n).keys()];
}


interface RadarChartSettings extends SummarySettingsBase {
  // radius: number;
}

function getSettings(gc: DG.GridColumn): RadarChartSettings {
  return gc.settings ??= {
    ...getSettingsBase(gc),
    // ...{radius: 10,},
  };
}

export class RadarChartCellRender extends DG.GridCellRenderer {
  get name() { return 'radar ts'; }

  get cellType() { return 'radarchart'; }

  // getPreferredCellSize(col: DG.GridColumn) {
  //   return new Size(80,80);
  // }

  get defaultWidth(): number | null { return 80; }

  get defaultHeight(): number | null { return 80; }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    const df = gridCell.grid.dataFrame;

    if (w < 20 || h < 10) return;

    const settings = getSettings(gridCell.gridColumn);
    const box = new DG.Rect(x, y, w, h).fitSquare().inflate(-2, -2);
    const row = gridCell.cell.row.idx;
    const cols = df.columns.byNames(settings.columnNames);

    g.strokeStyle = 'lightgray';

    // axes' point calculator
    const p = (col: number, ratio: number) => new DG.Point(
      box.midX + ratio * box.width * Math.cos(2 * Math.PI * col / (cols.length)) / 2,
      box.midY + ratio * box.width * Math.sin(2 * Math.PI * col / (cols.length)) / 2);

    // axes
    for (let i = 0; i < cols.length; i++)
      g.line(box.midX, box.midY, p(i, 1).x, p(i, 1).y, DG.Color.lightGray);

    // Grid
    for (let i = 1; i <= 4; i++) {
      g.setStrokeStyle(DG.Color.toRgb(DG.Color.lightGray))
        .polygon(it.range(cols.length).map((col) => p(col, i / 4)))
        .stroke();
    }

    const path = it.range(cols.length)
      .map((i) => p(i, !cols[i].isNone(row) ? cols[i].scale(row) : 0));
    g.setFillStyle('#c0ffc0')
      .polygon(path)
      .fill();
  }

  renderSettings(gc: DG.GridColumn): Element {
    gc.settings ??= getSettings(gc);
    const settings = gc.settings;

    return ui.inputs([
      ui.columnsInput('Radar columns', gc.grid.dataFrame, (columns) => {
        settings.columnNames = names(columns);
        gc.grid.invalidate();
      }, {
        available: names(gc.grid.dataFrame.columns.numerical),
        checked: settings?.columnNames ?? names(gc.grid.dataFrame.columns.numerical),
      }),
    ]);
  }
}
