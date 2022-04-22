import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {GridCell, GridColumn, Rect} from 'datagrok-api/src/grid';
import {getSettingsBase, names, SummarySettingsBase} from './base';


interface PiechartSettings extends SummarySettingsBase {
  radius: number;
}

function getSettings(gc: DG.GridColumn): PiechartSettings {
  return gc.settings ??= {
    ...getSettingsBase(gc),
    ...{normalize: true},
  };
}

export class PiechartCellRenderer extends DG.GridCellRenderer {
  get name() {return 'pie ts';}

  get cellType() {return 'piechart_ts';}

  // getPreferredCellSize(col: DG.GridColumn) {
  //   return new Size(80,80);
  // }

  get defaultWidth(): number | null {return 80;}

  get defaultHeight(): number | null {return 80;}

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: GridCell, cellStyle: DG.GridCellStyle,
  ) {
    const df = gridCell.grid.dataFrame;

    if (w < 5 || h < 5 || df === void 0) return;

    const settings = getSettings(gridCell.gridColumn);
    const row: number = gridCell.cell.row.idx;
    const cols = df.columns.byNames(settings.columnNames);
    const box = new Rect(x, y, w, h).fitSquare().inflate(-2, -2);

    for (let i = 0; i < cols.length; i++) {
      if (cols[i].isNone(row))
        continue;


      const r = cols[i].scale(row) * box.width / 2;
      g.beginPath();
      g.moveTo(box.midX, box.midY);
      g.arc(box.midX, box.midY, r,
        2 * Math.PI * i / cols.length, 2 * Math.PI * (i + 1) / cols.length);
      g.closePath();

      g.fillStyle = DG.Color.toRgb(DG.Color.getCategoricalColor(i));
      g.fill();
    }

    /*
    part of d4;

    class PieBarChartCellRenderer extends GridCellRenderer {
      static const String CELL_TYPE = "pie bar chart";

      String get cellType => CELL_TYPE;
      Size getPreferredCellSize(GridColumn col) => new Size(80, 80);

      void render(CanvasRenderingContext2D g, Rect bounds, GridCell gridCell, GridCellStyle cellStyle, GridLook look) {
        if (bounds.width < 5 || bounds.height < 5) return;

        var cols = gridCell.cell.dataFrame.columns.numerical.toList();
        var box = bounds.fitSquare().inflate(-2, -2);

        for (int i = 0; i < cols.length; i++) {
          if (cols[i].isNone(gridCell.cell.row))
            continue;

          var r = cols[i].scale(gridCell.cell.row) * box.width / 2;
          g.beginPath();
          g.moveTo(box.midX, box.midY);
          g.arc(box.midX, box.midY, r, 2 * math.PI * i / cols.length, 2 * math.PI * (i + 1) / cols.length);
          g.closePath();
          setFillColor(g, Color.getDefault(i));
          g.fill();
        }
      }
    }

    /**/
  }

  renderSettings(gc: GridColumn): Element {
    gc.settings ??= getSettings(gc);
    const settings = gc.settings;

    return ui.inputs([
      ui.columnsInput('Radar columns', gc.grid.dataFrame, (columns) => {
        settings.columnNames = names(columns);
        gc.grid.invalidate();
        console.log(JSON.stringify(gc.settings));
      }, {
        available: names(gc.grid.dataFrame.columns.numerical),
        checked: settings?.columnNames ?? names(gc.grid.dataFrame.columns.numerical),
      }),
    ]);
  }
}
