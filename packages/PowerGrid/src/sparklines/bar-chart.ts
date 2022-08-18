import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {getSettingsBase, names, SummarySettingsBase} from './shared';


interface BarChartSettings extends SummarySettingsBase {
  // normalize: boolean;
}

function getSettings(gc: DG.GridColumn): BarChartSettings {
  return gc.settings ??= {
    ...getSettingsBase(gc),
    ...{minH: 0.05},
    // ...{normalize: true},
  };
}

export class BarChartCellRenderer extends DG.GridCellRenderer {
  get name() { return 'bar ts'; }

  get cellType() { return 'barchart'; }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent | any): void {
    const settings: any = getSettings(gridCell.gridColumn);
    const df = gridCell.grid.dataFrame;
    const cols = df.columns.byNames(settings.columnNames);
    const row = gridCell.cell.row.idx;
    const b = new DG.Rect(gridCell.bounds.x, gridCell.bounds.y, gridCell.bounds.width, gridCell.bounds.height).inflate(-2, -2);
    const width = b.width / cols.length;
    const activeColumn = Math.floor((e.layerX - b.left) / width);
    if ((activeColumn > cols.length) || (activeColumn < 0)) {
      ui.tooltip.hide();
      return;
    }
    const bb = b
      .getLeftPart(cols.length, activeColumn)
      .getBottomScaled(cols[activeColumn].scale(row) > settings.minH ? cols[activeColumn].scale(row) : settings.minH)
      .inflateRel(0.9, 1);
    const maxHeight = b.bottom - b.top;
    if (e.layerY >= bb.top) {
      let arr = [];
      // create tooltip data
      for (let i = 0; i < cols.length; i++) {
        arr.push(ui.divH([ui.divText(`${cols[i].name}:`, {
              style: {
                margin: '0 10px 0 0',
                fontWeight: (activeColumn == i) ? 'bold' : 'normal',
              }
            }), ui.divText(`${Math.floor(cols[i].get(row) * 100) / 100}`, {
              style: {
                fontWeight: (activeColumn == i) ? 'bold' : 'normal',
              }
            })]
          )
        );
      }
      ui.tooltip.show(ui.divV(arr), e.x + 16, e.y + 16);
    } else {
      ui.tooltip.hide();
    }

  }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    const df = gridCell.grid.dataFrame;

    if (w < 20 || h < 10 || df === void 0) return;

    const settings: any = getSettings(gridCell.gridColumn);
    const b = new DG.Rect(x, y, w, h).inflate(-2, -2);
    const row = gridCell.cell.row.idx;
    const cols = df.columns.byNames(settings.columnNames);

    // g.fillRect(g, bb, Color.lightGray);
    g.fillStyle = '#8080ff';

    for (let i = 0; i < cols.length; i++) {
      if (!cols[i].isNone(row)) {
        const bb = b
          .getLeftPart(cols.length, i)
          .getBottomScaled(cols[i].scale(row) > settings.minH ? cols[i].scale(row) : settings.minH)
          .inflateRel(0.9, 1);

        g.fillRect(bb.left, bb.top, bb.width, bb.height);
      }
    }
  }

  renderSettings(gc: DG.GridColumn): Element {
    gc.settings ??= getSettings(gc);
    const settings = gc.settings;

    return ui.inputs([
      ui.columnsInput('Barchart columns', gc.grid.dataFrame, (columns) => {
        settings.columnNames = names(columns);
        gc.grid.invalidate();
      }, {
        available: names(gc.grid.dataFrame.columns.numerical),
        checked: settings?.columnNames ?? names(gc.grid.dataFrame.columns.numerical),
      }),
    ]);
  }
}
