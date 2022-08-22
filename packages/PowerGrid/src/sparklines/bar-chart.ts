import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {getSettingsBase, names, SummarySettingsBase} from './shared';
import {createTooltip, Hit} from './helper';


interface BarChartSettings extends SummarySettingsBase {
  // normalize: boolean;
  minH: number;
  colorCode: boolean;
}

function getSettings(gc: DG.GridColumn): BarChartSettings {
  return gc.settings ??= {
    ...getSettingsBase(gc),
    ...{minH: 0.05},
    ...{colorCode: false},
    // ...{normalize: true},
  };
}

function onHit(gridCell: DG.GridCell, e: MouseEvent | any): Hit {
  const settings = getSettings(gridCell.gridColumn);
  const df = gridCell.grid.dataFrame;
  const cols = df.columns.byNames(settings.columnNames);
  const row = gridCell.cell.row.idx;
  const b = new DG.Rect(gridCell.bounds.x, gridCell.bounds.y, gridCell.bounds.width, gridCell.bounds.height).inflate(-2, -2);
  const width = b.width / cols.length;
  const activeColumn = Math.floor((e.layerX - b.left) / width);
  let answer: Hit = {
    isHit: false,
    activeColumn: activeColumn,
    row: row,
    cols: cols,
  };
  if ((activeColumn > cols.length) || (activeColumn < 0)) {
    return answer;
  }
  const bb = b
    .getLeftPart(cols.length, activeColumn)
    .getBottomScaled(cols[activeColumn].scale(row) > settings.minH ? cols[activeColumn].scale(row) : settings.minH)
    .inflateRel(0.9, 1);
  answer.isHit = (e.layerY >= bb.top);
  return answer;
}


export class BarChartCellRenderer extends DG.GridCellRenderer {
  get name() { return 'bar ts'; }

  get cellType() { return 'barchart'; }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent | any): void {
    const hitData = onHit(gridCell, e);
    if (hitData.isHit) {
      ui.tooltip.show(ui.divV(createTooltip(hitData.cols, hitData.activeColumn, hitData.row)), e.x + 16, e.y + 16);
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


    for (let i = 0; i < cols.length; i++) {
      if (!cols[i].isNone(row)) {
        let color = settings.colorCode ? DG.Color.getCategoricalColor(i) : DG.Color.blue;
        g.setFillStyle(DG.Color.toRgb(color));
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

    const colorCodeScaleProp = DG.Property.js('colorCode', DG.TYPE.BOOL, {
      description: 'Activates color rendering'
    });

    const colorCodeNormalizeInput = DG.InputBase.forProperty(colorCodeScaleProp, settings);
    colorCodeNormalizeInput.onChanged(() => { gc.grid.invalidate(); });

    return ui.inputs([
      ui.columnsInput('Barchart columns', gc.grid.dataFrame, (columns) => {
        settings.columnNames = names(columns);
        gc.grid.invalidate();
      }, {
        available: names(gc.grid.dataFrame.columns.numerical),
        checked: settings?.columnNames ?? names(gc.grid.dataFrame.columns.numerical),
      }),
      colorCodeNormalizeInput
    ]);
  }
}
