import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {
  getSettingsBase,
  names,
  SparklineType,
  SummarySettingsBase,
  createTooltip,
  Hit
} from '../sparklines/shared';
import {GridColumn} from "datagrok-api/src/grid";


interface FormSettings extends SummarySettingsBase {
  colorCode: boolean;
}

function getSettings(gc: DG.GridColumn): FormSettings {
  gc.settings ??= getSettingsBase(gc);
  gc.settings.colorCode ??= true;
  return gc.settings;
}


export class FormCellRenderer extends DG.GridCellRenderer {
  get name() { return SparklineType.Form; }

  get cellType() { return SparklineType.Form; }

  getDefaultSize(gridColumn: GridColumn): { width?: number | null; height?: number | null; } {
    return {
      width: 200,
      height: getSettings(gridColumn).columnNames.length * 20
    };
  }

  render(g: CanvasRenderingContext2D,
         x: number, y: number, w: number, h: number,
         gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    const df = gridCell.grid.dataFrame;
    const settings: FormSettings = getSettings(gridCell.gridColumn);
    const b = new DG.Rect(x, y, w, h).inflate(-2, -2);
    const row = gridCell.cell.row.idx;
    const cols = df.columns.byNames(settings.columnNames);

    g.font = gridCell.grid.props.defaultCellFont;

    const maxNameWidth = Math.max(...settings.columnNames.map((name) => g.measureText(name).width));
    const showColumnNames = maxNameWidth < b.width / 2;
    const columnNamesWidth = showColumnNames ? maxNameWidth + 10 : 0;

    g.textBaseline = 'top';
    g.textAlign = 'left';

    for (let i = 0; i < cols.length; i++) {
      const col = cols[i];
      const numColor = col.meta.colors.getType() == DG.COLOR_CODING_TYPE.OFF ? gridCell.grid.props.cellTextColor : col.meta.colors.getColor(row);
      const valueColor = DG.Color.toHtml(numColor);

      if (h < cols.length * 20) {
        // render in one row
        const r = b.getLeftPart(cols.length, i);
        g.save();
        g.rect(r.x, r.y, r.width, r.height);
        g.clip();
        g.fillStyle = valueColor;
        if (r.width > 30)
          g.fillText(col.getString(row), r.x + 2, r.y + 2);
        else
          DG.Paint.marker(g, DG.MARKER_TYPE.CIRCLE, r.x + r.width / 2, r.y + r.height / 2, numColor, 5);
        g.restore();
      }
      else {
        // render in a column
        g.fillStyle = 'darkgrey';
        if (showColumnNames)
          g.fillText(col.name, b.x + columnNamesWidth - g.measureText(col.name).width - 5, b.y + i * 20);
        g.fillStyle = valueColor;
        g.fillText(col.getString(row), b.x + columnNamesWidth, b.y + i * 20);
      }
    }
  }
}