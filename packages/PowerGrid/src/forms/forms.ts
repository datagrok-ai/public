import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {getSettingsBase, SparklineType, SummarySettingsBase} from '../sparklines/shared';
import {GridColumn} from "datagrok-api/src/grid";
import {LabelElement, Scene} from "./scene";


interface FormSettings extends SummarySettingsBase {
  colorCode: boolean;
}

function getSettings(gc: DG.GridColumn): FormSettings {
  gc.settings ??= getSettingsBase(gc);
  gc.settings.colorCode ??= true;
  return gc.settings;
}


/** Returns approximate length in characters of the longest value in the column. */
function getMaxValueWidth(column: DG.Column): number {
  if (column.type == DG.TYPE.INT)
    return column.max.toString().length;
  else if (column.type == DG.TYPE.FLOAT || column.type == DG.TYPE.QNUM)
    return 50;
  else if (column.type == DG.TYPE.STRING) {
    const values = column.categories;
    if (values.length < 50)
      return Math.min(...values.map(v => v.length));
    return 100;
  }
  return 100;
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
    const row = gridCell.cell.row.idx;
    let cols = df.columns.byNames(settings.columnNames);
    let b = new DG.Rect(x, y, w, h).inflate(-2, -2);

    const scene = new Scene(b);
    g.font = gridCell.grid.props.defaultCellFont;

    // molecules first
    const molCol = cols.find(c => c.semType == DG.SEMTYPE.MOLECULE);
    if (molCol != null && b.width > 50 && b.height > 50) {
      const r = b.width / b.height > 1.5 ? b.getLeftScaled(0.5) : b.getTopScaled(0.5);
      b = b.width / b.height > 1.5 ? b.getRightScaled(0.5) : b.getBottomScaled(0.5);
      cols = cols.filter(c => c.semType !== DG.SEMTYPE.MOLECULE);
      let cell = gridCell.grid.cell(molCol.name, gridCell.gridRow);
      cell.renderer.render(g, r.x, r.y, r.width, r.height, cell, cell.style);
    }

    const maxNameWidth = Math.min(200, Math.max(...cols.map(c => g.measureText(c.name).width)));
    const maxValueWidth = Math.min(100, Math.max(...cols.map(c => getMaxValueWidth(c) * 8)));
    const showColumnNames = maxNameWidth + maxValueWidth + 10 < b.width;
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
        //scene.elements.push(new LabelElement(r, col.getString(row)));
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
        const r = new DG.Rect(b.x, b.y + i * 20, b.width, 20);
        // render in a column
        //g.fillStyle = 'darkgrey';
        if (showColumnNames)
          scene.elements.push(new LabelElement(r.getLeft(columnNamesWidth), col.name, {horzAlign: 'right'}));
          //g.fillText(col.name, b.x + columnNamesWidth - g.measureText(col.name).width - 5, b.y + i * 20);
        g.fillStyle = valueColor;

        if (col.semType == DG.SEMTYPE.MOLECULE) {
          let cell = gridCell.grid.cell(col.name, gridCell.gridRow);
          cell.renderer.render(g, b.x + columnNamesWidth, r.y, 50, 30, cell, cell.style);
        }
        else
          g.fillText(col.getString(row), b.x + columnNamesWidth, r.y);
      }
    }

    scene.render(g);
  }
}