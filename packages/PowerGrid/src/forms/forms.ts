import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {MARKER_TYPE} from 'datagrok-api/dg';
import {getSettingsBase, SparklineType, SummarySettingsBase} from '../sparklines/shared';
import {GridCell, GridColumn} from "datagrok-api/src/grid";
import {GridCellElement, LabelElement, MarkerElement, Scene} from "./scene";


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

  static makeScene(gridCell: DG.GridCell, b?: DG.Rect): Scene {

    const g = gridCell.grid.canvas.getContext("2d")!;
    const df = gridCell.grid.dataFrame;
    const settings: FormSettings = getSettings(gridCell.gridColumn);
    const row = gridCell.cell.row.idx;
    let cols = df.columns.byNames(settings.columnNames);
    b ??= gridCell.bounds.inflate(-2, -2);
    const scene = new Scene(b);

    // molecules first
    const molCol = cols.find(c => c.semType == DG.SEMTYPE.MOLECULE);
    if (molCol != null && b.width > 50 && b.height > 50) {
      const r = b.width / b.height > 1.5 ? b.getLeftScaled(0.5) : b.getTopScaled(0.5);
      b = b.width / b.height > 1.5 ? b.getRightScaled(0.5) : b.getBottomScaled(0.5);
      cols = cols.filter(c => c.semType !== DG.SEMTYPE.MOLECULE);
      let cell = gridCell.grid.cell(molCol.name, gridCell.gridRow);
      scene.elements.push(new GridCellElement(r, cell));

      //cell.renderer.render(g, r.x, r.y, r.width, r.height, cell, cell.style);
    }

    const maxNameWidth = Math.min(200, Math.max(...cols.map(c => g.measureText(c.name).width)));
    const maxValueWidth = Math.min(100, Math.max(...cols.map(c => getMaxValueWidth(c) * 8)));
    const showColumnNames = maxNameWidth + maxValueWidth + 10 < b.width;
    const columnNamesWidth = showColumnNames ? maxNameWidth + 10 : 0;
    const colHeight = Math.min(20, b.height / cols.length);

    for (let i = 0; i < cols.length; i++) {
      const col = cols[i];
      const cell = gridCell.grid.cell(col.name, gridCell.gridRow);
      let numColor = col.meta.colors.getType() == DG.COLOR_CODING_TYPE.OFF ? gridCell.grid.props.cellTextColor : col.meta.colors.getColor(row);
      const valueColor = DG.Color.toHtml(numColor);

      if (b.height < cols.length * 20 * 0.7) {
        // render in one row
        const r = b.getLeftPart(cols.length, i);
        const e = new GridCellElement(r, cell);
        // const e = r.width > 30
        //   ? new LabelElement(r, col.getString(row), { color: valueColor, horzAlign: 'center', vertAlign: 'center' })
        //   : new MarkerElement(new DG.Rect(r.midX - 5, r.midY - 5, 10, 10), MARKER_TYPE.CIRCLE, numColor);
        scene.elements.push(e);
      }
      else {
        // render in a column
        const r = new DG.Rect(b.x, b.y + i * colHeight, b.width, colHeight);
        if (showColumnNames)
          scene.elements.push(new LabelElement(r.getLeft(columnNamesWidth), col.name, {horzAlign: 'right'}));
        scene.elements.push(new LabelElement(r.cutLeft(columnNamesWidth).move(5, 0), col.getString(row), {horzAlign: 'left'}));
      }
    }

    return scene;
  }

  render(g: CanvasRenderingContext2D,
         x: number, y: number, w: number, h: number,
         gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    const scene = FormCellRenderer.makeScene(gridCell);
    scene.render(g);
  }

  makeBestScene(gridCell: DG.GridCell): Scene {
    const height = getSettings(gridCell.gridColumn).columnNames
      .map(name => gridCell.tableColumn?.semType === DG.SEMTYPE.MOLECULE ? 150 : 20)
      //.map(name => gridCell.grid.cell(name, gridCell.gridRow).renderer.getDefaultSize(gridCell.gridColumn).height)
      .reduce((sum, height) => sum! + height!, 0)!;
    return FormCellRenderer
      .makeScene(gridCell, new DG.Rect(0, 0, 200, height));
  }

  onMouseEnter(gridCell: DG.GridCell, e: MouseEvent): void {
    const scene = FormCellRenderer.makeScene(gridCell);
    if (scene.elements.every(e => e.bounds.width < 25 && e.bounds.height < 25)) {
      setTimeout(() => ui.tooltip.show(this.makeBestScene(gridCell).toCanvas(), gridCell.documentBounds.right, gridCell.documentBounds.top), 200);
    }
  }

  onMouseDown(gridCell: GridCell, e: MouseEvent): void {
    grok.shell.o = this.makeBestScene(gridCell).toCanvas();
  }
}