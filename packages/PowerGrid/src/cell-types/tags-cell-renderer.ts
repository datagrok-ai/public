import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from "datagrok-api/ui";
import {names} from '../sparklines/shared';


export const TAGS_CELL_TYPE = 'Tags';

export interface TagsColumnSettings {
  columnNames: string[];
}

const _measures: {[key: string]: number} = {};

function measure(g: CanvasRenderingContext2D, tag: string): number {
  return _measures[tag] ??= g.measureText(tag).width;
}

/**
 * Renders a comma-separated string value as checkboxes with options retrieved
 * from the column's `.choices` tag.
 * */
export class TagsCellRenderer extends DG.GridCellRenderer {
  get name() { return TAGS_CELL_TYPE; }

  get cellType() { return TAGS_CELL_TYPE; }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    // somehow dart comes null here, will need to investigate it and fix it, now just a workaround
    if (!gridCell.gridColumn.dart)
      return;

    const settings = this.getSettings(gridCell.gridColumn);
    const values: string[]
      = settings ? settings.columnNames.filter((colName) =>
        gridCell.tableRow?.table.col(colName) != null && gridCell.tableRow?.get(colName))
      : gridCell.cell.valueString.split(',').map((s) => s.trim());

    let cx = 2;
    let dy = h > 16 ? 5 : 2
    let cy = dy;
    g.font = cellStyle.font;
    g.textBaseline = 'top';
    g.textAlign = 'left';

    const getColor = (tag: string) => {
      const colors = gridCell.gridColumn.temp['catColors'] ??= {};
      return colors[tag] ??= DG.Color.getCategoricalColor(Object.keys(colors).length);
    };

    const len = values.map((tag) => measure(g, tag) + 3).reduce((total, num) => total + num, 0);
    const fits = len < w * Math.round(h / 20);

    for (const tag of values) {
      const width = fits ? measure(g, tag) : 4;
      const drawTag = () => {
        const color = getColor(tag);
        g.fillStyle = DG.Color.toHtml(color);
        g.roundRect(x + cx, y + cy, width + 4, Math.min(h - 2 * dy, 17), 4);
        g.fill();

        if (fits) {
          g.fillStyle = DG.Color.toHtml(DG.Color.getContrastColor(color));
          g.fillText(tag, x + cx + 2, y + cy + 1);
        }
      };

      if (cx + width <= w || cy + 16 > h)
        drawTag();
      else {
        cx = 2;
        cy += 19;
        drawTag();
      }

      cx += width + 8;

    }
  }

  getSettings(gridColumn: DG.GridColumn): TagsColumnSettings | undefined {
    // if we have values contained in this column as comma-separated list, we don't need settings
    if (gridColumn.column)
      return undefined;

    // backwards compat: migrate old flat settings format to nested under cell type key
    if (gridColumn.settings?.columnNames != undefined && gridColumn.settings[TAGS_CELL_TYPE] == undefined)
      gridColumn.settings = {[TAGS_CELL_TYPE]: {columnNames: gridColumn.settings.columnNames}};

    return gridColumn?.settings?.[TAGS_CELL_TYPE]?.columnNames !== undefined
      ? gridColumn.settings[TAGS_CELL_TYPE]
      : (gridColumn.settings ??= {})[TAGS_CELL_TYPE] = {
        columnNames: [...gridColumn.grid.dataFrame.columns.boolean].map((c) => c.name)
      };
  }

  renderSettings(gridColumn: DG.GridColumn): HTMLElement | null {
    const settings = this.getSettings(gridColumn);
    if (!settings)
      return null;

    return ui.inputs([
      ui.input.columns('Columns', {
        value: gridColumn.grid.dataFrame.columns.byNames(settings.columnNames),
        table: gridColumn.grid.dataFrame,
        onValueChanged: (value) => {
          settings.columnNames = names(value);
          gridColumn.grid.invalidate();
        },
        available: names(gridColumn.grid.dataFrame.columns.boolean)
      })
    ]);
  }
}
