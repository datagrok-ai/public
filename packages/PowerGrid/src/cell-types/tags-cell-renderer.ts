import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {names} from '../sparklines/shared';


export const TAGS_CELL_TYPE = 'Tags';

export interface TagsColumnSettings {
  columnNames: string[];
}

type TagMode = 'full' | 'short' | 'min';

interface TagLayoutItem {
  tag: string;
  display: string;
  x: number;
  y: number;
  width: number;
  height: number;
  mode: TagMode;
}

const _measures: {[key: string]: number} = {};

function measure(g: CanvasRenderingContext2D, tag: string): number {
  return _measures[tag] ??= g.measureText(tag).width;
}

/**
 * Heuristic shortening: cut at the first separator, else first 3 chars + ellipsis.
 * @param {string} tag full tag string
 * @return {string} shortened representation (may equal the input if already short)
 */
function shorten(tag: string): string {
  const m = tag.match(/^([^\s_\-:.|/\\,;])([^\s_\-:.|/\\,;]*?)[\s_\-:.|/\\,;]/);
  if (m) {
    const head = m[1] + m[2];
    if (head.length > 0 && head.length < Math.min(tag.length - 1, 6)) return head + '…';
  }
  if (tag.length > 4) return tag.slice(0, 3) + '…';
  return tag;
}

/**
 * Renders a comma-separated string value as checkboxes with options retrieved
 * from the column's `.choices` tag.
 * */
export class TagsCellRenderer extends DG.GridCellRenderer {
  private _hoverLayout: TagLayoutItem[] | null = null;

  get name() { return TAGS_CELL_TYPE; }

  get cellType() { return TAGS_CELL_TYPE; }

  private getValues(gridCell: DG.GridCell): string[] {
    const settings = this.getSettings(gridCell.gridColumn);
    const customSeparator = gridCell.tableColumn?.getTag(DG.Tags.MultiValueSeparator);
    return settings ?
      settings.columnNames.filter((colName) =>
        gridCell.tableRow?.table.col(colName) != null && gridCell.tableRow?.get(colName)) :
      // Empty tokens are dropped: an empty cell splits into [''] and used to
      // paint a zero-width chip that reads as data where there is none (and a
      // trailing separator would add one more).
      (gridCell.cell.valueString?.split(customSeparator ?? ',')
        .map((s) => s.trim()).filter((s) => s.length > 0) ?? []);
  }

  /**
   * Computes tag positions for the given cell bounds. Tries three render modes
   * per tag: full text, heuristically shortened text, or a minimized marker.
   * Optimizes by maximizing the number of full tags first, then the number of
   * shortened tags among the remaining. Positions are in cell-local coordinates.
   * @param {CanvasRenderingContext2D} g canvas context used for text measurement (font is taken from g)
   * @param {number} w available cell width in pixels
   * @param {number} h available cell height in pixels
   * @param {string[]} values tag strings in render order
   * @return {TagLayoutItem[]} ordered list of placed tags with their bounds and mode
   */
  calculateLayout(g: CanvasRenderingContext2D, w: number, h: number, values: string[]): TagLayoutItem[] {
    if (values.length === 0)
      return [];

    const n = values.length;
    const dy = h > 16 ? 5 : 2;
    const tagHeight = Math.min(h - 2 * dy, 17);
    const rowStride = 19;
    const maxRows = Math.max(1, Math.floor((h - dy) / rowStride) + (((h - dy) % rowStride >= tagHeight) ? 1 : 0));
    const minVisibleWidth = 4;
    const fullPadding = 4;
    const fullGap = 4;
    const minGap = 3;

    const shortText = values.map((t) => shorten(t));

    const tryFit = (a: number, b: number): TagLayoutItem[] | null => {
      const items: TagLayoutItem[] = [];
      let cx = 2;
      let cy = dy;
      let row = 1;

      const place = (i: number, mode: TagMode): boolean => {
        const display = mode === 'full' ? values[i] : mode === 'short' ? shortText[i] : '';
        const tw = mode === 'min' ? minVisibleWidth : measure(g, display) + fullPadding;
        if (tw + 2 > w) return false;
        if (cx > 2 && cx + tw > w) {
          cx = 2;
          cy += rowStride;
          row++;
          if (row > maxRows) return false;
        }
        items.push({tag: values[i], display, x: cx, y: cy, width: tw, height: tagHeight, mode});
        cx += tw + (mode === 'min' ? minGap : fullGap);
        return true;
      };

      for (let i = 0; i < a; i++)
        if (!place(i, 'full')) return null;
      for (let i = a; i < b; i++)
        if (!place(i, 'short')) return null;
      for (let i = b; i < n; i++)
        if (!place(i, 'min')) return null;
      return items;
    };

    // Minimize minimized chips first (more informative), then maximize full count.
    for (let minCount = 0; minCount <= n; minCount++) {
      const b = n - minCount;
      for (let a = b; a >= 0; a--) {
        const result = tryFit(a, b);
        if (result) return result;
      }
    }

    // Fallback: place as many minimized markers as fit, drop the rest.
    const items: TagLayoutItem[] = [];
    let cx = 2;
    let cy = dy;
    let row = 1;
    for (let i = 0; i < n; i++) {
      if (cx > 2 && cx + minVisibleWidth > w) {
        cx = 2;
        cy += rowStride;
        row++;
        if (row > maxRows) break;
      }
      items.push({tag: values[i], display: '', x: cx, y: cy, width: minVisibleWidth, height: tagHeight, mode: 'min'});
      cx += minVisibleWidth + minGap;
    }
    return items;
  }

  /** Tag → index in the exact list the platform's categorical color editor shows
   * (split by the multi-value separator, trimmed, deduped, sorted when the separator
   * tag is set; raw category order otherwise), so default colors match the editor.
   * Cached by column version. */
  private tagIndexes(gridColumn: DG.GridColumn, col: DG.Column): Map<string, number> {
    const cached = gridColumn.temp['tagIndexes'];
    if (cached?.version === col.version)
      return cached.indexes;
    let sep = col.getTag(DG.Tags.MultiValueSeparator);
    let values: string[];
    if (!sep)
      values = col.categories;
    else {
      if (sep === '\\n') sep = '\n';
      const set = new Set<string>();
      for (const cat of col.categories)
        for (const s of cat.split(sep).map((x: string) => x.trim()))
          if (s) set.add(s);
      values = [...set].sort();
    }
    const indexes = new Map(values.map((v, i) => [v, i]));
    gridColumn.temp['tagIndexes'] = {version: col.version, indexes};
    return indexes;
  }

  private getColor(gridCell: DG.GridCell, tag: string): number {
    const gridColumn = gridCell.gridColumn;
    const defined: number | string | undefined = gridColumn.categoryColors?.[tag];
    if (defined != null)
      return typeof defined === 'string' ? DG.Color.fromHtml(defined) : defined;
    const col = gridCell.tableColumn;
    const idx = col ? this.tagIndexes(gridColumn, col).get(tag) : undefined;
    if (idx !== undefined)
      return DG.Color.getCategoricalColor(idx);
    const colors = gridColumn.temp['catColors'] ??= {};
    return colors[tag] ??= DG.Color.getCategoricalColor(Object.keys(colors).length);
  }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    // somehow dart comes null here, will need to investigate it and fix it, now just a workaround
    if (!gridCell.gridColumn.dart)
      return;

    g.font = cellStyle.font;
    g.textBaseline = 'top';
    g.textAlign = 'left';

    const values = this.getValues(gridCell);
    const layout = this.calculateLayout(g, w, h, values);

    for (const item of layout) {
      const color = this.getColor(gridCell, item.tag);
      g.fillStyle = DG.Color.toHtml(color);
      g.roundRect(x + item.x, y + item.y, item.width, item.height, 4);
      g.fill();

      if (item.mode !== 'min') {
        g.fillStyle = DG.Color.toHtml(DG.Color.getContrastColor(color));
        g.fillText(item.display, x + item.x + 2, y + item.y + 3);
      }
    }
  }

  private computeHoverLayout(gridCell: DG.GridCell): TagLayoutItem[] {
    const canvas = gridCell.grid.canvas;
    const g = canvas.getContext('2d')!;
    const b = gridCell.bounds;
    return this.calculateLayout(g, b.width, b.height, this.getValues(gridCell));
  }

  onMouseEnter(gridCell: DG.GridCell, _e: MouseEvent): void {
    if (!gridCell.gridColumn.dart) return;
    this._hoverLayout = this.computeHoverLayout(gridCell);
  }

  onMouseLeave(_gridCell: DG.GridCell, _e: MouseEvent): void {
    this._hoverLayout = null;
    ui.tooltip.hide();
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    if (!gridCell.gridColumn.dart) return;
    if (!this._hoverLayout)
      this._hoverLayout = this.computeHoverLayout(gridCell);

    const b = gridCell.bounds;
    const lx = e.offsetX - b.left;
    const ly = e.offsetY - b.top;

    for (const item of this._hoverLayout) {
      if (lx >= item.x && lx <= item.x + item.width && ly >= item.y && ly <= item.y + item.height) {
        if (item.mode !== 'full')
          ui.tooltip.show(ui.divText(item.tag, {style: {fontSize: '14px', fontWeight: 'bold'}}), e.x + 16, e.y + 16);
        else
          ui.tooltip.hide();
        return;
      }
    }
    ui.tooltip.hide();
  }

  getSettings(gridColumn: DG.GridColumn): TagsColumnSettings | undefined {
    // if we have values contained in this column as comma-separated list, we don't need settings
    if (gridColumn.column)
      return undefined;

    // backwards compat: migrate old flat settings format to nested under cell type key
    if (gridColumn.settings?.columnNames != undefined && gridColumn.settings[TAGS_CELL_TYPE] == undefined)
      gridColumn.settings = {[TAGS_CELL_TYPE]: {columnNames: gridColumn.settings.columnNames}};

    return gridColumn?.settings?.[TAGS_CELL_TYPE]?.columnNames !== undefined ?
      gridColumn.settings[TAGS_CELL_TYPE] :
      (gridColumn.settings ??= {})[TAGS_CELL_TYPE] = {
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
