/* The rows of the catalog grid by the model a row shows (and of the Methods grid by its method), read from the grid's
   status (`text of cell <n> of name`: the frame's column names, not the captions). No steps here: the step modules
   import it. */
import type {Page} from '@playwright/test';
import {type ElementRef, expect, viewers} from '@datagrok-libraries/bdd/runtime';

declare const DG: any;

/** The 1-based row of the grid [target] whose [column] cell shows [name], or the names it shows. */
export function catalogRow(page: Page, target: ElementRef, name: string, column = 'name'): Promise<number | string> {
  return viewers.onViewer(page, target, (el, [n, c]) => {
    const values: Record<string, unknown> = DG.Widget.find(el).getWidgetStatus().values ?? {};
    const names: string[] = [];
    for (const [key, value] of Object.entries(values)) {
      const row = /^text of cell (\d+) of (.+)$/.exec(key);
      if (row === null || row[2] !== c)
        continue;
      if (value === n)
        return Number(row[1]);
      names.push(String(value));
    }
    return names.join(', ') || 'no rows';
  }, [name, column] as [string, string]);
}

/** Waits for the row of [name] (in its [column] cell) in [target]; its 1-based number. */
export async function catalogRowOf(page: Page, target: ElementRef, name: string, column = 'name'): Promise<number> {
  let row: number | string = '';
  await expect.poll(async () => {
    row = await catalogRow(page, target, name, column);
    return typeof row === 'number' ? 'listed' : `not listed; the rows show ${row}`;
  }, {message: `a row "${name}" in ${target.phrase}`}).toBe('listed');
  return Number(row);
}

export interface DrawnTexts {isFriendlyNameShown: boolean; headers: (string | null)[]; cells: string[]}

/** What the grid [target] draws as text: the caption of each visible column in order (the text of the header cell the
 * grid's hit test finds at its `header <column>` area; null where no header is drawn) and its drawn cells' texts. With
 * showFriendlyName on, the grid draws the name over the caption, so a caption claim needs it off. */
export function drawnTexts(page: Page, target: ElementRef): Promise<DrawnTexts> {
  return viewers.onViewer(page, target, (el) => {
    const grid = DG.Widget.find(el);
    const status = grid.getWidgetStatus();
    const areas: Record<string, {x: number; y: number; width: number; height: number}> = status.hitAreas ?? {};
    const headers: (string | null)[] = [];
    for (let i = 0; i < grid.columns.length; i++) {
      const column = grid.columns.byIndex(i);
      if (!column.visible || !column.column)
        continue;
      const area = areas[`header ${column.name}`];
      const cell = area === undefined ? null : grid.hitTest(area.x + area.width / 2, area.y + area.height / 2);
      headers.push(cell?.isColHeader ? String(cell.customText) : null);
    }
    const cells = Object.entries(status.values ?? {}).filter(([key]) => /^text of cell \d+ of /.test(key))
      .map(([, value]) => String(value));
    return {isFriendlyNameShown: grid.props.showFriendlyName === true, headers, cells};
  }, undefined);
}

/** The drawn text of the [column] cell in the row whose [keyColumn] cell shows [key], from one snapshot of the grid;
 * a sentence naming what is missing otherwise. */
export function cellOfRow(page: Page, target: ElementRef, keyColumn: string, key: string, column: string):
  Promise<string> {
  return viewers.onViewer(page, target, (el, [k, kc, c]) => {
    const values: Record<string, unknown> = DG.Widget.find(el).getWidgetStatus().values ?? {};
    const row = Object.keys(values).map((x) => /^text of cell (\d+) of (.+)$/.exec(x))
      .find((m) => m !== null && m[2] === kc && values[m[0]] === k);
    if (row === undefined || row === null)
      return `no drawn row of "${k}"`;
    const cell = `text of cell ${row[1]} of ${c}`;
    return cell in values ? String(values[cell]) : `(no drawn ${cell})`;
  }, [key, keyColumn, column] as [string, string, string]);
}

/** The centre of the area [area] of the row whose [column] cell shows [name], once the grid has settled. */
export async function modelArea(page: Page, target: ElementRef, name: string, area: (row: number) => string,
  column = 'name'): Promise<{x: number; y: number}> {
  const row = await catalogRowOf(page, target, name, column);
  return viewers.centerOf(await viewers.hitArea(page, target, area(row), true));
}
