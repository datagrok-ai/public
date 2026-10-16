/* The rows of the catalog grid by the model a row shows, read from the grid's status (`text of cell <n> of name`: the
   frame's column names, not the captions). No steps here: the step modules import it. */
import type {Page} from '@playwright/test';
import {type ElementRef, expect, viewers} from '@datagrok-libraries/bdd/runtime';

declare const DG: any;

/** The 1-based row of the catalog grid [target] whose name cell shows [name], or the names it shows. */
export function catalogRow(page: Page, target: ElementRef, name: string): Promise<number | string> {
  return viewers.onViewer(page, target, (el, n) => {
    const values: Record<string, unknown> = DG.Widget.find(el).getWidgetStatus().values ?? {};
    const names: string[] = [];
    for (const [key, value] of Object.entries(values)) {
      const row = /^text of cell (\d+) of name$/.exec(key);
      if (row === null)
        continue;
      if (value === n)
        return Number(row[1]);
      names.push(String(value));
    }
    return names.join(', ') || 'no rows';
  }, name);
}

/** Waits for the row of [name] in [target]; its 1-based number. */
export async function catalogRowOf(page: Page, target: ElementRef, name: string): Promise<number> {
  let row: number | string = '';
  await expect.poll(async () => {
    row = await catalogRow(page, target, name);
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

/** The centre of the area [area] of the row that shows [name], once the grid has settled. */
export async function modelArea(page: Page, target: ElementRef, name: string, area: (row: number) => string):
  Promise<{x: number; y: number}> {
  const row = await catalogRowOf(page, target, name);
  return viewers.centerOf(await viewers.hitArea(page, target, area(row), true));
}
