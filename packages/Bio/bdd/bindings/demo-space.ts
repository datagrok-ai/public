// Steps of the space demo clips (features/demos/16b): a selection drawn around a group of a scatter plot's markers.
import type {Page} from '@playwright/test';
import {When} from '@datagrok-libraries/bdd';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';
import {gestures, viewers} from '@datagrok-libraries/bdd/runtime';

type Box = Awaited<ReturnType<typeof viewers.hitArea>>;

/** The first of the points the plot itself takes a press at: inside its view, clear of the strips by its axes where
 * their range sliders show (10 px either side of an axis; the slider takes the press), and with the plot's own canvas
 * on top there (no legend over it). -1 when none is. */
async function firstPressable(page: Page, target: ElementRef, points: [number, number][]): Promise<number> {
  const areas = await viewers.hitAreas(page, target);
  const grown = (b: Box | undefined, dx: number, dy: number): Box | undefined =>
    b && {x: b.x - dx, y: b.y - dy, width: b.width + 2 * dx, height: b.height + 2 * dy};
  const inside = (b: Box | undefined, [x, y]: [number, number]): boolean =>
    !!b && x >= b.x && x <= b.x + b.width && y >= b.y && y <= b.y + b.height;
  const sliders = [grown(areas['y axis'], 11, 0), grown(areas['x axis'], 0, 11)];
  const onCanvas: boolean[] = await viewers.onViewer(page, target, (el, pts: [number, number][]) => pts.map(([x, y]) => {
    const hit = document.elementFromPoint(x, y);
    return hit instanceof HTMLCanvasElement && el.contains(hit);
  }), points);
  return points.findIndex((p, i) => inside(areas['view'], p) && !sliders.some((s) => inside(s, p)) && onCanvas[i]);
}

/** A person selects a group of points by drawing a rectangle around it with Shift held. The rectangle is the one
 * around the markers of those rows (the plot's `marker of row <n>` areas, rows from 1, the first 20 only), a little
 * larger; the plot selects every marker inside it, those of other rows too. The drag starts at a corner the plot
 * takes a press at (`firstPressable`), the top left one where it can, and ends at the opposite one: the plot follows
 * a drag over the whole page, so the end may lie anywhere. */
export const boxAroundMarkers = When('user drags a selection box around the markers of rows {int} to {int} of {widget}',
  async (page: Page, from: number, to: number, target: ElementRef) => {
    if (from < 1 || to < from)
      throw new Error(`rows ${from} to ${to}: a run of rows from 1`);
    const boxes = [];
    for (let row = from; row <= to; row++)
      boxes.push(await viewers.hitArea(page, target, `marker of row ${row}`, row === from));
    const pad = 12;
    const left = Math.min(...boxes.map((b) => b.x)) - pad;
    const top = Math.min(...boxes.map((b) => b.y)) - pad;
    const right = Math.max(...boxes.map((b) => b.x + b.width)) + pad;
    const bottom = Math.max(...boxes.map((b) => b.y + b.height)) + pad;
    // opposite corners pair up as 0-1 and 2-3
    const corners: [number, number][] = [[left, top], [right, bottom], [right, top], [left, bottom]];
    const start = await firstPressable(page, target, corners);
    if (start < 0) {
      throw new Error(`${target.phrase} takes a press at no corner of the box around the markers of rows ${from} to ${to} ` +
        `(${corners.map(([x, y]) => `${Math.round(x)},${Math.round(y)}`).join('; ')}): each is outside its view, by an axis or under something else`);
    }
    const [x0, y0] = corners[start];
    const [x1, y1] = corners[start ^ 1];
    await gestures.withKeys(page, ['Shift'], async () => {
      await page.mouse.move(x0, y0);
      await page.mouse.down();
      await page.mouse.move((x0 + x1) / 2, (y0 + y1) / 2);
      await page.mouse.move(x1, y1, {steps: 3});
      await page.mouse.up();
    });
    await viewers.settleAll(page);
  }, {tier: 'ui', description: 'a Shift-drag between opposite corners of the box 12 px around those markers, from the top left one unless the plot takes no press there (outside its view, by an axis where its range slider shows, or under its legend): the platform\'s rectangle selection'});
