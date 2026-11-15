/* Steps of the space demo clips (features/demos/16b): a selection drawn around a group of a scatter plot's markers. */
import type {Page} from '@playwright/test';
import {When} from '@datagrok-libraries/bdd';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';
import {gestures, viewers} from '@datagrok-libraries/bdd/runtime';

/** A person selects a group of points by drawing a rectangle around it with Shift held. The rectangle is the one
 * around the markers of those rows (the plot's `marker of row <n>` areas, rows from 1, the first 20 only), a little
 * larger; the plot selects every marker inside it, those of other rows too. */
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
    await gestures.withKeys(page, ['Shift'], async () => {
      await page.mouse.move(left, top);
      await page.mouse.down();
      await page.mouse.move((left + right) / 2, (top + bottom) / 2);
      await page.mouse.move(right, bottom, {steps: 3});
      await page.mouse.up();
    });
    await viewers.settleAll(page);
  }, {tier: 'ui', description: 'a Shift-drag from above left of those markers to below right of them, 12 px around: the platform\'s rectangle selection'});
