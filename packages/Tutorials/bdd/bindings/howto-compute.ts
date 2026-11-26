/* What the "Scientific computing" how-to clips (features/howto/11-compute) need besides the library's vocabulary:
   Diff Studio's model view docks the panel with the model's inputs (or, with Edit on, its equations editor) to the
   left of the charts at a fifth of the view, too narrow to read a model's lines on film; a person drags the dock's
   divider to widen it. The fitting form puts each parameter's switch beside its input, and the ribbon items that
   open the fitting and sensitivity views leave the pointer over the new view's menu bar. */
import type {Page} from '@playwright/test';
import {element, When} from '@datagrok-libraries/bdd';
import {type ElementRef, expect, gestures, guideHop, locate} from '@datagrok-libraries/bdd/runtime';

const BAR = 'contains(concat(" ", normalize-space(@class), " "), " splitbar-vertical ")';

/* dock-spawn puts the bar right after the panel docked to the left (or after the column of panels it sits in) */
element('equations editor divider', {
  selector: `xpath=//*[contains(concat(" ", normalize-space(@class), " "), " cm-editor ")]` +
    `/ancestor::*[following-sibling::*[1][${BAR}]][1]/following-sibling::*[1]`,
  aliases: ['model editor divider'],
  description: 'the dock splitter right of Diff Studio\'s equations editor, between it and the charts: dragged to widen the editor',
});

/** The switch the fitting form puts before a parameter's input, lit as a stop of the walk before the press. Written
 * when the library's `user switches on {element}` lit the input instead (fixed in the library on 2026-11-22; this
 * step stays for the fitting clip filmed with it). */
export const switchBeside = When('user switches on the switch beside {element}', async (page: Page, target: ElementRef) => {
  const sw = await gestures.switchOf(page, target);
  if (await gestures.readSwitch(sw) === true)
    return;
  await guideHop(page, sw);
  await sw.click();
  await expect.poll(() => gestures.readSwitch(sw), {message: `the switch beside ${target.phrase} after switching it on`})
    .toBe(true);
}, {tier: 'ui', description: 'the switch a fitting or sensitivity form puts beside a parameter, lit and pressed; idempotent'});

/** A ribbon item that opens a view of its own (Diff Studio's Fit and Sensitivity) leaves the pointer where the new
 * view's menu bar comes up, and the menu under it opens on film: the pointer steps off the ribbon at once. */
export const clickAndStepOff = When('user clicks on {element} and steps off the ribbon', async (page: Page, target: ElementRef) => {
  const box = await (await locate(page, target)).filter({visible: true}).first().boundingBox();
  if (!box)
    throw new Error(`${target.phrase} has no box to click`);
  const x = box.x + box.width / 2;
  const y = box.y + box.height / 2;
  await page.mouse.move(x, y, {steps: 8});
  await page.mouse.click(x, y);
  // down and well to the right, past the dropdown the hover opened (straight down stays inside it, which keeps it
  // open; up onto the view tabs brings a tab's tooltip)
  await page.mouse.move(x + 350, y + 150, {steps: 4});
}, {tier: 'ui', description: 'a click, then the pointer moved down and to the right, out of the menu bar the click brings up ' +
  'and its dropdown (a menu of the bar opens on hover and closes when the pointer leaves it and its dropdown)'});
