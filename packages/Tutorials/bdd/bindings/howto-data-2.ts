/* Steps and data the second round of "Work with data" how-to clips (features/howto/05-data, 16 on) and the Chem
   menu tour (features/howto/07-molecules/10-chem-menu-tour) need and the library's vocabulary lacks. The first
   round's steps are in howto-data.ts. */
import type {Page} from '@playwright/test';
import {dataset, When} from '@datagrok-libraries/bdd';
import {menuNames, openTopMenu, pollMs} from '@datagrok-libraries/bdd/runtime';

// Two 96-well plates of one assay in two files, the same eight columns; each file is one plate (its barcode).
dataset('plate-96-1', {path: 'System:DemoFiles/hts/plate-96-1.csv',
  description: '96 wells of plate 1737398: barcode, row, col, pos, volume, concentration, activity, role'});
dataset('plate-96-2', {path: 'System:DemoFiles/hts/plate-96-2.csv',
  description: '96 wells of plate 5899277, the same columns and layout as plate-96-1'});
// Its image column holds wikimedia URLs, which the grid would fetch: a clip removes it off camera.
dataset('tennis-players', {path: 'System:DemoFiles/tennis_players.csv',
  description: '814 tennis players from DBpedia: name, height, plays, birthDate, careerPrizeMoney, image, country, ' +
    'birthPlace, abstract; the view is "tennis_players"'});
dataset('titanic', {path: 'System:DemoFiles/titanic.csv',
  description: '1310 Titanic passengers; name reads "Surname, Title. Given names" ("Allen, Miss. Elisabeth Walton")'});

type Point = {x: number; y: number};
type Box = Point & {width: number; height: number};

/** Where the page last saw the pointer: a listener the step installs once, since the next move starts from there. */
async function pointerAt(page: Page): Promise<Point | null> {
  return page.evaluate(() => {
    const w = window as any;
    if (!w.__howtoPointer) {
      w.__howtoPointer = {at: null};
      document.addEventListener('mousemove', (e) => w.__howtoPointer.at = {x: e.clientX, y: e.clientY}, true);
    }
    return w.__howtoPointer.at;
  });
}

/** The library's `guide.waypoint`, which its runtime does not export: the last move is kept as a corner of the
 * pointer's way, so the video goes round it as the pointer did rather than straight across the menu. */
async function keepCorner(page: Page): Promise<void> {
  await page.evaluate(() => {
    const log = (window as any).__bddGuidePointer as {type: string; held?: boolean; kept?: boolean}[] | undefined;
    const last = log?.[log.length - 1];
    if (last?.type === 'move' && !last.held)
      last.kept = true;
  }).catch(() => undefined);
}

async function labelBox(page: Page, name: string): Promise<Box | null> {
  const item = page.locator(`[name="${name}"]`).filter({visible: true}).first();
  const label = item.locator(':scope > .d4-menu-item-label').first();
  return (await label.boundingBox().catch(() => null)) ?? item.boundingBox().catch(() => null);
}

/** A tour of a top menu that goes the way a person's hand does. The Dart menu hides a group's flyout when the
 * pointer leaves the group at a steep angle, and opens a group only on a second move inside it: so the pointer
 * enters a flyout along its group's row and then goes down it, leaves a flyout along a row back into the
 * dropdown, and settles on a group with two moves. The menu is opened by the library's walk when the item is not
 * showing. Nothing is picked; a group's flyout is awaited. */
export const restOnTopMenuItem = When('user rests the pointer on {string} in the top menu', async (page: Page, path: string) => {
  const {segments, names} = menuNames(path);
  const name = names[names.length - 1];
  const item = page.locator(`[name="${name}"]`).filter({visible: true}).first();
  let from = await pointerAt(page);
  if (!await item.isVisible()) {
    await openTopMenu(page, segments.slice(0, -1).join(' > '), false);
    from = null;
  }
  await item.waitFor({state: 'visible', timeout: pollMs(5000)});
  const box = (await labelBox(page, name))!;
  const to = {x: box.x + Math.min(40, box.width / 2), y: box.y + box.height / 2};
  const parent = names.length > 2 ? await labelBox(page, names[names.length - 2]) : null;
  // the dropdown or flyout the item is in: its group's child container (role menu)
  const column = await item.evaluate((e) => {
    const r = (e.parentElement?.closest('[role="menu"]') ?? e.parentElement!).getBoundingClientRect();
    return {left: r.left, right: r.right, top: r.top, bottom: r.bottom};
  });
  from ??= parent ? {x: parent.x + Math.min(40, parent.width / 2), y: parent.y + parent.height / 2} : null;
  if (from && (from.x < column.left || from.x > column.right)) {
    if (from.x > column.right) {
      // out of a flyout on the right: up or down it to a row the dropdown has too, then left along that row
      const y = Math.min(Math.max(from.y, column.top + 6), column.bottom - 6);
      if (y !== from.y) {
        await page.mouse.move(from.x, y, {steps: 4});
        await keepCorner(page);
      }
      await page.mouse.move(to.x, y, {steps: 5});
    }
    else
      // into a flyout: along the group's row first
      await page.mouse.move(to.x, from.y, {steps: 5});
    await keepCorner(page);
  }
  await page.mouse.move(to.x, to.y, {steps: 6});
  await page.mouse.move(to.x + 2, to.y);
  await page.mouse.move(to.x, to.y);
  if (await item.evaluate((e) => !!e.querySelector(':scope > .d4-menu-item-container-fixed > [role="menu"] > .d4-menu-item')))
    await page.locator(`[name^="${name}---"]`).filter({visible: true}).first()
      .waitFor({state: 'visible', timeout: pollMs(3000)}).catch(() => undefined);
}, {tier: 'ui', description: 'opens the menu if it is not (the library\'s walk), then moves the pointer the way a ' +
  'hand does (into a flyout along its group\'s row, out of it along a row) onto the item and lets a group open; picks nothing'});
