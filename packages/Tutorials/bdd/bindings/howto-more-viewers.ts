/* Names and gestures the how-to clips of the More viewers section need and the shell does not give: the tree map's
   split selectors are bare native <select>s with no name, one per level plus an empty one that adds the next
   level, so a clip reaches them by position ("second tree map split selector"); and the tile viewer's drag of a
   card into another lane, which UsageAnalysis defines for its own features only. */
import {type Page} from '@playwright/test';
import {When, element} from '@datagrok-libraries/bdd';
import {type ElementRef, viewers} from '@datagrok-libraries/bdd/runtime';

element('tree map split selector', {selector: 'select.d4-column-selector-tree-map',
  aliases: ['tree map level selector'],
  description: 'a split column of the tree map, in the strip above the map: the first is the outer level, the last (empty) one adds a level'});

/* The drop zone is a body-level overlay the first move past five pixels creates (`ui.makeDroppable`), and moves
   queued while the page is busy are coalesced: the drag starts, crosses, and moves once more inside the lane while
   the zone is there, so the release never lands before the zone exists. */
export const dragTileToLane = When('user drags the card of row {int} of {widget} to the {string} lane',
  async (page: Page, row: number, target: ElementRef, lane: string) => {
    const from = viewers.centerOf(await viewers.hitArea(page, target, `tile of row ${row}`, true));
    const to = viewers.centerOf(await viewers.hitArea(page, target, `lane content "${lane}"`));
    await page.mouse.move(from.x, from.y);
    await page.mouse.down();
    await page.mouse.move(from.x + 12, from.y + 12);
    await page.mouse.move(to.x, to.y, {steps: 6});
    await page.mouse.move(to.x + 1, to.y + 1);
    await page.mouse.up();
    await viewers.settle(page, target, 1000);
  }, {tier: 'ui', description: 'a card picked up and dropped on another lane: the row takes that lane\'s category'});
