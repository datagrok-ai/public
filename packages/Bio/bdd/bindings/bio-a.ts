/* The step the analyze features needed beyond the library's: the cliff count an activity-cliffs
   plot reports. */
import {Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {expect, pollMs} from '@datagrok-libraries/bdd/runtime';

/** The "N cliffs" button the activity-cliffs analysis puts on its scatter plot
 * (`libraries/ml/src/viewers/activity-cliffs.ts`) counts the pairs it drew as cliffs; the `sali`
 * column and the plot exist whether or not any pair was found. */
export const cliffCount = Then('the activity cliffs plot should report at least {int} cliff(s)', async (page: Page, n: number) => {
  let text = '';
  await expect.poll(async () => {
    text = await page.evaluate(() => Array.from(document.querySelectorAll('.cliffs_grid'))
      .filter((b) => (b as HTMLElement).offsetParent != null).map((b) => b.textContent ?? '').pop() ?? '');
    const m = /^(\d+) cliffs?$/.exec(text.trim());
    return m != null && Number(m[1]) >= n;
  }, {timeout: pollMs(10000), message: `the cliffs button of the newest activity-cliffs plot reads "${text}"`}).toBe(true);
}, {description: 'the newest visible "N cliffs" button on an activity-cliffs scatter plot'});
