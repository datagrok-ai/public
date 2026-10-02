/* The one step only the histogram needs. Its two range inputs are bare `<input>` elements the
   viewer appends to its own root — no `ui-input-root` host, so no `input` kind matches them and
   the library's element phrases cannot reach them. The histogram reports their rectangles as the
   `range min input` / `range max input` hit areas, which is enough to type into them
   (`user enters "30" into the "range min input" area of histogram viewer`); reading the text back
   is what this step adds, and it is the claim that shows a bound the slider clamped is still
   displayed as the user typed it. Everything else in the histogram features is the library's
   `viewers` tier and the platform's data steps (`grok-bdd list-steps`). */
import {expect, Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';

export const rangeInputReads = Then('the range {word} input of histogram viewer should read {string}',
  async (page: Page, which: string, text: string) => {
    if (which !== 'min' && which !== 'max')
      throw new Error(`the histogram has a range min and a range max input, not "${which}"`);
    const input = page.locator(`[name="viewer-Histogram"] .d4-filter-input-${which}`).first();
    await expect(input, `the range ${which} input of histogram viewer`).toHaveValue(text, {timeout: 5000});
  }, {description: 'the text of the range bound as it is shown, which an out-of-range value keeps while the slider clamps'});
