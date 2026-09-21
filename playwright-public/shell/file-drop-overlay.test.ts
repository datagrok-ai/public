import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

// GROK-20507: the "Incoming!" overlay that greets an external file drag must always be
// dismissible. A drag whose end the page never hears about (dragleave/drop lost while the
// UI was busy) used to leave it over the whole workspace.
const overlay = '.grok-file-drop-overlay';
const gridCanvas = '.grok-table-view .d4-grid canvas';

function fileDragEvent(type: string, target: Element): void {
  const dt = new DataTransfer();
  dt.items.add(new File(['a,b\n1,2\n'], 'drop.csv', {type: 'text/csv'}));
  target.dispatchEvent(new DragEvent(type, {bubbles: true, cancelable: true, dataTransfer: dt}));
}

async function dragFileIn(page: Page): Promise<void> {
  await page.evaluate(`(${fileDragEvent})('dragenter', document.querySelector('${gridCanvas}'))`);
  await expect(page.locator(overlay)).toHaveCount(1);
}

test('File drop overlay: Esc, click, mouse move and drop all dismiss it', async ({page}) => {
  test.setTimeout(180_000);
  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});

  await test.step('Esc while the grid is focused (its own Esc handler stops propagation)', async () => {
    const box = (await page.locator(gridCanvas).first().boundingBox())!;
    await page.mouse.click(box.x + 200, box.y + 60);
    await dragFileIn(page);
    await page.keyboard.press('Escape');
    await expect(page.locator(overlay)).toHaveCount(0);
  });

  await test.step('Click', async () => {
    await dragFileIn(page);
    await page.evaluate((sel) => (document.querySelector(sel) as HTMLElement).click(), overlay);
    await expect(page.locator(overlay)).toHaveCount(0);
  });

  await test.step('Mouse move: no mouse events arrive during a drag, so a move means it ended', async () => {
    await dragFileIn(page);
    await page.mouse.move(600, 500);
    await page.mouse.move(620, 510);
    await expect(page.locator(overlay)).toHaveCount(0);
  });

  await test.step('Drop still reaches the overlay after dismissals', async () => {
    await dragFileIn(page);
    await page.evaluate(`(${fileDragEvent})('drop', document.querySelector('${overlay}'))`);
    await expect(page.locator(overlay)).toHaveCount(0);
  });
});
