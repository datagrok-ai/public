/* ---
realizes: []
--- */

import {test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

// The server lane of the mixed pc-plot scenario: To Script renders its output through a
// platform function that local mode does not serve, so the balloon never appears there.
// The rest of the scenario is client behaviour and stays in pc-plot-spec.ts.

test('PC Plot tests — menu ribbon and To Script', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  const pcPlotPresent = () => page.evaluate(() => !!grok.shell.tv.viewers.find(v => v.type === 'PC Plot'));

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-pc-plot"]');
    if (icon) (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-PC-Plot"]').waitFor({timeout: 10000});

  await v.installEventWaits(page);

  await softStep('Menu Ribbon and To Script', async () => {
    const viewerExists = await page.evaluate(() => !!grok.shell.tv.viewers.find(v => v.type === 'PC Plot'));
    expect(viewerExists).toBe(true);

    await page.evaluate(async () => {
      const w = window as any;
      const canvas = document.querySelector('[name="viewer-PC-Plot"] canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2}));

      await w.__poll(() => document.querySelectorAll('.d4-menu-item-label').length,
        (n: number) => n > 0, 500);
      document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      await w.__poll(() => document.querySelectorAll('.d4-menu-popup').length,
        (n: number) => n === 0, 800);
    });

    const balloon = await page.evaluate(async () => {
      const w = window as any;
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();

      const reset = async () => {
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
        await w.__poll(() => document.querySelectorAll('.d4-menu-popup').length,
          (n: number) => n === 0, 500);
      };
      const attempt = async () => {
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
        }));
        const opened = await w.__menuLeaf('To Script', 'To JavaScript').catch(() => false);
        if (!opened) { await reset(); return ''; }
        const balloonText = await w.__poll(
          () => (document.querySelector('.d4-balloon') as HTMLElement | null)?.innerText ?? '',
          (t: string) => t.length > 0, 6000, 250);
        if (balloonText.length > 0) return balloonText;
        await reset();
        return '';
      };

      let text = '';
      for (let a = 0; a < 5 && !text; a++) text = await attempt();
      return {present: text.length > 0, text};
    });

    expect(balloon.present).toBe(true);
    expect(balloon.text).toContain('addViewer');

    await page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot');
      if (pc) pc.close();
    });
    await v.pollValue(() => page.evaluate(() => !grok.shell.tv.viewers.find(v => v.type === 'PC Plot')),
      (gone) => gone, 500, 100);
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-pc-plot"]');
      if (icon) (icon as HTMLElement).click();
    });
    const reopened = await v.pollValue(pcPlotPresent, (present) => present, 1000, 100);
    expect(reopened).toBe(true);
  });

  expect(pageErrors).toEqual([]);

  v.finishSpec();
});
