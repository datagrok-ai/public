import {test, expect} from '@playwright/test';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';

test.use(specTestOptions);

// The legend is rebuilt inside the viewer's render(), so every user action must cost
// exactly one committed legend frame. These specs assert that arithmetic directly —
// they are the regression net against a timer-driven refresh storm.

test('Legend lifecycle — one frame per action, no double layout', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', withFilterPanel: true});
  await v.addLegendViewers(page, {column: 'RACE', viewers: ['Scatter plot']});

  const settled = await v.waitForLegendIdle(page, 'Scatter plot');
  expect(settled, 'legend never settled').not.toBeNull();

  await softStep('No double layout', async () => {
    expect(settled!.layoutPasses).toBe(settled!.frame);
  });

  await softStep('A filter costs exactly one frame', async () => {
    const before = settled!.frame;
    await page.evaluate(async () => {
      (window as any).grok.shell.tv.dataFrame.rows.match({'SEX': 'F'}).filter();
      await new Promise((r) => setTimeout(r, 1500));
    });
    const after = await v.waitForLegendIdle(page, 'Scatter plot');
    expect(after!.frame - before).toBe(1);
  });

  await softStep('A legend click never resizes or rebuilds the legend', async () => {
    const before = await v.waitForLegendIdle(page, 'Scatter plot');
    const keysBefore = await v.getLegendItemKeys(page, 'Scatter plot');
    await page.locator('[name="viewer-Scatter-plot"] [data-item-key]').first().click();
    const after = await v.waitForLegendIdle(page, 'Scatter plot');

    // Two commits: the click repaints the selection, then the filter it applies re-renders
    // the viewer. The corner ladder may re-run on that repaint (it did in the legacy path
    // too, and it is what keeps the legend off the data), but the content and the size are
    // signature-gated and must not move.
    expect(after!.frame - before!.frame).toBeLessThanOrEqual(2);
    expect(after!.box!.width).toBe(before!.box!.width);
    expect(after!.box!.height).toBe(before!.box!.height);
    expect(after!.mode).toBe(before!.mode);
    expect(await v.getLegendItemKeys(page, 'Scatter plot')).toEqual(keysBefore);
  });

  await softStep('Cleanup', async () => { await v.resetFilters(page); await v.cleanupShell(page); });
  v.finishSpec();
});

test('Legend lifecycle — the legend does not jump after settling', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv'});
  await v.addLegendViewers(page, {column: 'RACE', viewers: ['Scatter plot']});
  await v.waitForLegendIdle(page, 'Scatter plot');

  // The direct "must not jump" assertion: once frame stops advancing, nothing may move.
  const first = await v.getLegendState(page, 'Scatter plot');
  await page.waitForTimeout(1000);
  const second = await v.getLegendState(page, 'Scatter plot');

  await softStep('Frame did not advance while idle', async () => {
    expect(second!.frame).toBe(first!.frame);
  });

  await softStep('Geometry is unchanged while idle', async () => {
    expect(second!.box).toEqual(first!.box);
  });

  await softStep('Cleanup', async () => { await v.cleanupShell(page); });
  v.finishSpec();
});
