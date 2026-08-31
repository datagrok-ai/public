import {test, expect} from '@playwright/test';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';

test.use(specTestOptions);

// Placement thresholds are two-sided, so sweeping the viewer through the switch band and
// back may change the slot at most once each way, at different sizes. A one-sided threshold
// oscillates and produces many more. The exact-two-transition arithmetic is pinned on the
// pure resolver in d4/test/legend/; this asserts the same property through the real DOM.

test('Legend hysteresis — a sweep through the band does not oscillate', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv'});
  await v.addLegendViewers(page, {column: 'RACE', viewers: ['Scatter plot']});
  await v.waitForLegendIdle(page, 'Scatter plot');

  const widths: number[] = [];
  for (let w = 600; w >= 400; w -= 20) widths.push(w);
  for (let w = 420; w <= 600; w += 20) widths.push(w);

  const slots: string[] = [];
  for (const w of widths) {
    await v.resizeViewer(page, 'Scatter plot', w, 1000 - w);
    const state = await v.getLegendState(page, 'Scatter plot');
    slots.push(state?.slot ?? '');
  }

  const transitions: number[] = [];
  for (let i = 1; i < slots.length; i++)
    if (slots[i] !== slots[i - 1]) transitions.push(i);

  await softStep('At most two slot changes across the sweep', async () => {
    expect(transitions.length, `slots: ${slots.join(',')}`).toBeLessThanOrEqual(2);
  });

  await softStep('The down-sweep and up-sweep switch at different widths', async () => {
    if (transitions.length === 2)
      expect(widths[transitions[0]]).not.toBe(widths[transitions[1]]);
  });

  await softStep('Cleanup', async () => { await v.cleanupShell(page); });
  v.finishSpec();
});

test('Legend hysteresis — the mini-icon boundary does not flicker', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv'});
  await v.addLegendViewers(page, {column: 'RACE', viewers: ['Scatter plot']});
  await v.waitForLegendIdle(page, 'Scatter plot');

  const modes: string[] = [];
  for (let s = 400; s >= 200; s -= 20) {
    await v.resizeViewer(page, 'Scatter plot', s, s);
    modes.push((await v.getLegendState(page, 'Scatter plot'))?.mode ?? '');
  }
  for (let s = 220; s <= 400; s += 20) {
    await v.resizeViewer(page, 'Scatter plot', s, s);
    modes.push((await v.getLegendState(page, 'Scatter plot'))?.mode ?? '');
  }

  await softStep('At most two mode changes over the whole sweep', async () => {
    let changes = 0;
    for (let i = 1; i < modes.length; i++)
      if (modes[i] !== modes[i - 1]) changes++;

    expect(changes, `modes: ${modes.join(',')}`).toBeLessThanOrEqual(2);
  });

  await softStep('Cleanup', async () => { await v.cleanupShell(page); });
  v.finishSpec();
});
