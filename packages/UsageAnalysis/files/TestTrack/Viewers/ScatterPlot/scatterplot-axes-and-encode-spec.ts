/* ---
realizes: [scatterplot.cp.axes-and-encode, viewers.scatter-plot]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as sp from './scatterplot-shared';

declare const grok: any;

// Scenarios 1, 3 and 4 on the local lane. Scenario 2 needs the Formula Lines dialog, which the
// PowerPack package contributes, so it lives with the round-trips (Scenario 5) in
// scatterplot-axes-and-encode-server-spec.ts.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const readConfig = (page: Page) => page.evaluate(() => {
  const s = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
  const label = (r: string) => (s.root.querySelector(
    `[name="div-column-combobox-${r}"] .d4-column-selector-column`)?.textContent ?? '').trim();
  return {
    x: s.props.xColumnName, y: s.props.yColumnName, color: s.props.colorColumnName,
    size: s.props.sizeColumnName, markers: s.props.markersColumnName,
    labels: {x: label('x'), y: label('y'), color: label('color'), size: label('size')},
  };
});

async function setAxisBound(page: Page, rowName: string, value: string | null): Promise<void> {
  const cell = page.locator(`[name="${rowName.replace(/^prop-/, 'prop-view-')}"]`);
  await cell.scrollIntoViewIfNeeded();
  await cell.click();
  await v.pollValue(() => page.evaluate(() => document.activeElement instanceof HTMLInputElement), (f) => f, 2000, 25);
  await page.keyboard.press('Control+A');
  if (value === null) await page.keyboard.press('Delete');
  else await page.keyboard.type(value);
  await page.keyboard.press('Enter');
  await sp.propIs(page, sp.rowProp(rowName), value === null ? null : Number(value), 2500);
}

async function toggleInvertX(page: Page, value: boolean): Promise<void> {
  await page.locator('[name="prop-invert-x-axis"] input[type="checkbox"]').click();
  expect(await sp.propIs(page, 'invertXAxis', value, 2500)).toBe(value);
}

test('Scatter Plot — Axes and Encodings', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const errors = sp.trackErrors(page);
  const errCount = errors.count;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await sp.addScatterPlot(page);

  await softStep('Scenario 1 — Set the axes and the encodings through the on-viewer selectors (GROK-18411)', async () => {
    const errBefore = errCount();

    await sp.pickOnViewer(page, 'x', 'AGE');
    await sp.pickOnViewer(page, 'y', 'HEIGHT');

    await sp.pickOnViewer(page, 'color', 'RACE');
    const labelHit = await v.pickColumnViaSelectorTrusted(page, {
      role: 'color', columnName: 'RACE', target: 'column', requirePopup: false,
    });
    expect(labelHit.popupOpened).toBe(true);

    await sp.pickOnViewer(page, 'size', 'WEIGHT');

    await sp.pickPanelColumn(page, 'prop-markers', 'div-column-combobox-markers', 'marker', 'SEX', 'markersColumnName');

    await sp.pickOnViewer(page, 'x', 'WEIGHT');
    await sp.pickOnViewer(page, 'x', 'AGE');

    const cfg = await v.pollValue(() => readConfig(page), (c) => c.labels.x === 'AGE', 2000, 50);
    expect(cfg.x).toBe('AGE');
    expect(cfg.y).toBe('HEIGHT');
    expect(cfg.color).toBe('RACE');
    expect(cfg.size).toBe('WEIGHT');
    expect(cfg.markers).toBe('SEX');
    expect(cfg.labels.x).toBe('AGE');
    expect(cfg.labels.y).toBe('HEIGHT');
    expect(cfg.labels.color).toBe('RACE');
    expect(cfg.labels.size).toBe('WEIGHT');
    expect(errors.all().slice(errBefore)).toEqual([]);
  });

  await softStep('Scenario 3 — Logarithmic and inverted axis with a reversed range window (GROK-13110)', async () => {
    const errBefore = errCount();
    await sp.setChoiceProp(page, 'prop-x-axis-type', 'x-axis', 'logarithmic');
    await toggleInvertX(page, true);

    await setAxisBound(page, 'prop-x-min', '60');
    await setAxisBound(page, 'prop-x-max', '20');
    await v.pollValue(() => page.evaluate(() => {
      const s = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
      const vp = s?.viewport;
      return !!vp && Number.isFinite(vp.width) && vp.width > 0 && Number.isFinite(vp.height) && vp.height > 0;
    }), (ok) => ok, 3000, 50);

    const state = await page.evaluate(() => {
      const s = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
      const vp = s.viewport;
      return {
        axisType: s.props.xAxisType, invert: s.props.invertXAxis,
        xMin: s.props.xMin, xMax: s.props.xMax,
        attached: document.body.contains(s.root),
        vp: vp && {x: vp.x, y: vp.y, width: vp.width, height: vp.height},
      };
    });
    expect(state.axisType).toBe('logarithmic');
    expect(state.invert).toBe(true);
    expect(state.xMin).toBe(60);
    expect(state.xMax).toBe(20);
    expect(state.attached).toBe(true);
    expect(state.vp).not.toBeNull();
    for (const k of ['x', 'y', 'width', 'height'] as const)
      expect(Number.isFinite((state.vp as any)[k])).toBe(true);
    expect(state.vp!.width).toBeGreaterThan(0);
    expect(state.vp!.height).toBeGreaterThan(0);
    const raised = errors.all().slice(errBefore);
    expect(raised.filter((e) => /Wrong range/i.test(e))).toEqual([]);
    expect(raised).toEqual([]);

    await setAxisBound(page, 'prop-x-min', null);
    await setAxisBound(page, 'prop-x-max', null);
    await toggleInvertX(page, false);
    await sp.setChoiceProp(page, 'prop-x-axis-type', 'x-axis', 'linear');
    const back = await page.evaluate(() => {
      const s = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
      return {axisType: s.props.xAxisType, invert: s.props.invertXAxis,
        xMin: s.props.xMin ?? null, xMax: s.props.xMax ?? null};
    });
    expect(back.axisType).toBe('linear');
    expect(back.invert).toBe(false);
    expect(back.xMin).toBeNull();
    expect(back.xMax).toBeNull();
  });

  await softStep('Scenario 4 — Axis type control disabled for a datetime axis (GROK-20395)', async () => {
    const errBefore = errCount();
    await sp.openSettings(page);
    await sp.revealPropEditor(page, '[name="prop-view-x-axis-type"]', 'x-axis');

    await sp.pickOnViewer(page, 'x', 'STARTED');
    expect(await v.pollValue(() => sp.rowOpacity(page, 'prop-x-axis-type'), (o) => o === '0.5', 3000, 50)).toBe('0.5');
    expect(await sp.rowOpacity(page, 'prop-x-map')).toBe('1');

    await sp.pickOnViewer(page, 'x', 'AGE');
    expect(await v.pollValue(() => sp.rowOpacity(page, 'prop-x-axis-type'), (o) => o === '1', 3000, 50)).toBe('1');
    expect(await sp.rowOpacity(page, 'prop-x-map')).toBe('0.5');

    const cfg = await readConfig(page);
    expect(cfg.x).toBe('AGE');
    expect(errors.all().slice(errBefore)).toEqual([]);
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
