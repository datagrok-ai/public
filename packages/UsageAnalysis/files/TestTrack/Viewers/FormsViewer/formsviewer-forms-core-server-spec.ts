/* ---
realizes: [formsviewer.int.pinned-rows-persist-by-value]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as projects from '../../helpers/projects';
import {
  HOST, ORDINARY, PINNED, PINNED_PANE,
  cardFieldValue, drawnLabelNames, fieldValuesByPosition, waitForOrderStable, withConsoleErrorCount,
  sortIndicatorLabels, cardContextMenu,
} from '../../helpers/forms';

declare const grok: any;

// The server lane of the mixed forms-core scenario: Steps 7a-7c, whose subject is a layout and a
// project surviving a round-trip through the server. The ladder that proves the viewer itself
// (fields, sort, pinning) is formsviewer-forms-core-spec.ts on the local lane; the configured
// state is set up here directly rather than re-driven step by step.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

interface FormsState { labels: string[]; fields: string[]; pinnedValues: string[]; indicator: string[]; }

async function formsState(page: Page): Promise<FormsState> {
  const [labels, indicator, props] = await Promise.all([
    drawnLabelNames(page),
    sortIndicatorLabels(page),
    page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return {
        fields: Array.from(vw.props.fieldsColumnNames as string[]),
        pinnedValues: Array.from(vw.props.pinnedRowValues as string[]),
      };
    }),
  ]);
  return {labels, indicator, ...props};
}

async function expectFormsState(page: Page, pre: FormsState, timeout: number): Promise<void> {
  await expect.poll(() => drawnLabelNames(page), {timeout}).toEqual(pre.labels);
  await expect.poll(() => page.evaluate(() => {
    const vw = grok.shell.tv?.viewers?.find((x: any) => x.type === 'FormsViewer');
    return vw ? Array.from(vw.props.fieldsColumnNames as string[]) : null;
  }), {timeout}).toEqual(pre.fields);
  await expect.poll(() => sortIndicatorLabels(page), {timeout}).toEqual(pre.indicator);
  await expect.poll(() => page.evaluate((sel) => Array.from(document.querySelectorAll(sel))
    .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value), PINNED),
  {timeout}).toEqual(pre.pinnedValues);
}

test('Forms viewer — layout and project persistence', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Setup — Forms with a reversed field subset, sorted by AGE, one row pinned by value', async () => {
    await v.addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator('.d4-multi-form').first().waitFor({timeout: 30_000});

    const chosenFields = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      const all = Array.from(vw.props.fieldsColumnNames as string[]);
      const must = ['USUBJID', 'AGE'];
      const rest = all.filter((n) => !must.includes(n));
      const picked = [...must, ...rest.slice(0, Math.max(1, rest.length - 2))].reverse();
      vw.setOptions({fieldsColumnNames: picked, sortByColumnName: 'AGE', showMouseOverRow: false});
      const df = grok.shell.t;
      df.mouseOverRowIdx = -1;
      df.currentRowIdx = 0;
      df.selection.setAll(false);
      df.selection.set(5, true); df.selection.set(10, true); df.selection.set(20, true);
      return picked;
    });
    await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(chosenFields);
    await expect.poll(() => sortIndicatorLabels(page), {timeout: 20_000}).toEqual(['div-AGE']);
    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(4);
    await waitForOrderStable(page);

    const pos = (await fieldValuesByPosition(page, 'USUBJID')).findIndex((val, i) => i >= 1 && val !== null);
    expect(pos).toBeGreaterThanOrEqual(0);
    const anchor = await cardFieldValue(page, pos, 'USUBJID');
    await cardContextMenu(page, ORDINARY, pos, 'div-Pin-Row', 'USUBJID');
    await expect.poll(async () =>
      page.evaluate((sel) => getComputedStyle(document.querySelector(sel) as HTMLElement).display, PINNED_PANE),
    {timeout: 15_000}).not.toBe('none');
    await expect.poll(() => page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer').props.pinnedRowValues), {timeout: 15_000})
      .toEqual([anchor]);
  });

  await softStep('Step 7a / Step 7b — Re-applying the saved layout over a deliberately corrupted view restores the field set, sort-label identity and pinned row by value, and drops a foreign viewer not in the layout', async () => {
    const pre = await formsState(page);
    expect(pre.labels.length).toBeGreaterThan(0);
    expect(pre.pinnedValues.length).toBeGreaterThan(0);
    expect(pre.indicator.length).toBeGreaterThan(0);

    // awaiting layouts.save IS the completion signal, the same call the Save to Gallery menu makes
    const layoutId: string = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      layout.name = 'zz-formsviewer-core-' + Date.now();
      return String((await grok.dapi.layouts.save(layout)).id);
    });

    try {
      await page.evaluate(() => {
        const tv = grok.shell.tv;
        tv.viewers.find((x: any) => x.type === 'FormsViewer')?.close();
        tv.addViewer('Histogram');
      });
      await expect.poll(() => page.locator('[name="viewer-Histogram"]').count(), {timeout: 20_000}).toBe(1);
      await expect.poll(() => page.locator(HOST).count(), {timeout: 20_000}).toBe(0);

      const applyErrTexts: string[] = [];
      const applyErrCount = await withConsoleErrorCount(page, async () => {
        await page.evaluate(async (id) => {
          grok.shell.tv.loadLayout(await grok.dapi.layouts.find(id));
        }, layoutId);
        await page.locator(HOST).first().waitFor({timeout: 30_000});
      }, undefined, applyErrTexts);
      expect(applyErrCount, `layout-apply console errors: ${JSON.stringify(applyErrTexts)}`).toBe(0);

      await expect.poll(() => page.locator('[name="viewer-Histogram"]').count(), {timeout: 20_000}).toBe(0);
      expect(await page.evaluate(() => grok.shell.tv.viewers
        .filter((x: any) => x.type === 'Histogram').length)).toBe(0);
      await expectFormsState(page, pre, 20_000);
    } finally {
      await page.evaluate(async (id) => {
        const saved = await grok.dapi.layouts.find(id);
        if (saved) await grok.dapi.layouts.delete(saved);
      }, layoutId);
    }
  });

  await softStep('Step 7c — A project round-trip preserves the field set and the pinned row across a session boundary', async () => {
    const pre = await formsState(page);
    expect(pre.labels.length).toBeGreaterThan(0);
    expect(pre.fields.length).toBeGreaterThan(0);
    expect(pre.pinnedValues.length).toBeGreaterThan(0);
    expect(pre.indicator.length).toBeGreaterThan(0);

    let projectId: string | null = null;
    try {
      projectId = (await projects.saveProjectViaApi(page, `zz-formsviewer-core-${Date.now()}`)).projectId;

      await page.evaluate(async (id) => {
        grok.shell.closeAll();
        for (let i = 0; i < 60 && Array.from(grok.shell.tableViews).length > 0; i++)
          await new Promise((r) => setTimeout(r, 50));
        await (await grok.dapi.projects.find(id)).open();
      }, projectId);
      await page.locator(HOST).first().waitFor({timeout: 30_000});

      await expectFormsState(page, pre, 30_000);
    } finally {
      if (projectId)
        await projects.deleteProjectWithCleanup(page, {projectId});
    }
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
