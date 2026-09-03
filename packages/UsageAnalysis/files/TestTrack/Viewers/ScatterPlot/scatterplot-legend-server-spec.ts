/* ---
realizes: [scatterplot.cp.legend, viewers.scatter-plot]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';
import * as sp from './scatterplot-shared';

declare const grok: any;

// The server lane of the legend scenario: the layout and project round-trips at the peak
// configuration (Color RACE, Markers SEX). The legend ladder itself is scatterplot-legend-spec.ts
// on the local lane; the peak state is set up here directly through the API.
test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const RACE_CATEGORIES = ['Asian', 'Black', 'Caucasian', 'Other'];
const SEX_CATEGORIES = ['F', 'M'];

let inProjectSaveWindow = false;
const isBenignError = (text: string) => sp.isAmbientError(text) || (inProjectSaveWindow && (
  /Unable to find element in cloned iframe/.test(text) || /Stack trace [A-Za-z]+/.test(text) ||
  /NullError: method not found: '\w+' on null/.test(text)));

interface Legend { present: boolean; coloring: number; extra: number; colorLabels: string[]; glyphLabels: string[]; }

const readLegend = (page: Page): Promise<Legend> => page.evaluate(() => {
  const root = grok.shell.tv?.viewers?.find((x: any) => x.type === 'Scatter plot')?.root as HTMLElement | undefined;
  const host = root?.querySelector('[name="legend"]') as HTMLElement | null;
  const present = !!host && getComputedStyle(host).display !== 'none';
  const items = present ? [...root!.querySelectorAll('[name="legend"] .d4-legend-item')] as HTMLElement[] : [];
  const txt = (i: Element) => (i.querySelector('.d4-legend-value')?.textContent ?? '').trim();
  const colorItems = items.filter((i) => !i.classList.contains('d4-legend-item-extra'));
  return {
    present,
    coloring: colorItems.length,
    extra: items.length - colorItems.length,
    colorLabels: colorItems.map(txt),
    glyphLabels: items.filter((i) => i.querySelector('i[name="legend-icon-color-picker"]')).map(txt),
  };
});

const viewerColumns = (page: Page) => page.evaluate(() => {
  const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return {color: s?.props.colorColumnName, markers: s?.props.markersColumnName};
});

const peakLegend = (l: Legend) =>
  l.present && l.coloring === RACE_CATEGORIES.length && l.extra === SEX_CATEGORIES.length;

test('Scatter Plot — Legend Persistence', async ({page}: {page: Page}) => {
  test.setTimeout(300_000);

  const errors = sp.trackErrors(page, isBenignError);

  await openDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});

  await softStep('Setup — a scatter plot at the peak configuration, Color RACE and Markers SEX', async () => {
    await page.evaluate(() => grok.shell.tv.addViewer('Scatter plot',
      {xColumnName: 'WEIGHT', yColumnName: 'HEIGHT', colorColumnName: 'RACE', markersColumnName: 'SEX'}));
    await page.locator('[name="viewer-Scatter-plot"]').waitFor({timeout: 15_000});
    await v.installEventWaits(page);
    const legend = await v.pollValue(() => readLegend(page), peakLegend, 10_000, 100);
    expect(legend.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    expect(legend.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    expect(await viewerColumns(page)).toEqual({color: 'RACE', markers: 'SEX'});
  });

  await softStep('Layout and project persistence at peak configuration — layout round-trip', async () => {
    const errBefore = errors.count();
    const layoutId: string = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      layout.name = 'zz-scatterplot-legend-' + Date.now();
      return String((await grok.dapi.layouts.save(layout)).id);
    });
    try {
      await page.evaluate(() => grok.shell.tv.addViewer('Histogram'));
      await expect.poll(() => page.locator('[name="viewer-Histogram"]').count(), {timeout: 10_000}).toBe(1);
      const [clearedColor] = await v.setViewerProps(page, sp.SP_TYPE,
        [{set: {colorColumnName: ''}, wait: 1000, read: 'colorColumnName'}]);
      expect(clearedColor).toBe('');

      await page.evaluate(async (id: string) => {
        grok.shell.tv.loadLayout(await grok.dapi.layouts.find(id));
      }, layoutId);

      const result = await v.pollValue(() => page.evaluate(() => {
        const tv = grok.shell.tv;
        const restored = tv?.viewers.find((x: any) => x.type === 'Scatter plot') as any;
        return {
          hasScatter: !!restored,
          hasHistogram: !!tv?.viewers.some((x: any) => x.type === 'Histogram'),
          color: restored?.props.colorColumnName, markers: restored?.props.markersColumnName,
        };
      }), (r) => r.hasScatter && !r.hasHistogram && r.color === 'RACE', 6000, 100);
      expect(result).toEqual({hasScatter: true, hasHistogram: false, color: 'RACE', markers: 'SEX'});

      const legend = await v.pollValue(() => readLegend(page), peakLegend, 5000, 100);
      expect(legend.present).toBe(true);
      expect(legend.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
      expect(legend.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
      expect(errors.count()).toBe(errBefore);
    } finally {
      await page.evaluate(async (id: string) => {
        const saved = await grok.dapi.layouts.find(id);
        if (saved) await grok.dapi.layouts.delete(saved);
      }, layoutId);
    }
  });

  await softStep('Layout and project persistence at peak configuration — project save / Close All / reopen',
    async () => {
      let projectId: string | null = null;
      inProjectSaveWindow = true;
      try {
        projectId = (await saveProjectViaApi(page, 'zz-scatterplot-legend-probe-' + Date.now())).projectId;
        expect(projectId).toBeTruthy();

        await v.closeAllAndWait(page);
        await page.evaluate(async (id: string) => (await grok.dapi.projects.find(id)).open(), projectId);
        await page.locator('[name="viewer-Scatter-plot"]').waitFor({timeout: 30_000});
        await v.installEventWaits(page);

        const result = await v.pollValue(() => page.evaluate(() => {
          let s: any = null;
          for (const view of grok.shell.tableViews)
            for (const vw of view.viewers)
              if (vw.type === 'Scatter plot') s = vw;
          const items = s ? [...s.root.querySelectorAll('[name="legend"] .d4-legend-item')] : [];
          const txt = (i: Element) => (i.querySelector('.d4-legend-value')?.textContent ?? '').trim();
          return {
            found: !!s,
            color: s?.props.colorColumnName, markers: s?.props.markersColumnName,
            legendPresent: !!s?.root.querySelector('[name="legend"]'),
            colorLabels: items.filter((i) => !i.classList.contains('d4-legend-item-extra')).map(txt),
            glyphLabels: items.filter((i) => i.querySelector('i[name="legend-icon-color-picker"]')).map(txt),
          };
        }), (r) => r.found && r.colorLabels.length === RACE_CATEGORIES.length &&
          r.glyphLabels.length === SEX_CATEGORIES.length, 20_000, 200);

        expect(result.found).toBe(true);
        expect(result.color).toBe('RACE');
        expect(result.markers).toBe('SEX');
        expect(result.legendPresent).toBe(true);
        expect(result.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
        expect(result.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
      } finally {
        inProjectSaveWindow = false;
        if (projectId) await deleteProjectWithCleanup(page, {projectId});
      }
    });

  await v.cleanupShell(page);
  v.finishSpec();
});
