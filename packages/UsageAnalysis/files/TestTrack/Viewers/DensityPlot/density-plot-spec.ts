/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';

const isBenignError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Unable to find element in cloned iframe/.test(text) ||
  /Stack trace [A-Za-z]+/.test(text) ||
  /NullError: method not found: '\w+' on null/.test(text);

async function togglePropCheckbox(page: Page, rowName: string) {
  await page.locator(`[name="${rowName}"] input[type="checkbox"]`)
    .waitFor({state: 'attached', timeout: 15000});
  const before = await page.evaluate((rn: string) =>
    (document.querySelector(`[name="${rn}"] input[type="checkbox"]`) as HTMLInputElement)?.checked, rowName);
  await page.evaluate((rn: string) => {
    (document.querySelector(`[name="${rn}"] input[type="checkbox"]`) as HTMLInputElement)?.click();
  }, rowName);
  await page.waitForFunction(({rn, b}) => {
    const cb = document.querySelector(`[name="${rn}"] input[type="checkbox"]`) as HTMLInputElement;
    return !!cb && cb.checked === !b;
  }, {rn: rowName, b: before}, {timeout: 5000});
  await v.waitForViewerRendered(page, 'Density plot', 200);
}

const readProp = (page: Page, prop: string) => page.evaluate((p) => {
  const dp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
  return dp?.props[p];
}, prop);

test('Density plot tests', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error' && !isBenignError(m.text()))
      consoleErrors.push(m.text());
  });
  const errCount = () => pageErrors.length + consoleErrors.length;

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'density-plot', 'Density-plot');

  const propRowPresent = () => page.evaluate(() =>
    !!document.querySelector('[name="prop-show-color-scale"] input[type="checkbox"]'));
  let panelOpen = false;
  for (let i = 0; i < 8 && !panelOpen; i++) {
    panelOpen = await propRowPresent();
    if (panelOpen) break;
    await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Density-plot"]') as HTMLElement;
      const panel = (viewer?.closest('.panel-base') ?? viewer?.parentElement?.parentElement) as HTMLElement | null;
      const gear = panel?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      gear?.click();
    });
    panelOpen = await v.pollValue(propRowPresent, (open) => open, 1200, 100);
  }
  expect(panelOpen, 'Density plot settings panel did not open').toBe(true);

  await softStep('Show/Hide Color Scale — prop round-trip over the error-free floor', async () => {
    const errBefore = errCount();
    await togglePropCheckbox(page, 'prop-show-color-scale');
    const off = await readProp(page, 'showColorScale');
    await togglePropCheckbox(page, 'prop-show-color-scale');
    const on = await readProp(page, 'showColorScale');

    expect(off).toBe(false);
    expect(on).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Axis Visibility — X/Y show props round-trip over the error-free floor', async () => {
    const errBefore = errCount();
    await togglePropCheckbox(page, 'prop-show-x-axis');
    await togglePropCheckbox(page, 'prop-show-y-axis');
    const hidden = await page.evaluate(() => {
      const dp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      return {x: dp.props.showXAxis, y: dp.props.showYAxis};
    });
    await togglePropCheckbox(page, 'prop-show-x-axis');
    await togglePropCheckbox(page, 'prop-show-y-axis');
    const shown = await page.evaluate(() => {
      const dp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      return {x: dp.props.showXAxis, y: dp.props.showYAxis};
    });

    expect(hidden).toEqual({x: false, y: false});
    expect(shown).toEqual({x: true, y: true});
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Show/Hide Selectors and Bin Slider — computed visibility flips both ways', async () => {
    const errBefore = errCount();
    const selVisibility = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Density-plot"]')!;
      const x = root.querySelector('[name="div-column-combobox-x"]') as HTMLElement | null;
      const y = root.querySelector('[name="div-column-combobox-y"]') as HTMLElement | null;
      const bin = root.querySelector('input[type="range"]') as HTMLElement | null;
      return {
        x: x ? getComputedStyle(x).visibility : 'missing',
        y: y ? getComputedStyle(y).visibility : 'missing',
        binPresent: !!bin,
      };
    });
    const before = await selVisibility();
    await togglePropCheckbox(page, 'prop-show-x-selector');
    await togglePropCheckbox(page, 'prop-show-y-selector');
    await togglePropCheckbox(page, 'prop-show-bin-selector');
    const hidden = await selVisibility();
    const binPropOff = await readProp(page, 'showBinSelector');
    await togglePropCheckbox(page, 'prop-show-x-selector');
    await togglePropCheckbox(page, 'prop-show-y-selector');
    await togglePropCheckbox(page, 'prop-show-bin-selector');
    const shown = await selVisibility();
    const binPropOn = await readProp(page, 'showBinSelector');

    expect(before.x).toBe('visible');
    expect(before.y).toBe('visible');
    expect(hidden.x).toBe('hidden');
    expect(hidden.y).toBe('hidden');
    expect(shown.x).toBe('visible');
    expect(shown.y).toBe('visible');

    expect(before.binPresent).toBe(true);
    expect(binPropOff).toBe(false);
    expect(binPropOn).toBe(true);
    expect(shown.binPresent).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Title and Description — title in the panel titlebar, description at the configured position, clearing removes both', async () => {
    const errBefore = errCount();

    await v.setViewerProps(page, 'Density plot', [{set: {
      showTitle: true,
      title: 'Density Distribution',
      description: 'AGE vs HEIGHT density',
      descriptionVisibilityMode: 'Always',
      descriptionPosition: 'Bottom',
    }, wait: 700}]);
    const set = await v.pollValue(() => page.evaluate(() => {
      const dp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      const root = dp.root as HTMLElement;
      const panelBase = root.closest('.panel-base') as HTMLElement | null;
      const titlebar = panelBase
        ? Array.from(panelBase.querySelectorAll('.panel-titlebar-text'))
          .map((e) => (e as HTMLElement).innerText).join(' | ')
        : '';

      const rootRect = root.getBoundingClientRect();
      let descEl: HTMLElement | null = null;
      const walker = document.createTreeWalker(root, NodeFilter.SHOW_ELEMENT);
      while (walker.nextNode()) {
        const el = walker.currentNode as HTMLElement;
        if (el.children.length === 0 && (el.innerText || '').includes('AGE vs HEIGHT density')) {
          descEl = el;
          break;
        }
      }
      const descRect = descEl ? descEl.getBoundingClientRect() : null;
      const descInBottomHalf = descRect
        ? (descRect.top + descRect.height / 2) > (rootRect.top + rootRect.height / 2)
        : false;
      return {
        titlebar,
        descInDom: (root.innerText || '').includes('AGE vs HEIGHT density'),
        descInBottomHalf,
      };
    }), (r) => r.titlebar.includes('Density Distribution') && r.descInDom && r.descInBottomHalf, 700, 100);
    await v.setViewerProps(page, 'Density plot', [{set: {title: '', description: ''}, wait: 700}]);
    const cleared = await v.pollValue(() => page.evaluate(() => {
      const dp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      const root = dp.root as HTMLElement;
      const panelBase = root.closest('.panel-base') as HTMLElement | null;
      const titlebar = panelBase
        ? Array.from(panelBase.querySelectorAll('.panel-titlebar-text'))
          .map((e) => (e as HTMLElement).innerText).join(' | ')
        : '';
      return {
        titlebar,
        descInDom: (root.innerText || '').includes('AGE vs HEIGHT density'),
      };
    }), (r) => !r.titlebar.includes('Density Distribution') && !r.descInDom, 700, 100);

    expect(set.titlebar).toContain('Density Distribution');
    expect(set.descInDom).toBe(true);
    expect(set.descInBottomHalf).toBe(true);
    expect(cleared.titlebar).not.toContain('Density Distribution');
    expect(cleared.descInDom).toBe(false);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Row Source Filtering — filter round-trip and a table switch with no errors', async () => {
    const errBefore = errCount();

    const [filterSet] = await v.setViewerProps(page, 'Density plot',
      [{set: {filter: '${AGE} > 30'}, wait: 400, read: 'filter'}]);
    const [filterCleared] = await v.setViewerProps(page, 'Density plot',
      [{set: {filter: ''}, wait: 400, read: 'filter'}]);
    expect(filterSet).toBe('${AGE} > 30');
    expect(filterCleared).toBe('');

    await page.evaluate(async (path: string) => {
      const df2 = await grok.dapi.files.readCsv(path);
      df2.name = 'spgi-100';
      grok.shell.addTableView(df2);
      await new Promise((resolve) => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
        setTimeout(resolve, 5000);
      });
    }, spgiPath);
    await page.evaluate(() => {
      const demogView = Array.from(grok.shell.tableViews)
        .find((tv: any) => tv.dataFrame?.name === 'demog' || tv.dataFrame?.rowCount === 5850) as any;
      if (demogView) grok.shell.v = demogView;
    });
    await v.pollValue(() => page.evaluate(() => grok.shell.tv?.dataFrame?.rowCount),
      (n) => n === 5850, 600, 50);

    const switched = await page.evaluate(async () => {
      const dp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      const spgi = Array.from(grok.shell.tables).find((t: any) => t.rowCount === 100) as any;
      const demog = Array.from(grok.shell.tables).find((t: any) => t.rowCount === 5850) as any;
      if (!dp || !spgi || !demog) return {ok: false, dp: !!dp, spgi: !!spgi, demog: !!demog};
      const bound = async (t: any) => {
        for (let i = 0; i < 14 && dp.dataFrame?.rowCount !== t.rowCount; i++)
          await new Promise((r) => setTimeout(r, 50));
      };
      dp.props.table = spgi.name;
      await bound(spgi);
      const boundToSpgi = dp.dataFrame?.rowCount;
      const rootAttachedSpgi = document.body.contains(dp.root);
      dp.props.table = demog.name;
      await bound(demog);
      const boundToDemog = dp.dataFrame?.rowCount;
      const rootAttachedDemog = document.body.contains(dp.root);
      return {ok: true, boundToSpgi, boundToDemog, rootAttachedSpgi, rootAttachedDemog};
    });

    expect(switched.ok).toBe(true);
    expect(switched.boundToSpgi).toBe(100);
    expect(switched.rootAttachedSpgi).toBe(true);
    expect(switched.boundToDemog).toBe(5850);
    expect(switched.rootAttachedDemog).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await page.evaluate(() => grok.shell.closeAll());

  v.finishSpec();
});
