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

/** Toggle a property-grid boolean checkbox row. Rows in collapsed categories are
 * visibility:hidden but present — a JS-driven .click() works on them, so waits are
 * for the ATTACHED state (not visible) and the flip is confirmed via the checked
 * property, which reads regardless of visibility. */
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
  await page.waitForTimeout(200);
}

/** Read a density-plot prop off the live viewer handle. */
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

  // Open the settings panel and confirm its rows rendered. The gear toggles the
  // panel, so click it only while the panel is closed (detected by the absence
  // of a known prop row) and retry until a checkbox row is present.
  let panelOpen = false;
  for (let i = 0; i < 8 && !panelOpen; i++) {
    panelOpen = await page.evaluate(() =>
      !!document.querySelector('[name="prop-show-color-scale"] input[type="checkbox"]'));
    if (panelOpen) break;
    await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Density-plot"]') as HTMLElement;
      const panel = (viewer?.closest('.panel-base') ?? viewer?.parentElement?.parentElement) as HTMLElement | null;
      const gear = panel?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      gear?.click();
    });
    await page.waitForTimeout(1200);
  }
  expect(panelOpen, 'Density plot settings panel did not open').toBe(true);

  await softStep('Show/Hide Color Scale — prop round-trip over the error-free floor', async () => {
    const errBefore = errCount();
    await togglePropCheckbox(page, 'prop-show-color-scale');
    const off = await readProp(page, 'showColorScale');
    await togglePropCheckbox(page, 'prop-show-color-scale');
    const on = await readProp(page, 'showColorScale');
    // The color scale is canvas-drawn (no DOM signal), so the assertable signal
    // is the prop round-trip plus the error-free repaint.
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
    // Axis lines/labels are canvas-drawn — prop round-trip plus no-error floor.
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
    // Hiding a selector only flips computed visibility (the element stays in DOM
    // with a nonzero rect) — presence/offsetParent checks are false-green.
    expect(before.x).toBe('visible');
    expect(before.y).toBe('visible');
    expect(hidden.x).toBe('hidden');
    expect(hidden.y).toBe('hidden');
    expect(shown.x).toBe('visible');
    expect(shown.y).toBe('visible');
    // The bin slider stays in DOM through the toggle; its signal is the prop.
    expect(before.binPresent).toBe(true);
    expect(binPropOff).toBe(false);
    expect(binPropOn).toBe(true);
    expect(shown.binPresent).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Title and Description — title in the panel titlebar, description at the configured position, clearing removes both', async () => {
    const errBefore = errCount();
    // The property-grid text editors do not fire the Dart change listener
    // reliably, so title/description are driven via the props; the title lands
    // in the enclosing panel-base titlebar and the description in a DOM overlay
    // inside the viewer root — both DOM-readable.
    await page.evaluate(() => {
      const dp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      dp.props.showTitle = true;
      dp.props.title = 'Density Distribution';
      dp.props.description = 'AGE vs HEIGHT density';
      dp.props.descriptionVisibilityMode = 'Always';
      dp.props.descriptionPosition = 'Bottom';
    });
    await page.waitForTimeout(700);
    const set = await page.evaluate(() => {
      const dp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      const root = dp.root as HTMLElement;
      const panelBase = root.closest('.panel-base') as HTMLElement | null;
      const titlebar = panelBase
        ? Array.from(panelBase.querySelectorAll('.panel-titlebar-text'))
          .map((e) => (e as HTMLElement).innerText).join(' | ')
        : '';
      // Locate the leaf element carrying the description text and compare its
      // vertical midpoint to the viewer rect's midpoint (Bottom => lower half).
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
    });
    await page.evaluate(() => {
      const dp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      dp.props.title = '';
      dp.props.description = '';
    });
    await page.waitForTimeout(700);
    const cleared = await page.evaluate(() => {
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
    });
    // The Title renders as the enclosing panel-base titlebar text; clearing the
    // Title field removes it. The Description renders inside the viewer element
    // at the configured position (Bottom => lower half of the viewer rect);
    // clearing the Description removes it.
    expect(set.titlebar).toContain('Density Distribution');
    expect(set.descInDom).toBe(true);
    expect(set.descInBottomHalf).toBe(true);
    expect(cleared.titlebar).not.toContain('Density Distribution');
    expect(cleared.descInDom).toBe(false);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Row Source Filtering — filter round-trip and a table switch with no errors', async () => {
    const errBefore = errCount();
    // Filter formula: the render effect is a canvas-only repaint, so the filter
    // itself is a prop round-trip over the no-error floor.
    await page.evaluate(() => {
      const dp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      dp.props.filter = '${AGE} > 30';
    });
    await page.waitForTimeout(400);
    const filterSet = await readProp(page, 'filter');
    await page.evaluate(() => {
      const dp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      dp.props.filter = '';
    });
    await page.waitForTimeout(400);
    const filterCleared = await readProp(page, 'filter');
    expect(filterSet).toBe('${AGE} > 30');
    expect(filterCleared).toBe('');

    // Open spgi-100 as a second table view, then return to the demog view.
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
    await page.waitForTimeout(600);

    // Rebind the viewer to spgi-100 through the Table prop and back to demog.
    const switched = await page.evaluate(async () => {
      const dp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      const spgi = Array.from(grok.shell.tables).find((t: any) => t.rowCount === 100) as any;
      const demog = Array.from(grok.shell.tables).find((t: any) => t.rowCount === 5850) as any;
      if (!dp || !spgi || !demog) return {ok: false, dp: !!dp, spgi: !!spgi, demog: !!demog};
      dp.props.table = spgi.name;
      await new Promise((r) => setTimeout(r, 700));
      const boundToSpgi = dp.dataFrame?.rowCount;
      const rootAttachedSpgi = document.body.contains(dp.root);
      dp.props.table = demog.name;
      await new Promise((r) => setTimeout(r, 700));
      const boundToDemog = dp.dataFrame?.rowCount;
      const rootAttachedDemog = document.body.contains(dp.root);
      return {ok: true, boundToSpgi, boundToDemog, rootAttachedSpgi, rootAttachedDemog};
    });
    // The table switch re-binds the viewer to the new table's rows and keeps the
    // viewer alive (root attached) with no new errors, then restores the demog binding.
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
