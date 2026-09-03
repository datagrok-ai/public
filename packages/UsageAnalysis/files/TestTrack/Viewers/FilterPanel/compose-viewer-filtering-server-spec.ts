/* ---
realizes: [filters.cp.compose-with-viewer-filtering]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {expectHeaderCounter, trueCount} from '../../helpers/filter-panel';
import {addViewer, FULL, PANEL_CATEGORY, raceSelectedCategories, seedPanelCriterion, viewerCanvasRect,
  zoomScatterPlot} from './compose-viewer-shared';

declare const grok: any;

// Scenario 2 Step 5: the layout saved through the real View > Layout > Save to Gallery leaf while a
// panel criterion and a scatter-plot zoom both filter, and re-applied over the API (GROK-18281).
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const runTag = `compose-${Date.now()}-${Math.floor(Math.random() * 1e9)}`;

async function closeFilterPanel(page: Page): Promise<void> {
  await page.locator('[name="viewer-Filters"]').first().hover();
  let clicked = false;
  for (const icon of ['icon-times', 'Close']) {
    try {
      await v.clickViewerTitlebarIcon(page, 'Filters', icon);
      clicked = true;
      break;
    } catch (_) { /* try the other title-bar close control */ }
  }
  expect(clicked, 'the Filters panel exposes no title-bar close control').toBe(true);
  await expect.poll(async () => page.locator('[name="viewer-Filters"]').count(),
    {timeout: 10_000, intervals: [300, 600, 1200]}).toBe(0);
}

async function myLayoutsCarrying(page: Page, marker: string): Promise<string[]> {
  return page.evaluate(async (m: string) => {
    const me = String(grok.shell.user.id);
    const ls = (await grok.dapi.layouts.getApplicable(grok.shell.tv.dataFrame)) ?? [];
    return ls.filter((l: any) => !l.author || !l.author.id || String(l.author.id) === me)
      .filter((l: any) => {
        const name = String(l.friendlyName ?? l.name ?? '');
        return name === m || name.startsWith(`${m} (`);
      })
      .map((l: any) => String(l.id));
  }, marker);
}

test('Filters — layout round-trip with panel and viewer filtering combined', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, withFilterPanel: true});

  const total = await trueCount(page);
  expect(total).toBe(FULL);
  const truePanel = await seedPanelCriterion(page, PANEL_CATEGORY);
  expect(truePanel).toBeLessThan(FULL);
  expect(truePanel).toBeGreaterThan(0);
  await expectHeaderCounter(page, '1',
    'the seeded RACE criterion is the only thing filtering, so the header counter must read 1');

  try {
    await softStep('Scenario 2 Step 5 Layout round-trip with combined filtering → no error, count restored', async () => {
      let savedLayout: string | null = null;
      try {
        await addViewer(page, 'Scatter plot');
        const rect = await viewerCanvasRect(page, 'Scatter plot');
        expect(rect).not.toBeNull();
        await zoomScatterPlot(page, rect!);
        await expect.poll(async () => trueCount(page),
          {timeout: 15_000, intervals: [400, 800, 1500]}).toBeLessThan(truePanel);
        const zoomed = await trueCount(page);
        const beforeSave = await page.evaluate(() => {
          const tv = grok.shell.tv;
          const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
          return {
            filtered: tv.dataFrame.filter.trueCount,
            zoomAndFilter: sp?.props?.zoomAndFilter ?? null,
            xMin: sp?.props?.xMin ?? null, xMax: sp?.props?.xMax ?? null,
            yMin: sp?.props?.yMin ?? null, yMax: sp?.props?.yMax ?? null,
            spFilter: sp ? String(sp.props.filter ?? '') : null,
          };
        });

        const pageErrors: string[] = [];
        const onPageError = (e: Error) => pageErrors.push(e.message);
        page.on('pageerror', onPageError);
        await page.evaluate(() => {
          const w = window as any;
          w.__errBalloons = [];
          const seen = new WeakSet<Element>();
          const record = (el: Element) => {
            if (seen.has(el)) return;
            seen.add(el);
            w.__errBalloons.push((el.textContent ?? '').trim());
          };
          const scan = (n: Node) => {
            if (!(n instanceof Element)) return;
            if (n.matches('.d4-balloon.error')) record(n);
            for (const el of Array.from(n.querySelectorAll('.d4-balloon.error'))) record(el);
          };
          w.__errBalloonObs = new MutationObserver((records: MutationRecord[]) => {
            for (const r of records) {
              if (r.type === 'attributes') {
                const t = r.target as Element;
                if (t.matches('.d4-balloon.error')) record(t);
              } else
                for (const n of Array.from(r.addedNodes)) scan(n);
            }
          });
          w.__errBalloonObs.observe(document.body,
            {childList: true, subtree: true, attributes: true, attributeFilter: ['class']});
        });

        // the table name is what the gallery save names the layout after, so a run-unique marker
        // identifies ours without listing the layouts that existed before
        const stamped = await page.evaluate((name: string) => {
          grok.shell.tv.dataFrame.name = name;
          return String(grok.shell.tv.dataFrame.name);
        }, runTag);
        expect(stamped).toBe(runTag);
        expect(await v.driveTopMenuLeaf(page, ['View', 'Layout', 'Save to Gallery'])).toBe(true);
        let fresh: string[] = [];
        await expect.poll(async () => {
          fresh = await myLayoutsCarrying(page, runTag);
          return fresh.length;
        }, {timeout: 25_000, intervals: [500, 1000, 2000, 3000]}).toBeGreaterThanOrEqual(1);
        expect(fresh.length, `expected exactly 1 new layout, got ${fresh.length}`).toBe(1);
        savedLayout = fresh[0];

        await closeFilterPanel(page);

        await page.evaluate(async (layoutId: string) => {
          const saved = await grok.dapi.layouts.find(layoutId);
          const applied = new Promise<void>((resolve) => {
            const sub = grok.events.onViewLayoutApplied.subscribe(() => { sub.unsubscribe(); resolve(); });
            setTimeout(resolve, 8000);
          });
          grok.shell.tv.loadLayout(saved);
          await applied;
        }, savedLayout);

        await expect.poll(async () => page.locator('[name="viewer-Filters"] .d4-filter').count(),
          {timeout: 20_000, intervals: [500, 1000, 2000, 3000]}).toBeGreaterThanOrEqual(1);
        await expect.poll(async () => raceSelectedCategories(page),
          {timeout: 15_000, intervals: [500, 1000, 2000]}).toEqual([PANEL_CATEGORY]);
        const afterLayout = await page.evaluate(() => {
          const tv = grok.shell.tv;
          const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
          return {
            viewers: Array.from(tv.viewers).map((x: any) => x.type),
            filtered: tv.dataFrame.filter.trueCount,
            zoomAndFilter: sp?.props?.zoomAndFilter ?? null,
            xMin: sp?.props?.xMin ?? null, xMax: sp?.props?.xMax ?? null,
            yMin: sp?.props?.yMin ?? null, yMax: sp?.props?.yMax ?? null,
            spFilter: sp ? String(sp.props.filter ?? '') : null,
          };
        });
        expect(afterLayout.viewers,
          'the re-applied layout did not bring the Scatter Plot back, so the row count settling at ' +
          'the panel-only value below says nothing about a dropped zoom — it would be the count of a ' +
          `view with no scatter plot in it at all; viewers: [${afterLayout.viewers.join(', ')}]`)
          .toContain('Scatter plot');
        expect(afterLayout.zoomAndFilter,
          'the restored Scatter Plot is no longer routing its zoom to the dataframe filter, so the ' +
          'panel-only row count below would be explained by the viewer being disarmed rather than by ' +
          'the zoom itself being dropped from the layout').toBe('filter by zoom');
        expect(zoomed, 'the pre-save zoom was not narrower than the panel criterion alone, so the ' +
          'round-trip bound below could not tell a restored zoom from a lost one').toBeLessThan(truePanel);
        await expect.poll(async () => trueCount(page), {
          timeout: 30_000,
          intervals: [500, 1000, 2000, 3000],
          message: 'the re-applied layout did not settle at the panel-only row count: ' +
            `before=${JSON.stringify(beforeSave)} after=${JSON.stringify(afterLayout)}`,
        }).toBe(truePanel);
        const restored = await trueCount(page);
        expect(restored,
          'the layout round-trip no longer settles at the panel-only row count. Below it means the ' +
          'scatter-plot zoom is now restored with the layout (scope_reductions SR-06 is obsolete — ' +
          'tighten this to toBeLessThan(truePanel)); above it means the panel criterion itself was ' +
          `lost. before=${JSON.stringify(beforeSave)} after=${JSON.stringify(afterLayout)}`)
          .toBe(truePanel);

        const errBalloons = await page.evaluate(async () => {
          const w = window as any;
          let last = (w.__errBalloons ?? []).length;
          let since = Date.now();
          while (Date.now() - since < 2500) {
            await new Promise((r) => setTimeout(r, 200));
            const n = (w.__errBalloons ?? []).length;
            if (n !== last) { last = n; since = Date.now(); }
          }
          w.__errBalloonObs?.disconnect();
          return (w.__errBalloons ?? []) as string[];
        });
        page.off('pageerror', onPageError);
        expect(errBalloons,
          `error balloons appeared during the layout round-trip — GROK-18281: ${errBalloons.join(' | ')}`)
          .toEqual([]);
        expect(pageErrors, `layout re-apply raised page errors — GROK-18281: ${pageErrors.join('; ')}`).toEqual([]);
      } finally {
        await page.evaluate(() => { try { (window as any).__errBalloonObs?.disconnect(); } catch (_) {} });
        if (savedLayout) {
          await page.evaluate(async (layoutId: string) => {
            try { const s = await grok.dapi.layouts.find(layoutId); await grok.dapi.layouts.delete(s); } catch (_) {}
          }, savedLayout);
        }
        await page.evaluate(() => { try { grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot')?.close(); } catch (_) {} });
      }
    });
  } finally {
    await page.evaluate(() => {
      try {
        grok.shell.tv?.dataFrame?.resetFilter();
        grok.shell.tv?.dataFrame?.selection?.setAll(false);
      } catch (_) {}
    });
    await v.cleanupShell(page);
  }

  v.finishSpec();
});
