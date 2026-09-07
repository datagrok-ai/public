/* ---
realizes: [filters.cp.cloned-view-sync]
--- */
import {expect} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {activateView, clickCardCheckboxIn, FULL, installViewResolver, trueCountOf} from './cloned-view-shared';

declare const grok: any;

// Scenario 4: a layout carrying one disabled card and one filtering card is saved to the server
// and re-applied over a perturbed view.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Filters — Cloned View Synchronization: saved layout restores per-card state', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, withFilterPanel: true});
  await installViewResolver(page);

  try {
    await softStep('Scenario 4 - Step 3 saved layout restores panel, cards and trueCount', async () => {
      let layoutId: string | null = null;
      const LAYOUT_VIEW = 'LayoutHost';
      try {
        await page.evaluate(async (vn: string) => {
          const view = grok.shell.tv;
          view.dataFrame.name = vn;
          view.name = vn;
          await (window as any).__poll(() => (window as any).__tv(vn), (v2: any) => !!v2, 1500, 25);
        }, LAYOUT_VIEW);
        await v.resetFilters(page);
        await activateView(page, LAYOUT_VIEW);
        await v.applyCategoricalFilter(page, 'RACE', ['Asian']);
        const withBothCards = await v.applyNumericFilter(page, 'AGE', 30, 60);

        await clickCardCheckboxIn(page, LAYOUT_VIEW, 'RACE');
        await expect.poll(async () => page.evaluate((vn: string) =>
          (window as any).__tv(vn).getFiltersGroup().getStates('RACE', 'categorical')[0]?.active, LAYOUT_VIEW),
        {message: 'unticking the RACE card\'s own checkbox did not switch that card off, so the layout '
          + 'about to be saved does not carry the one-card-off state the restore is measured against',
        timeout: 20_000, intervals: [30, 60, 120, 250, 500, 1000]}).toBe(false);

        const saved = await page.evaluate(async (vn: string) => {
          const view = (window as any).__tv(vn);
          const fg = view.getFiltersGroup();
          const raceEl = (window as any).__card(vn, 'RACE') as HTMLElement;
          const ageEl = (window as any).__card(vn, 'AGE') as HTMLElement;
          const layout = view.saveLayout();
          await grok.dapi.layouts.save(layout);
          const serverFound = !!(await (window as any).__findSaved(
            () => grok.dapi.layouts.find(layout.id), 10_000));
          return {
            serverFound,
            savedTrueCount: view.dataFrame.filter.trueCount,
            captions: ((window as any).__cards(vn) as HTMLElement[])
              .map((e) => (e.querySelector('.d4-filter-column-name')?.textContent ?? '').trim()),
            layoutId: layout.id,
            raceActive: fg.getStates('RACE', 'categorical')[0]?.active,
            ageActive: fg.getStates('AGE', 'histogram')[0]?.active,
            raceDisabledClass: raceEl.classList.contains('d4-filter-disabled'),
            ageDisabledClass: ageEl.classList.contains('d4-filter-disabled'),
          };
        }, LAYOUT_VIEW);
        layoutId = saved.layoutId;
        expect(saved.serverFound,
          'the saved layout cannot be fetched back from the server — saveLayout() stamps the id '
          + 'client-side before the round-trip, so the re-apply below would restore nothing')
          .toBe(true);
        expect(saved.raceActive).toBe(false);
        expect(saved.raceDisabledClass).toBe(true);
        expect(saved.ageActive).toBe(true);
        expect(saved.ageDisabledClass).toBe(false);
        expect(saved.savedTrueCount).toBeGreaterThan(withBothCards);
        expect(saved.savedTrueCount).toBeGreaterThan(0);
        expect(saved.savedTrueCount).toBeLessThan(FULL);
        expect(saved.captions.length).toBeGreaterThan(0);

        await v.applyNumericFilter(page, 'AGE', 0, 200);
        await page.evaluate(async (vn: string) => {
          const w = window as any;
          const view = w.__tv(vn);
          const root = view.root as HTMLElement;
          view.getFiltersGroup().close();
          await w.__poll(() => root.querySelectorAll('[name="viewer-Filters"]').length,
            (n: number) => n === 0, 800, 25);
          const added = view.addViewer('Bar chart');
          await w.__poll(() => Array.from(view.viewers).includes(added),
            (there: boolean) => there, 1500, 25);
        }, LAYOUT_VIEW);
        expect(await page.locator('[name="viewer-Filters"]').count()).toBe(0);
        const perturbed = await trueCountOf(page, LAYOUT_VIEW);
        expect(perturbed,
          'the perturbation left df.filter.trueCount at the saved value — the restore assert below could not fail')
          .not.toBe(saved.savedTrueCount);

        const restored = await page.evaluate(async ({id, vn, cap}) => {
          const w = window as any;
          const view = w.__tv(vn);
          const s = await grok.dapi.layouts.find(id);
          if (!s) throw new Error(`the saved layout ${id} is not on the server, so nothing can be re-applied`);
          const applied = new Promise<void>((resolve) => {
            const sub = grok.events.onViewLayoutApplied.subscribe(() => { sub.unsubscribe(); resolve(); });
            setTimeout(resolve, cap);
          });
          view.loadLayout(s);
          await applied;
          const root = view.root as HTMLElement;
          await w.__poll(() => root.querySelector('[name="viewer-Filters"] .d4-filter') != null,
            (there: boolean) => there, 5000, 50);
          const panelOpen = root.querySelector('[name="viewer-Filters"]') != null;
          const cards = w.__cards(vn) as HTMLElement[];
          const captions = cards.map((e) => (e.querySelector('.d4-filter-column-name')?.textContent ?? '').trim());
          const cardEl = (caption: string) =>
            cards.find((x) => x.querySelector('.d4-filter-column-name')?.textContent?.trim() === caption);
          const raceEl = cardEl('RACE');
          const ageEl = cardEl('AGE');
          const fg = view.getFiltersGroup();
          return {
            panelOpen,
            captions,
            trueCount: view.dataFrame.filter.trueCount,
            raceActive: fg.getStates('RACE', 'categorical')[0]?.active,
            ageActive: fg.getStates('AGE', 'histogram')[0]?.active,
            raceDisabledClass: raceEl ? raceEl.classList.contains('d4-filter-disabled') : null,
            ageDisabledClass: ageEl ? ageEl.classList.contains('d4-filter-disabled') : null,
          };
        }, {id: layoutId, vn: LAYOUT_VIEW, cap: 5000});

        expect(restored.panelOpen).toBe(true);
        expect(restored.captions).toEqual(saved.captions);
        expect(restored.raceActive).toBe(false);
        expect(restored.raceDisabledClass).toBe(true);
        expect(restored.ageActive).toBe(true);
        expect(restored.ageDisabledClass,
          'the AGE card came back painted disabled — the re-apply did not restore per-card state')
          .toBe(false);
        expect(restored.trueCount).toBe(saved.savedTrueCount);
      } finally {
        if (layoutId) {
          await page.evaluate(async (id: string) => {
            try { const s = await grok.dapi.layouts.find(id); await grok.dapi.layouts.delete(s); } catch (_) { /* */ }
          }, layoutId);
        }
      }
    });
  } finally {
    await v.cleanupShell(page);
  }

  v.finishSpec();
});
