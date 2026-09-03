/* ---
realizes: [filters.cp.hierarchical-and-combined-boolean]
--- */
import {expect} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';
import {addHierarchicalCard, applyBoolState, applyHierarchyState, closeFilterPanelInPage, hierCaption,
  hierNode, openDemogWithSexBool, ROW_COUNT, trueCountOf} from './hierarchical-shared';

declare const grok: any;
declare const window: any;

// Steps 10, 11, 17 and 18: the hierarchical criterion and the combined boolean card through a
// layout round-trip, and the combined boolean card through a project round-trip (GROK-16488).
// The tree gestures themselves are proven in hierarchical-and-combined-boolean-spec.ts.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Filter Panel — Hierarchical and Combined Boolean Filters: layout and project round-trips', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const isAmbientNoise = (text: string) =>
    /ProjectMeta\.publish/.test(text) ||
    /project_meta\.dart/.test(text) ||
    /could not be cloned/i.test(text) ||
    /Failed to load resource/.test(text) ||
    /favicon/.test(text);
  const pageErrors: string[] = [];
  const onConsole = (msg: import('@playwright/test').ConsoleMessage) => {
    if (msg.type() === 'error' && !isAmbientNoise(msg.text())) pageErrors.push(msg.text());
  };
  const onPageError = (err: Error) => {
    if (!isAmbientNoise(err.message)) pageErrors.push(`pageerror: ${err.message}`);
  };
  page.on('console', onConsole);
  page.on('pageerror', onPageError);

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, withFilterPanel: true});

  const rowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
  expect(rowCount).toBe(ROW_COUNT);
  const caucasianFemale = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const sex = df.col('SEX').toList();
    const race = df.col('RACE').toList();
    let n = 0;
    for (let i = 0; i < df.rowCount; i++) if (sex[i] === 'F' && race[i] === 'Caucasian') n++;
    return n;
  });
  expect(caucasianFemale).toBeGreaterThan(0);
  expect(caucasianFemale).toBeLessThan(rowCount);

  let hlTrueCount = 0;
  await softStep('Step 6-restore: establish SEX / RACE + Caucasian criterion for the layout save', async () => {
    await addHierarchicalCard(page);
    expect(await hierCaption(page)).toBe('SEX / RACE');
    await hierNode(page, ['F'], 'expand');
    await expect.poll(async () => (await hierNode(page, ['F'], 'read')).childCaptions,
      {message: 'the RACE children of F never rendered after the expander click',
        timeout: 10_000, intervals: [200, 400, 800]}).toContain('Caucasian');
    await hierNode(page, ['F', 'Caucasian'], 'toggle');
    await expect.poll(async () => trueCountOf(page),
      {message: 'the criterion about to be saved is not the Caucasian-female one derived from the raw '
        + `SEX / RACE columns (${caucasianFemale})`,
      timeout: 10_000, intervals: [200, 400, 800]}).toBe(caucasianFemale);
    hlTrueCount = await trueCountOf(page);
  });

  let hlLayoutId = '';
  await softStep('Step 10: Save the layout with the hierarchical criterion active', async () => {
    hlLayoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id;
    });
    expect(hlLayoutId).toBeTruthy();
    await expect.poll(async () => page.evaluate(async (id: string) => {
      try { return !!(await grok.dapi.layouts.find(id)); }
      catch (_) { return false; }
    }, hlLayoutId),
    {message: 'the saved layout cannot be fetched back from the server — saveLayout() stamps the id '
      + 'client-side before the round-trip, so a dapi.layouts.save that silently failed leaves the id '
      + 'just as non-empty and Step 11 would re-apply nothing',
    timeout: 20_000, intervals: [500, 1000, 2000]}).toBe(true);
  });

  await softStep('Step 9: GROK-16528 — reorder columns to RACE / SEX before the re-apply', async () => {
    await applyHierarchyState(page, {colNames: ['RACE', 'SEX'], allEnabled: true});
    await expect.poll(() => hierCaption(page), {timeout: 10_000, intervals: [100, 200, 400]}).toBe('RACE / SEX');
  });

  try {
    await softStep('Step 11: Layout round-trip (hierarchical) — close panel, re-apply saved layout', async () => {
      const beforeClose = await page.evaluate(async () => {
        const was = grok.shell.tv.dataFrame.filter.trueCount;
        const card = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
          .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.includes('/'));
        if (!card)
          throw new Error('no hierarchical filter card is painted in the Filter Panel — its caption is the only '
            + 'one carrying "/", and no card in the panel carries one');
        let clicked = false;
        for (const n of card.querySelectorAll('.d4-tree-view-node')) {
          const val = n.querySelector('.d4-hierarchical-filter-caption-value');
          if (val?.textContent?.trim() === 'Caucasian') {
            const cb = n.querySelector('input.d4-hierarchical-filter-checkbox') as HTMLElement | null;
            if (!cb)
              throw new Error('the RACE root node "Caucasian" carries no checkbox, so it could not be ticked');
            cb.click();
            clicked = true;
            break;
          }
        }
        if (!clicked) throw new Error('RACE root node "Caucasian" not found in the hierarchical card');
        return (window as any).__moved(() => grok.shell.tv.dataFrame.filter.trueCount, was, 900);
      });
      expect(beforeClose).toBeGreaterThan(0);
      expect(beforeClose).toBeLessThan(rowCount);

      await page.evaluate(closeFilterPanelInPage);
      const result = await page.evaluate(async (id: string) => {
        const w = window as any;
        await w.__poll(() => document.querySelectorAll('[name="viewer-Filters"]').length,
          (n: number) => n === 0, 1000, 50);
        const afterClose = grok.shell.tv.dataFrame.filter.trueCount;
        const panelsAfterClose = document.querySelectorAll('[name="viewer-Filters"]').length;

        const saved = await grok.dapi.layouts.find(id);
        grok.shell.tv.loadLayout(saved);
        grok.shell.tv.getFiltersGroup();
        await w.__poll(() => Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name'))
          .some((c) => (c.textContent ?? '').includes('/')), (there: boolean) => there, 3500, 50);
        const afterRestore = await w.__moved(
          () => grok.shell.tv.dataFrame.filter.trueCount, afterClose, 1000);

        let caption = '';
        for (const c of document.querySelectorAll('[name="viewer-Filters"] .d4-filter')) {
          const cn = c.querySelector('.d4-filter-column-name');
          if (cn && cn.textContent!.includes('/')) { caption = cn.textContent!.trim(); break; }
        }
        return {afterClose, panelsAfterClose, afterRestore, caption};
      }, hlLayoutId);
      expect(result.afterClose).toBe(rowCount);
      expect(result.afterClose).toBeGreaterThan(beforeClose);
      expect(result.panelsAfterClose).toBe(0);
      expect(result.afterRestore).toBe(hlTrueCount);
      expect(result.caption).toBe('SEX / RACE');
    });
  } finally {
    await page.evaluate(async (id: string) => {
      const saved = await grok.dapi.layouts.find(id);
      if (saved) await grok.dapi.layouts.delete(saved);
    }, hlLayoutId).catch(() => {});
  }

  await softStep('Step 12 / Step 13: a fresh demog view with SEX_bool, the Combined Boolean card auto-added', async () => {
    const result = await openDemogWithSexBool(page, datasetPath);
    expect(result.type).toBe('bool');
    expect(result.hasControl).toBe(true);
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('.d4-bool-combined-filter').waitFor({timeout: 10000});
    expect(await page.evaluate(() => document.querySelectorAll('.d4-bool-combined-filter').length)).toBe(1);
  });

  const probeProject = 'zz-grok16488-' + Date.now();
  let boolProjectId = '';
  try {
    await softStep('Step 17: GROK-16488 — save project, reopen, remove combined boolean card, count returns to full, no console error', async () => {
      await page.evaluate(async () => {
        const was = grok.shell.tv.dataFrame.filter.trueCount;
        const fg = grok.shell.tv.getFiltersGroup();
        for (const f of fg.filters) {
          if (f.filterType === 'bool-columns') {
            window.grok_GridFilterBase_ApplyState(f.dart ?? f, {'true': [true, true], 'false': [true, true], mode: 'OR'});
            grok.shell.tv.dataFrame.rows.requestFilter();
            break;
          }
        }
        await (window as any).__moved(() => grok.shell.tv.dataFrame.filter.trueCount, was, 500);
      });
      expect(await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount),
        'the combined boolean reset left the table filtered — the project would be saved narrowed')
        .toBe(ROW_COUNT);

      const saved = await saveProjectViaApi(page, probeProject);
      boolProjectId = saved.projectId;

      await page.evaluate(async (id: string) => {
        grok.shell.closeAll();
        const proj = await grok.dapi.projects.find(id);
        await proj.open();
        await (window as any).__tableReady(6000);
        grok.shell.tv.getFiltersGroup();
      }, boolProjectId);

      await page.locator('.d4-bool-combined-filter').waitFor({timeout: 15000});
      const afterReopen = await page.evaluate(() => ({
        boolCards: document.querySelectorAll('.d4-bool-combined-filter').length,
        filterPanel: document.querySelectorAll('[name="viewer-Filters"]').length,
      }));
      expect(afterReopen.filterPanel).toBeGreaterThan(0);
      expect(afterReopen.boolCards).toBe(1);

      const expectedRows = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const boolNames = df.columns.toList().filter((c: any) => c.type === 'bool').map((c: any) => c.name);
        return df.col(boolNames[0]).toList().filter((v: any) => v === true).length;
      });
      const preRemoval = await applyBoolState(page, {'true': [true, false], 'false': [false, false], mode: 'OR'});
      expect(expectedRows).toBeGreaterThan(0);
      expect(expectedRows).toBeLessThan(ROW_COUNT);
      expect(preRemoval).toBe(expectedRows);

      const errorsBefore = pageErrors.length;

      await page.evaluate(() => {
        const boolCard = document.querySelector('.d4-bool-combined-filter');
        let host: any = boolCard;
        while (host && !host.classList.contains('d4-filter')) host = host.parentElement;
        const x = host?.querySelector('[name="icon-times"]') as HTMLElement | null;
        if (!x) throw new Error('the combined boolean card carries no [name="icon-times"] remove icon — '
          + 'the removal this step is about never happened');
        x.click();
      });

      await expect.poll(async () => page.evaluate(() => ({
        boolCards: document.querySelectorAll('.d4-bool-combined-filter').length,
        trueCount: grok.shell.tv.dataFrame.filter.trueCount,
      })),
      {message: 'removing the combined boolean card did not both take the card away and release its '
        + `narrowing back to the full ${rowCount} rows`,
      timeout: 15_000, intervals: [200, 400, 800]}).toEqual({boolCards: 0, trueCount: rowCount});

      const removalSamples: string[] = [];
      for (let i = 0; i < 6; i++) {
        await page.waitForTimeout(400);
        removalSamples.push(JSON.stringify(pageErrors.slice(errorsBefore)));
      }
      expect(Array.from(new Set(removalSamples)),
        'GROK-16488 — removing the combined boolean card must not throw at any point of a 2.4s window '
        + `after the click, so an async throw arriving late is caught too; samples: ${removalSamples.join(' ; ')}`)
        .toEqual(['[]']);
    });
  } finally {
    await deleteProjectWithCleanup(page, {projectId: boolProjectId});
  }

  let boolLayoutId = '';
  try {
    await softStep('Step 18: Layout round-trip (combined boolean) — save active state, close, re-apply', async () => {
      const removeClicks = await page.evaluate(() => {
        let clicks = 0;
        for (const boolCard of document.querySelectorAll('.d4-bool-combined-filter')) {
          let host: any = boolCard;
          while (host && !host.classList.contains('d4-filter')) host = host.parentElement;
          const x = host?.querySelector('[name="icon-times"]') as HTMLElement | null;
          if (!x) throw new Error('a .d4-bool-combined-filter card carries no [name="icon-times"] remove icon');
          x.click();
          clicks++;
        }
        return clicks;
      });
      await expect.poll(async () => page.evaluate(() => ({
        cards: document.querySelectorAll('.d4-bool-combined-filter').length,
        filters: grok.shell.tv.getFiltersGroup().filters
          .filter((f: any) => f.filterType === 'bool-columns').length,
      })),
      {message: `the combined boolean card did not go away after ${removeClicks} remove-icon click(s) — ` +
        'the menu drive below would then be satisfied by a card it did not create',
      timeout: 15_000, intervals: [400, 800, 1500]}).toEqual({cards: 0, filters: 0});

      const beforeAdd = await page.evaluate(() => ({
        cards: document.querySelectorAll('.d4-bool-combined-filter').length,
        filters: grok.shell.tv.getFiltersGroup().filters
          .filter((f: any) => f.filterType === 'bool-columns').length,
      }));
      expect(beforeAdd,
        'a combined boolean filter survived the removal, so the Add Filter > Combined Boolean drive ' +
        `cannot be shown to have created anything; read ${JSON.stringify(beforeAdd)}`)
        .toEqual({cards: 0, filters: 0});

      await v.drivePanelMenuLeaf(page, 'Filters', 'Add Filter', 'Combined Boolean');
      await expect.poll(async () => page.evaluate(() => ({
        cards: document.querySelectorAll('.d4-bool-combined-filter').length,
        filters: grok.shell.tv.getFiltersGroup().filters
          .filter((f: any) => f.filterType === 'bool-columns').length,
      })),
      {message: 'Add Filter > Combined Boolean did not take the panel from 0 combined boolean cards to ' +
        'exactly 1 — either the leaf created nothing, or it fired twice and left two identical cards',
      timeout: 15_000, intervals: [400, 800, 1500]}).toEqual({cards: 1, filters: 1});
      const saved = await page.evaluate(async () => {
        let stage = 'toggling the combined boolean card before the save';
        try {
          const w = window as any;
          const countBeforeToggle = grok.shell.tv.dataFrame.filter.trueCount;
          const fg = grok.shell.tv.getFiltersGroup();
          for (const f of fg.filters) {
            if (f.filterType === 'bool-columns') {
              window.grok_GridFilterBase_ApplyState(f.dart ?? f, {'true': [true, false], 'false': [false, false], mode: 'OR'});
              grok.shell.tv.dataFrame.rows.requestFilter();
              break;
            }
          }
          const savedCount = await w.__moved(
            () => grok.shell.tv.dataFrame.filter.trueCount, countBeforeToggle, 700);

          stage = 'saving the layout';
          const layout = grok.shell.tv.saveLayout();
          await grok.dapi.layouts.save(layout);
          const id = layout.id;
          stage = 'fetching the saved layout back from the server';
          const serverFound = !!(await w.__findSaved(() => grok.dapi.layouts.find(id), 10_000));
          return {id, serverFound, savedCount};
        }
        catch (e: any) {
          throw new Error(`the combined boolean layout round-trip failed while ${stage}: ${e?.message ?? e}`);
        }
      });
      boolLayoutId = saved.id;
      await page.evaluate(closeFilterPanelInPage);
      const result = await page.evaluate(async (id: string) => {
        let stage = 'closing the Filter Panel';
        try {
          const w = window as any;
          await w.__poll(() => document.querySelectorAll('[name="viewer-Filters"]').length,
            (n: number) => n === 0, 1000, 50);
          const afterClose = grok.shell.tv.dataFrame.filter.trueCount;
          const cardsAfterClose = document.querySelectorAll('.d4-bool-combined-filter').length;

          stage = 're-applying the saved layout';
          const s = await grok.dapi.layouts.find(id);
          grok.shell.tv.loadLayout(s);
          grok.shell.tv.getFiltersGroup();
          await w.__poll(() => document.querySelectorAll('.d4-bool-combined-filter').length,
            (n: number) => n > 0, 3500, 50);
          const afterRestore = await w.__moved(
            () => grok.shell.tv.dataFrame.filter.trueCount, afterClose, 1000);
          const hasBoolCard = document.querySelectorAll('.d4-bool-combined-filter').length;
          return {afterClose, cardsAfterClose, afterRestore, hasBoolCard};
        }
        catch (e: any) {
          throw new Error(`the combined boolean layout round-trip failed while ${stage}: ${e?.message ?? e}`);
        }
      }, boolLayoutId);
      expect(saved.serverFound,
        'the saved layout could not be fetched back from the server — saveLayout() stamps the id '
        + 'client-side before the round-trip, so a dapi.layouts.save that silently failed would leave '
        + 'the re-apply below round-tripping nothing')
        .toBe(true);
      expect(saved.savedCount).toBeGreaterThan(0);
      expect(saved.savedCount).toBeLessThan(rowCount);
      expect(result.cardsAfterClose).toBe(0);
      expect(result.afterClose).toBe(rowCount);
      expect(result.hasBoolCard).toBeGreaterThan(0);
      expect(result.afterRestore).toBe(saved.savedCount);
    });
  } finally {
    await page.evaluate(async (id: string) => {
      if (!id) return;
      const saved = await grok.dapi.layouts.find(id);
      if (saved) await grok.dapi.layouts.delete(saved);
    }, boolLayoutId).catch(() => {});
    page.off('console', onConsole);
    page.off('pageerror', onPageError);
    await v.cleanupShell(page);
  }

  v.finishSpec();
});
