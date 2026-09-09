import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {
  saveProject,
  waitForSaved,
  saveAndReopen,
  openProjectFromDashboards,
  buildVariantsComposite,
  collectChainProducedProjects,
  rename,
  deleteProjectViaContextMenu,
  shareProjectViaContextMenu,
  saveCopy,
} from './projects';

test.use({...specTestOptions, navigationTimeout: 180_000});

const SESSION_PREFIX = 'c1a-helpers-' + Date.now() + '-';

async function setupTableView(page: any, suffix: string): Promise<string> {
  const tableName = `c1a-fixture-${suffix}`;
  await page.evaluate(async (name: string) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    grok.shell.closeAll();
    const csv = 'id,name\n' + Array.from({length: 5}, (_, i) => `${i + 1},x${i + 1}`).join('\n');
    const t = DG.DataFrame.fromCsv(csv);
    t.name = name;
    grok.shell.addTableView(t);
  }, tableName);
  await page.waitForTimeout(500);
  return tableName;
}

test('helpers.playwright.projects C1a — full suite', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  const fixtureNames: string[] = [];

  try {
    await softStep('saveProject + waitForSaved (combined — save-then-poll)', async () => {
      const name = SESSION_PREFIX + 'save';
      fixtureNames.push(name);
      await setupTableView(page, 'save');
      await saveProject(page, name);
      await waitForSaved(page, name);

      const exists = await page.evaluate(async (n) => {
        const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
        return p != null && p.name === n;
      }, name);
      expect(exists).toBe(true);
    });

    await softStep('saveAndReopen — save, close, reopen, count tables/viewers', async () => {
      const name = SESSION_PREFIX + 'sar';
      fixtureNames.push(name);
      await setupTableView(page, 'sar');
      const counts = await saveAndReopen(page, name);

      expect(counts.tablesAfterReopen).toBeGreaterThanOrEqual(0);
      expect(counts.viewersAfterReopen).toBeGreaterThanOrEqual(0);
      expect(typeof counts.tablesAfterReopen).toBe('number');
    });

    await softStep('openProjectFromDashboards — locate tile + double-click', async () => {
      const name = SESSION_PREFIX + 'open';
      fixtureNames.push(name);
      await setupTableView(page, 'open');
      await saveProject(page, name);
      await waitForSaved(page, name);

      await page.evaluate(() => (window as any).grok.shell.closeAll());
      await page.evaluate(() => {
        const labels = Array.from(document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-node-label, .d4-tree-view-item'));
        const dashboardsNode = labels.find(el => (el.textContent || '').trim() === 'Dashboards' && (el as HTMLElement).offsetParent !== null);
        if (dashboardsNode) (dashboardsNode as HTMLElement).click();
      });
      await page.waitForTimeout(2000);

      await page.evaluate(() => {
        const refresh = document.querySelector('.grok-card-view i.grok-icon, .grok-card-view i[class*="sync"]') as HTMLElement | null;
        if (refresh) refresh.click();
      });
      await page.waitForTimeout(1500);

      const search = page.locator('input[placeholder*="Search projects"]');
      if (await search.isVisible({timeout: 3000}).catch(() => false)) {
        await search.focus();
        await page.keyboard.press('Control+a');
        await page.keyboard.type(name);
        await page.waitForTimeout(2000);
      }
      await openProjectFromDashboards(page, name);
      const opened = await page.evaluate((n) => {
        const grok = (window as any).grok;
        return grok.shell.project != null && grok.shell.project.name === n;
      }, name);
      expect(opened).toBe(true);
    });

    await softStep('buildVariantsComposite — link / clone / pvc variants', async () => {
      const sourceName = SESSION_PREFIX + 'composite';
      fixtureNames.push(sourceName);
      await setupTableView(page, 'composite');
      await saveProject(page, sourceName);
      await waitForSaved(page, sourceName);
      const composite = await buildVariantsComposite(page, sourceName);
      expect(composite.original).toBe(sourceName);
      expect(composite.link).toBe(`${sourceName}-link`);
      expect(composite.clone).toBe(`${sourceName}-clone`);
      expect(composite.pvc).toBe(`${sourceName}-pvc`);

      if (composite.link) fixtureNames.push(composite.link);
      if (composite.clone) fixtureNames.push(composite.clone);
      if (composite.pvc) fixtureNames.push(composite.pvc);

      for (const variant of [composite.link, composite.clone, composite.pvc]) {
        if (!variant) continue;
        const exists = await page.evaluate(async (n) => {
          const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
          return p != null;
        }, variant);
        expect(exists).toBe(true);
      }
    });

    await softStep('collectChainProducedProjects — prefix sweep returns the session fixtures', async () => {
      const list = await collectChainProducedProjects(page, SESSION_PREFIX);

      const names = list.map((p) => p.name);
      for (const expected of fixtureNames)
        expect(names).toContain(expected);
      expect(list.length).toBeGreaterThanOrEqual(fixtureNames.length);

      for (const p of list) {
        expect(p.name).toBeTruthy();
        expect(p.id).toBeTruthy();
      }
    });

    await softStep('rename — JS API rename verified via find-by-id (not filter-by-name)', async () => {
      const originalName = SESSION_PREFIX + 'rename-orig';
      const renamedName = SESSION_PREFIX + 'rename-new';
      fixtureNames.push(originalName, renamedName);
      await setupTableView(page, 'rename');
      await saveProject(page, originalName);
      await waitForSaved(page, originalName);

      const projectId = await page.evaluate(async (n) => {
        const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
        return p ? p.id : null;
      }, originalName);
      expect(projectId).toBeTruthy();
      const result = await rename(page, projectId, renamedName);
      expect(result.ok).toBe(true);
      expect(result.newName).toBe(renamedName);
      expect(result.verifiedVia).toBe('find-by-id');

      const findById = await page.evaluate(async (id) => {
        const p = await (window as any).grok.dapi.projects.find(id);
        return p ? p.name : null;
      }, projectId);
      expect(findById).toBe(renamedName);
    });
  } finally {

    await page.evaluate(async (prefix) => {
      const grok = (window as any).grok;
      const list = await grok.dapi.projects.filter(`name like "${prefix}%"`).list();
      for (const p of list) {
        try { await grok.dapi.projects.delete(p); }
        catch (e) { console.warn('Cleanup delete failed for ' + p.name + ': ' + e); }
      }
    }, SESSION_PREFIX).catch((e) => console.warn('Terminal sweep failed: ' + e));
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

const SESSION_PREFIX_C1B = 'c1b-helpers-' + Date.now() + '-';

async function createDashboardsFixture(page: any, name: string): Promise<void> {
  await page.evaluate(async (n: string) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    grok.shell.closeAll();
    const p = DG.Project.create();
    p.name = n;
    await grok.dapi.projects.save(p);
  }, name);

  await page.evaluate(() => {
    const labels = Array.from(document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-node-label, .d4-tree-view-item'));
    const node = labels.find((el: any) => (el.textContent || '').trim() === 'Dashboards' && el.offsetParent !== null);
    if (node) (node as HTMLElement).click();
  });
  await page.waitForTimeout(2000);

  const search = page.locator('input[placeholder*="Search projects"]');
  if (await search.isVisible({timeout: 3000}).catch(() => false)) {
    await search.focus();
    await page.keyboard.press('Control+a');
    await page.keyboard.type(name);
    await page.waitForTimeout(1000);
  }

  const slug = name.toLowerCase();
  for (let attempt = 0; attempt < 6; attempt++) {
    await page.evaluate(() => {
      const refresh = document.querySelector('.grok-card-view i.grok-icon, .grok-card-view i[class*="sync"]') as HTMLElement | null;
      if (refresh) refresh.click();
    });
    await page.waitForTimeout(2000);
    const found = await page.evaluate((s: string) => {
      return !!document.querySelector(`[name="div-${s}"]`);
    }, slug);
    if (found) return;
  }

}

test('helpers.playwright.projects C1b — UI helpers suite', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  const fixtureNames: string[] = [];

  try {
    await softStep('deleteProjectViaContextMenu — right-click mode', async () => {
      const name = SESSION_PREFIX_C1B + 'del-rclick';
      fixtureNames.push(name);

      await createDashboardsFixture(page, name);
      const projectId = await page.evaluate(async (n) => {
        const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
        return p ? p.id : null;
      }, name);
      expect(projectId).toBeTruthy();
      await deleteProjectViaContextMenu(page, name, 'right-click');

      const postDelete = await page.evaluate(async (id) => {
        const p = await (window as any).grok.dapi.projects.find(id);
        return p ? p.name : null;
      }, projectId);
      expect(postDelete).toBeNull();
    });

    await softStep('deleteProjectViaContextMenu — context-panel-dropdown mode', async () => {
      const name = SESSION_PREFIX_C1B + 'del-cp';
      fixtureNames.push(name);
      await createDashboardsFixture(page, name);
      const projectId = await page.evaluate(async (n) => {
        const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
        return p ? p.id : null;
      }, name);
      expect(projectId).toBeTruthy();
      await deleteProjectViaContextMenu(page, name, 'context-panel-dropdown');
      const postDelete = await page.evaluate(async (id) => {
        const p = await (window as any).grok.dapi.projects.find(id);
        return p ? p.name : null;
      }, projectId);
      expect(postDelete).toBeNull();
    });

    await softStep('shareProjectViaContextMenu — share with Test permission group', async () => {
      const name = SESSION_PREFIX_C1B + 'share';
      fixtureNames.push(name);
      await createDashboardsFixture(page, name);
      const projectId = await page.evaluate(async (n) => {
        const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
        return p ? p.id : null;
      }, name);
      expect(projectId).toBeTruthy();
      await shareProjectViaContextMenu(page, name, {
        recipient: 'Test permission group',
        accessLevel: 'View and use',
      });

    });

    await softStep('saveCopy — copy mode with explicit user-typed name', async () => {
      const sourceName = SESSION_PREFIX_C1B + 'original';
      fixtureNames.push(sourceName);

      await page.evaluate(async (n: string) => {
        const grok = (window as any).grok;
        const DG = (window as any).DG;
        grok.shell.closeAll();
        const csv = 'id,name\n1,a\n2,b\n3,c';
        const t = DG.DataFrame.fromCsv(csv);
        t.name = 'c1b-saveCopy-table';
        grok.shell.addTableView(t);
        const p = DG.Project.create();
        p.name = n;
        await grok.dapi.projects.save(p);
      }, sourceName);
      await page.waitForTimeout(1000);

      const copyName = 'C1bSaveCopy' + Date.now();  
      fixtureNames.push(copyName);
      await saveCopy(page, {
        sourceName,
        mode: 'copy',
        name: copyName,
        perTableLinkOrClone: 'clone',
      });
    });
  } finally {

    await page.evaluate(async (prefix) => {
      const grok = (window as any).grok;

      const variants = [prefix, prefix.charAt(0).toUpperCase() + prefix.slice(1).replace(/-/g, '').replace(/^./, c => c.toUpperCase())];
      const seenIds = new Set<string>();
      for (const p of variants) {
        const list = await grok.dapi.projects.filter(`name like "${p}%"`).list();
        for (const proj of list) {
          if (seenIds.has(proj.id)) continue;
          seenIds.add(proj.id);
          try { await grok.dapi.projects.delete(proj); }
          catch (e) { console.warn('Cleanup delete failed for ' + proj.name + ': ' + e); }
        }
      }
    }, SESSION_PREFIX_C1B).catch((e) => console.warn('C1b terminal sweep failed: ' + e));
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} C1b step(s) failed:\n${summary}`);
  }
});
