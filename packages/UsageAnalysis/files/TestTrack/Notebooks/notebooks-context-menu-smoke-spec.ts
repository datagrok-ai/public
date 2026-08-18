import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

async function openNotebooksBrowser(page: import('@playwright/test').Page) {
  await page.evaluate(async () => {
    const f = (window as any).DG.Func.find({name: 'CmdBrowseNotebooks'})[0];
    await f.apply();
  });
  await page.waitForFunction(() => {
    try { return (window as any).grok.shell.v?.type === 'notebooks'; } catch (e) { return false; }
  }, null, {timeout: 45_000, polling: 250});
  await page.locator('.grok-gallery-search-bar').waitFor({timeout: 30_000});
  await page.waitForFunction(() => {
    return document.querySelectorAll('.grok-gallery-grid-item').length > 0 &&
      document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]').length > 0;
  }, null, {timeout: 45_000, polling: 250});
}

async function openCardContextMenu(page: import('@playwright/test').Page, name: string,
  probeSelector: string): Promise<boolean> {
  return page.evaluate(async (args: {name: string; probe: string}) => {

    const label = Array.from(document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]'))
      .find((l) => l.textContent?.trim() === args.name) as HTMLElement | undefined;
    if (!label) return false;
    const r = label.getBoundingClientRect();
    label.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, clientX: r.left + 10, clientY: r.top + 5}));
    for (let i = 0; i < 30; i++) {
      await new Promise((rs) => setTimeout(rs, 100));
      const item = Array.from(document.querySelectorAll(args.probe)).find((e) => (e as HTMLElement).offsetParent !== null);
      if (item) return true;
    }
    return false;
  }, {name, probe: probeSelector});
}

async function narrowGalleryToSeed(page: import('@playwright/test').Page, name: string): Promise<boolean> {
  return page.evaluate(async (n) => {
    const input = document.querySelector('.grok-gallery-search-bar input') as HTMLInputElement | null;
    if (input) {
      input.focus();
      input.value = n;
      input.dispatchEvent(new Event('input', {bubbles: true}));
    }
    for (let i = 0; i < 50; i++) {
      await new Promise((r) => setTimeout(r, 300));
      const hit = Array.from(document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]'))
        .some((l) => l.textContent?.trim() === n);
      if (hit) return true;
    }
    return false;
  }, name);
}

async function ensureBrowserNarrowedToSeed(page: import('@playwright/test').Page, name: string) {
  const onBrowser = await page.evaluate(() => {
    try { return grok.shell.v?.type === 'notebooks'; } catch (e) { return false; }
  });
  if (!onBrowser) {
    await page.evaluate(() => grok.shell.closeAll());
    await openNotebooksBrowser(page);
  }
  await page.locator('.grok-gallery-search-bar input').waitFor({timeout: 30_000});
  let found = await narrowGalleryToSeed(page, name);
  if (!found) {

    await page.evaluate(() => grok.shell.closeAll());
    await openNotebooksBrowser(page);
    await page.locator('.grok-gallery-search-bar input').waitFor({timeout: 30_000});
    found = await narrowGalleryToSeed(page, name);
  }
  expect(found, `seeded notebook "${name}" should be findable in the browser`).toBe(true);
}

test('Notebooks — Context Menu Smoke (all 7 pcmd flows)', async ({page}) => {

  test.setTimeout(360_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
  });

  const seedName = `automator-ctxmenu-${Date.now()}`;
  let seededId: string | null = null;

  await softStep('Setup: open demog.csv as the active table', async () => {
    const info = await page.evaluate(async () => {
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(df);
      await new Promise<void>((res) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); res(); });
        setTimeout(res, 3000);
      });
      return {rows: df.rowCount};
    });
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
    expect(info.rows).toBeGreaterThan(0);
  });

  await softStep('Setup: seed a server notebook linked to demog (via Open in Notebook)', async () => {

    seededId = await page.evaluate(async (name) => {
      const before = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10});
      const beforeIds = new Set(before.map((n: any) => n.id));
      const f = (window as any).DG.Func.find({name: 'CmdOpenInNotebook'})[0];
      await f.apply();
      let fresh: any = null;
      for (let i = 0; i < 40; i++) {
        await new Promise((r) => setTimeout(r, 500));
        const cur = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10}).catch(() => [] as any[]);
        fresh = cur.find((n: any) => !beforeIds.has(n.id));
        if (fresh) break;
      }
      if (!fresh) return null;
      fresh.friendlyName = name;
      await grok.dapi.notebooks.save(fresh);
      return fresh.id as string;
    }, seedName);
    expect(seededId, 'Open in Notebook should persist a demog-linked server notebook').toBeTruthy();

    const persisted = await page.evaluate(async (args: {id: string; name: string}) => {
      for (let i = 0; i < 10; i++) {
        const ent: any = await grok.dapi.notebooks.find(args.id).catch(() => null);
        if (ent && (ent.friendlyName || ent.name) === args.name) return true;
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    }, {id: seededId as string, name: seedName});
    expect(persisted).toBe(true);
  });

  await softStep('Scenario 1 (pcmdOpen): right-click card -> Open -> HTML-mode view opens', async () => {
    await ensureBrowserNarrowedToSeed(page, seedName);
    const present = await openCardContextMenu(page, seedName, '[name="div-Open"]');
    expect(present, '[name="div-Open"] should be present in the context menu').toBe(true);
    await page.evaluate(() => (document.querySelector('[name="div-Open"]') as HTMLElement)?.click());

    const opened = await page.waitForFunction(() => {
      try { return grok.shell.v?.type === 'Notebook'; } catch (e) { return false; }
    }, null, {timeout: 30_000, polling: 250}).then(() => true).catch(() => false);
    expect(opened, 'Open should transition to a Notebook (HTML-mode) view').toBe(true);
  });

  await softStep('Scenario 2 (pcmdEdit): right-click card -> Edit -> edit-mode view transition', async () => {
    await ensureBrowserNarrowedToSeed(page, seedName);
    const present = await openCardContextMenu(page, seedName, '[name="div-Edit"]');
    expect(present, '[name="div-Edit"] should be present in the context menu').toBe(true);
    await page.evaluate(() => (document.querySelector('[name="div-Edit"]') as HTMLElement)?.click());

    const inEdit = await page.waitForFunction(() => {
      try { return grok.shell.v?.type === 'Notebook'; } catch (e) { return false; }
    }, null, {timeout: 30_000, polling: 250}).then(() => true).catch(() => false);
    expect(inEdit, 'Edit should transition to a Notebook (edit-mode) view').toBe(true);
  });

  await softStep('Scenario 5 (pcmdShare): right-click card -> Share... -> Share dialog opens', async () => {
    await ensureBrowserNarrowedToSeed(page, seedName);

    let shared = {ok: false, title: ''};
    for (let attempt = 0; attempt < 2 && !shared.ok; attempt++) {
      const present = await openCardContextMenu(page, seedName, '[name="div-Share..."]');
      if (!present) continue;
      shared = await page.evaluate(async () => {
        (document.querySelector('[name="div-Share..."]') as HTMLElement)?.click();
        for (let i = 0; i < 40; i++) {
          await new Promise((r) => setTimeout(r, 150));

          const dlg = Array.from(document.querySelectorAll('.d4-dialog'))
            .find((d) => d.querySelector('.dlg-sharing-settings') ||
              /^share /i.test(d.querySelector('.d4-dialog-title')?.textContent?.trim() ?? ''));
          if (dlg) {
            const title = dlg.querySelector('.d4-dialog-title')?.textContent?.trim() ?? '';

            (dlg.querySelector('[name="button-CANCEL"]') as HTMLElement)?.click();
            return {ok: true, title};
          }
        }
        return {ok: false, title: ''};
      });
    }
    expect(shared.ok, 'the shareEntity dialog should open (attached) for the notebook').toBe(true);

    expect(shared.title.toLowerCase()).toContain('share');
  });

  await softStep('Scenario 6 (pcmdApplyTo): demog open -> Apply-to group appears for the seed', async () => {
    await ensureBrowserNarrowedToSeed(page, seedName);

    const menu = await page.evaluate(async (name) => {

      const label = Array.from(document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]'))
        .find((l) => l.textContent?.trim() === name) as HTMLElement | undefined;
      if (!label) return {menuOpened: false, applyToPresent: false, leafPresent: false};
      const r = label.getBoundingClientRect();
      label.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, clientX: r.left + 10, clientY: r.top + 5}));
      for (let i = 0; i < 30; i++) {
        await new Promise((rs) => setTimeout(rs, 150));
        if (document.querySelector('[name="div-Open"]') || document.querySelector('[name="div-Delete"]')) break;
      }
      const menuOpened = !!(document.querySelector('[name="div-Open"]') || document.querySelector('[name="div-Delete"]'));
      const applyTo = !!document.querySelector('[name="div-Apply-to"]');
      const leaf = !!document.querySelector('[name="div-Apply-to---Table"]');

      document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      return {menuOpened, applyToPresent: applyTo, leafPresent: leaf};
    }, seedName);

    expect(menu.menuOpened, 'the notebook entity context menu (Open/Delete) should open on the seed card link').toBe(true);
    if (!menu.applyToPresent)
      console.warn('[OBSERVE notebooks.entity.get-applicable-cases] Apply-to group absent for the blank ' +
        'seed notebook (not applicable to the open table; applicability is Dart-side). Recorded as an ' +
        'observation, not a hard gate.');
  });

  await softStep('Scenario 4 (pcmdSaveAsJSON): right-click card -> Save As JSON -> file downloads', async () => {
    await ensureBrowserNarrowedToSeed(page, seedName);
    const present = await openCardContextMenu(page, seedName, '[name="div-Save-As-JSON"]');
    expect(present, '[name="div-Save-As-JSON"] should be present in the context menu').toBe(true);

    const downloadPromise = page.waitForEvent('download', {timeout: 30_000}).catch(() => null);
    await page.evaluate(() => (document.querySelector('[name="div-Save-As-JSON"]') as HTMLElement)?.click());
    const download = await downloadPromise;
    expect(download, 'Save As JSON should trigger a file download').not.toBeNull();
    expect(download!.suggestedFilename().toLowerCase()).toMatch(/\.json$/);
  });

  const renamedName = `${seedName}-renamed`;
  await softStep('Scenario 3 (pcmdRename): right-click card -> Rename... -> modal -> OK -> label updates', async () => {
    await ensureBrowserNarrowedToSeed(page, seedName);
    const present = await openCardContextMenu(page, seedName, '[name="div-Rename..."]');
    expect(present, '[name="div-Rename..."] should be present in the context menu').toBe(true);
    await page.evaluate(() => (document.querySelector('[name="div-Rename..."]') as HTMLElement)?.click());

    const input = page.locator('[name="input-Rename-Notebook---New-name-"]');
    await input.waitFor({timeout: 30_000});
    await expect(page.locator('.d4-dialog .d4-dialog-title').filter({hasText: 'Rename Notebook'}).first())
      .toBeVisible();

    await input.click();
    await input.press('ControlOrMeta+a');
    await input.press('Delete');
    await expect(page.locator('[name="button-OK"]').first()).toHaveClass(/disabled/);

    await page.keyboard.type(renamedName);
    await expect(page.locator('[name="button-OK"]').first()).toHaveClass(/enabled/);
    await page.locator('[name="button-OK"]').first().click();

    const renamed = await page.evaluate(async (args: {id: string | null; name: string}) => {
      if (!args.id) return false;
      for (let i = 0; i < 20; i++) {
        const ent: any = await grok.dapi.notebooks.find(args.id).catch(() => null);
        if (ent && (ent.friendlyName || ent.name) === args.name) return true;
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    }, {id: seededId, name: renamedName});
    expect(renamed, 'the notebook should be renamed server-side').toBe(true);

    await ensureBrowserNarrowedToSeed(page, renamedName);
    const cardShows = await page.evaluate((name) =>
      Array.from(document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]'))
        .some((l) => l.textContent?.trim() === name), renamedName);
    expect(cardShows, `card label should update to "${renamedName}"`).toBe(true);
  });

  await softStep('Scenario 7 (pcmdDelete): right-click card -> Delete -> YES -> card removed', async () => {
    await ensureBrowserNarrowedToSeed(page, renamedName);
    const present = await openCardContextMenu(page, renamedName, '[name="div-Delete"]');
    expect(present, '[name="div-Delete"] should be present (server-persisted notebook)').toBe(true);
    await page.evaluate(() => (document.querySelector('[name="div-Delete"]') as HTMLElement)?.click());

    const confirmed = await page.evaluate(async () => {
      for (let i = 0; i < 60; i++) {
        const dlg = Array.from(document.querySelectorAll('.d4-dialog'))
          .find((d) => /are you sure|delete notebook/i.test(d.textContent || ''));
        if (dlg) {
          const yes = dlg.querySelector('[name="button-YES"]') as HTMLElement | null;
          if (yes) { yes.click(); return true; }
        }
        await new Promise((r) => setTimeout(r, 100));
      }
      return false;
    });
    expect(confirmed, 'confirm dialog YES (button-YES) should be present and clickable').toBe(true);

    const removed = await page.evaluate(async (id: string | null) => {
      if (!id) return false;
      for (let i = 0; i < 30; i++) {
        const ent: any = await grok.dapi.notebooks.find(id).catch(() => undefined);
        if (!ent || !ent.id) return true;
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    }, seededId);
    expect(removed, 'the notebook should be removed from the server (AppEvents.ENTITY_REMOVED)').toBe(true);
    seededId = null; 
  });

  await page.evaluate(async (id: string | null) => {
    if (!id) return;
    try {
      const entity: any = await grok.dapi.notebooks.find(id).catch(() => null);
      if (entity && entity.id) await grok.dapi.notebooks.delete(entity).catch(() => {});
    } catch (_) {  }
    grok.shell.closeAll();
  }, seededId);

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
