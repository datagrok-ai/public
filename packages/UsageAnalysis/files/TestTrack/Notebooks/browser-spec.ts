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
  }, null, {timeout: 60_000, polling: 250});
  await page.locator('.grok-gallery-search-bar').waitFor({timeout: 30_000});

  await page.waitForFunction(() => {
    return document.querySelectorAll('.grok-gallery-grid-item').length > 0 &&
      document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]').length > 0;
  }, null, {timeout: 60_000, polling: 250});
}

test('Notebooks / Browser (Integration): navigate, filter, context panel, apply-to, back, new', async ({page}) => {

  test.setTimeout(240_000);
  stepErrors.length = 0;

  let newNotebookId: string | null = null;

  await loginToDatagrok(page);
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });

  await page.waitForFunction(() => {
    try { return (window as any).grok.shell.v?.type !== 'notebooks'; } catch (e) { return false; }
  }, null, {timeout: 30_000, polling: 100});

  await softStep('Setup: open demog.csv as the active table', async () => {
    const info = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(df);
      await new Promise<void>((res) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); res(); });
        setTimeout(res, 3000);
      });
      return {rows: df.rowCount, cols: df.columns.length};
    });
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
    expect(info.rows).toBeGreaterThan(0);
  });

  await softStep('S1: open the Notebooks browser -> browser opens with a Demog card', async () => {
    await openNotebooksBrowser(page);

    const vtype = await page.evaluate(() => (window as any).grok.shell.v?.type);
    expect(vtype).toBe('notebooks');

    const demogCount = await page.evaluate(() => {
      const links = Array.from(document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]'));
      return links.filter((l) => (l.textContent ?? '').trim() === 'Demog').length;
    });
    expect(demogCount).toBeGreaterThan(0);
  });

  await softStep('S2: filter the gallery via search -> list narrows -> clear restores full list', async () => {
    const fullCount = await page.evaluate(() => {
      const t = document.querySelector('.grok-items-view-counts')?.textContent ?? '';
      return parseInt((t.split('/')[1] ?? t).trim(), 10) || document.querySelectorAll('.grok-gallery-grid-item').length;
    });

    const search = page.locator('.grok-gallery-search-bar input[placeholder*="notebook"]').first();
    await search.waitFor({timeout: 15_000});

    await search.click();
    const filtered = await page.evaluate(() => {
      const input = document.querySelector('.grok-gallery-search-bar input') as HTMLInputElement | null;
      return new Promise<{count: number; allMatch: boolean}>((resolve) => {
        if (!input) return resolve({count: 0, allMatch: false});
        input.focus();
        input.value = 'Demog';
        input.dispatchEvent(new Event('input', {bubbles: true}));
        let tries = 0;
        const tick = () => {
          const links = Array.from(document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]'));
          const count = document.querySelectorAll('.grok-gallery-grid-item').length;
          const allMatch = links.length > 0 && links.every((l) => (l.textContent ?? '').trim().toLowerCase().includes('demog'));
          if ((count > 0 && allMatch) || tries++ > 60) return resolve({count, allMatch});
          setTimeout(tick, 300);
        };
        setTimeout(tick, 300);
      });
    });
    const filteredItems = filtered.count;
    expect(filteredItems).toBeGreaterThan(0);
    expect(filteredItems).toBeLessThanOrEqual(fullCount);

    expect(filtered.allMatch).toBe(true);

    await page.evaluate(() => {
      const input = document.querySelector('.grok-gallery-search-bar input') as HTMLInputElement | null;
      if (input) { input.value = ''; input.dispatchEvent(new Event('input', {bubbles: true})); }
    });
    await page.waitForFunction((minItems) =>
      document.querySelectorAll('.grok-gallery-grid-item').length >= minItems,
    filteredItems, {timeout: 30_000, polling: 250});
    const restoredItems = await page.evaluate(() => document.querySelectorAll('.grok-gallery-grid-item').length);
    expect(restoredItems).toBeGreaterThanOrEqual(filteredItems);
  });

  await softStep('S3: select Demog card -> all 5 accordion panes present; Sharing renders (GROK-11693)', async () => {

    const selected = await page.evaluate(async () => {
      const links = Array.from(document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]'));
      const demog = links.find((l) => (l.textContent ?? '').trim() === 'Demog') as HTMLElement | undefined;
      if (!demog) return {ok: false, id: null};
      demog.click();
      const id = demog.getAttribute('data-link')?.split('/notebook/')[1] ?? null;
      for (let i = 0; i < 20; i++) {
        if (document.querySelector('[name="div-section--Sharing"]')) break;
        await new Promise((r) => setTimeout(r, 250));
      }
      return {ok: true, id};
    });
    expect(selected.ok).toBe(true);

    for (const section of ['div-section--Details', 'div-section--Actions', 'div-section--Activity', 'div-section--Sharing', 'div-section--Chats'])
      await expect(page.locator(`[name="${section}"]`).first()).toBeVisible();

    await page.locator('[name="div-section--Details"]').first().click();
    await page.waitForFunction(() => {
      const sec = document.querySelector('[name="div-section--Details"]');
      return (sec?.closest('.d4-accordion-pane')?.textContent ?? '').includes('Created');
    }, null, {timeout: 30_000, polling: 250});
    const detailsText = await page.evaluate(() => {
      const sec = document.querySelector('[name="div-section--Details"]');
      return sec?.closest('.d4-accordion-pane')?.textContent ?? '';
    });
    expect(detailsText).toContain('Created');

    const consoleErrors: string[] = [];
    const onConsole = (msg: import('@playwright/test').ConsoleMessage) => {
      if (msg.type() === 'error') consoleErrors.push(msg.text());
    };
    const onPageError = (err: Error) => consoleErrors.push(err.message);
    page.on('console', onConsole);
    page.on('pageerror', onPageError);

    await page.locator('[name="div-section--Sharing"]').first().click();

    await page.waitForFunction(() => {
      const sec = document.querySelector('[name="div-section--Sharing"]');
      const txt = (sec?.closest('.d4-accordion-pane')?.textContent ?? '').replace(/\s/g, '');
      return txt.length > 'Sharing'.length;
    }, null, {timeout: 30_000, polling: 250}).catch(() => {});

    page.off('console', onConsole);
    page.off('pageerror', onPageError);

    const sharingText = await page.evaluate(() => {
      const sec = document.querySelector('[name="div-section--Sharing"]');
      return sec?.closest('.d4-accordion-pane')?.textContent ?? '';
    });

    expect(sharingText.replace(/\s/g, '').length).toBeGreaterThan('Sharing'.length);

    const nullErr = consoleErrors.find((e) =>
      /NullError/i.test(e) && /notebook/i.test(e) || /method not found.*on null/i.test(e));
    expect(nullErr, `GROK-11693 regression: Sharing tab threw a notebook NullError: ${nullErr}`).toBeUndefined();
  });

  await softStep('S4: right-click Demog card -> entity context menu opens; Apply-to observed if applicable', async () => {

    await page.keyboard.press('Escape').catch(() => {});
    const demogCardLink = page.locator('.grok-gallery-grid .d4-link-label[data-link^="/notebook/"]', {hasText: 'Demog'}).first();
    await demogCardLink.scrollIntoViewIfNeeded().catch(() => {});
    await demogCardLink.click({button: 'right', timeout: 10_000}).catch(() => {});
    const menu = await page.evaluate(async () => {
      for (let i = 0; i < 30; i++) {
        if (document.querySelector('[name="div-Open"]') || document.querySelector('[name="div-Delete"]')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      const menuOpened = !!(document.querySelector('[name="div-Open"]') || document.querySelector('[name="div-Delete"]'));
      const leaf = document.querySelector('[name="div-Apply-to---Table"]');
      return {menuOpened, applyToPresent: !!document.querySelector('[name="div-Apply-to"]'), leafName: leaf ? 'div-Apply-to---Table' : null};
    });
    expect(menu.menuOpened, 'the notebook entity context menu (Open/Delete) should open on the card link').toBe(true);

    if (!menu.applyToPresent)
      console.warn('[OBSERVE notebooks.entity.get-applicable-cases] Apply-to group absent — the Demog ' +
        'notebook is not applicable to the open table on this instance (applicability is Dart-side, ' +
        'not JS-reachable). Recorded as an observation, not a hard gate.');

    if (menu.leafName) {
      await page.locator('[name="div-Apply-to---Table"]').first().click().catch(() => {});
      const opened = await page.evaluate(async () => {
        const grok = (window as any).grok;
        for (let i = 0; i < 30; i++) {
          try { if (grok.shell.v?.type === 'Notebook') return true; } catch (e) {  }
          await new Promise((r) => setTimeout(r, 1000));
        }
        return false;
      });
      if (!opened)
        console.warn('[OBSERVE] Apply-to did not open an HTML output view within ~30s — the apply path ' +
          'requires a live Jupyter container (notebooks.entity.apply is manual_only / environmental); ' +
          'covered as a behavioral remark, not a hard gate.');
    } else {

      await page.keyboard.press('Escape').catch(() => {});
    }
  });

  await softStep('S5: back-navigate to the Notebooks browser -> list restored (no active filter)', async () => {

    await openNotebooksBrowser(page);

    const restored = await page.evaluate(() => {
      const input = document.querySelector('.grok-gallery-search-bar input') as HTMLInputElement | null;
      const items = document.querySelectorAll('.grok-gallery-grid-item').length;
      return {searchValue: input?.value ?? '', items};
    });
    expect(restored.searchValue).toBe('');
    expect(restored.items).toBeGreaterThan(0);
  });

  await softStep('S6: click NEW NOTEBOOK in the ribbon -> a new notebook is created and opens', async () => {
    const btn = page.locator('[name="button-New-Notebook..."]').first();
    await btn.waitFor({timeout: 15_000});
    await btn.click();

    const res = await page.waitForFunction(() => {
      const grok = (window as any).grok;
      let v; let o;
      try { v = grok.shell.v; o = grok.shell.o; } catch (e) { return false; }
      if (v && v.type === 'Notebook' && o && o.id) return {oId: o.id};
      return false;
    }, null, {timeout: 90_000, polling: 250});
    const created = await res.jsonValue() as {oId: string | null};
    newNotebookId = created.oId;
    expect(newNotebookId).toBeTruthy();

    const persisted = await page.evaluate(async (id) => {
      const grok = (window as any).grok;
      if (!id) return false;
      for (let i = 0; i < 10; i++) {
        try { const nb = await grok.dapi.notebooks.find(id); if (nb && nb.id) return true; } catch (e) {  }
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    }, newNotebookId);
    expect(persisted).toBe(true);
  });

  await page.evaluate(async (id) => {
    const grok = (window as any).grok;
    try {
      if (id) {
        const nb = await grok.dapi.notebooks.find(id);
        if (nb) await grok.dapi.notebooks.delete(nb);
      }
      grok.shell.closeAll();
    } catch (e) {  }
  }, newNotebookId).catch(() => {});

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
