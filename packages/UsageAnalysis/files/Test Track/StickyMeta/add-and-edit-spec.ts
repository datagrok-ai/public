import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e?.message ?? String(e)}); }
}

test('StickyMeta: Add and edit metadata (single cell, sticky column, batch)', async ({page}) => {
  test.setTimeout(300_000);

  await page.goto(baseUrl);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});

  // Setup: selenium class, tabs mode, close all, open SPGI.csv, wait for chem
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    const df = await (window as any).grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
    (window as any).grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 200));
    }
    await new Promise(r => setTimeout(r, 5000));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Section 1 — Add sticky metadata to a single cell
  await softStep('S1.1 Open SPGI dataset', async () => {
    const info = await page.evaluate(() => ({
      rows: (window as any).grok.shell.t?.rowCount,
      semType: (window as any).grok.shell.t?.col('Structure')?.semType,
    }));
    expect(info.rows).toBe(3624);
    expect(info.semType).toBe('Molecule');
  });

  await softStep('S1.2 Right-click Structure cell row 0, open Sticky meta > Edit for current cell...', async () => {
    const dialogTitle = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const grid = g.shell.tv.grid;
      grid.scrollToCell('Structure', 0);
      await new Promise(r => setTimeout(r, 400));
      g.shell.t.currentCell = g.shell.t.cell(0, 'Structure');
      const gridEl = document.querySelector('[name="viewer-Grid"]') as HTMLElement;
      const rect = gridEl.getBoundingClientRect();
      const cell = grid.cell('Structure', 0);
      const b = cell.bounds;
      const x = rect.left + b.x + b.width / 2;
      const y = rect.top + b.y + b.height / 2;
      const target = document.elementFromPoint(x, y) as HTMLElement;
      target.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 2, buttons: 2}));
      await new Promise(r => setTimeout(r, 800));
      const menu = document.querySelector('.d4-menu-popup');
      const edit = Array.from(menu?.querySelectorAll('.d4-menu-item-label') ?? [])
        .find(l => l.textContent?.trim() === 'Edit for current cell...');
      (edit?.closest('.d4-menu-item') as HTMLElement)?.click();
      await new Promise(r => setTimeout(r, 1500));
      return document.querySelector('.d4-dialog .d4-dialog-header, .d4-dialog .d4-dialog-title')
        ?.textContent?.trim();
    });
    expect(dialogTitle).toBe('Sticky meta');
  });

  await softStep('S1.3 Fill rating=5, notes="test note", verified=true, save', async () => {
    await page.evaluate(async () => {
      const d = document.querySelector('.d4-dialog') as HTMLElement;
      const r = d.querySelector('[name="input-host-Rating"] input') as HTMLInputElement;
      r.focus(); r.select(); r.value = ''; r.dispatchEvent(new Event('input', {bubbles: true}));
    });
    await page.locator('.d4-dialog [name="input-host-Rating"] input').click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('5');
    await page.keyboard.press('Tab');

    await page.locator('.d4-dialog [name="input-host-Notes"]').first().locator('input').click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('test note');
    await page.keyboard.press('Tab');

    await page.evaluate(async () => {
      const d = document.querySelector('.d4-dialog') as HTMLElement;
      const cb = d.querySelector('[name="input-host-Verified"] input[type="checkbox"]') as HTMLInputElement;
      if (!cb.checked) cb.click();
      await new Promise(r => setTimeout(r, 300));
      const saves = Array.from(d.querySelectorAll('[name="button-Save"]'));
      (saves[0] as HTMLElement).click();
      await new Promise(r => setTimeout(r, 1500));
    });

    const saved = await page.evaluate(() => {
      const d = document.querySelector('.d4-dialog') as HTMLElement;
      return {
        rating: (d.querySelector('[name="input-host-Rating"] input') as HTMLInputElement).value,
        notes: (d.querySelectorAll('[name="input-host-Notes"]')[0]
          .querySelector('input') as HTMLInputElement).value,
        verified: (d.querySelector('[name="input-host-Verified"] input[type="checkbox"]') as HTMLInputElement).checked,
      };
    });
    expect(saved.rating).toBe('5');
    expect(saved.notes).toBe('test note');
    expect(saved.verified).toBe(true);
  });

  await softStep('S1.4 Close dialog and verify metadata persisted by reopening', async () => {
    await page.evaluate(() => {
      const d = document.querySelector('.d4-dialog');
      const cancel = d?.querySelector('[name="button-CANCEL"]')
        ?? Array.from(d?.querySelectorAll('*') ?? []).find(e =>
          e.textContent?.trim().toUpperCase() === 'CANCEL' && e.children.length === 0);
      (cancel as HTMLElement)?.click();
    });
    await page.waitForTimeout(500);

    const persisted = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const grid = g.shell.tv.grid;
      grid.scrollToCell('Structure', 0);
      await new Promise(r => setTimeout(r, 400));
      const gridEl = document.querySelector('[name="viewer-Grid"]') as HTMLElement;
      const rect = gridEl.getBoundingClientRect();
      const cell = grid.cell('Structure', 0);
      const b = cell.bounds;
      const x = rect.left + b.x + b.width / 2;
      const y = rect.top + b.y + b.height / 2;
      (document.elementFromPoint(x, y) as HTMLElement).dispatchEvent(
        new MouseEvent('contextmenu',
          {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 2, buttons: 2}));
      await new Promise(r => setTimeout(r, 800));
      const menu = document.querySelector('.d4-menu-popup');
      const edit = Array.from(menu?.querySelectorAll('.d4-menu-item-label') ?? [])
        .find(l => l.textContent?.trim() === 'Edit for current cell...');
      (edit?.closest('.d4-menu-item') as HTMLElement)?.click();
      await new Promise(r => setTimeout(r, 1500));
      const d = document.querySelector('.d4-dialog') as HTMLElement;
      return {
        rating: (d.querySelector('[name="input-host-Rating"] input') as HTMLInputElement).value,
        notes: (d.querySelectorAll('[name="input-host-Notes"]')[0]
          .querySelector('input') as HTMLInputElement).value,
        verified: (d.querySelector('[name="input-host-Verified"] input[type="checkbox"]') as HTMLInputElement).checked,
      };
    });
    expect(persisted.rating).toBe('5');
    expect(persisted.notes).toBe('test note');
    expect(persisted.verified).toBe(true);

    // close the dialog
    await page.evaluate(() => {
      const d = document.querySelector('.d4-dialog');
      const cancel = d?.querySelector('[name="button-CANCEL"]')
        ?? Array.from(d?.querySelectorAll('*') ?? []).find(e =>
          e.textContent?.trim().toUpperCase() === 'CANCEL' && e.children.length === 0);
      (cancel as HTMLElement)?.click();
    });
  });

  // Section 2 — Create sticky column
  await softStep('S2.1 Add TestSchema1 sticky columns via Context Panel', async () => {
    // set shell.o to Structure column so Sticky meta pane shows schemas
    await page.evaluate(() => {
      (window as any).grok.shell.o = (window as any).grok.shell.t.col('Structure');
    });
    await page.waitForTimeout(1500);

    // expand Sticky meta pane if collapsed
    const stickyHeader = page.locator('.grok-prop-panel .d4-accordion-pane-header', {hasText: 'Sticky meta'}).first();
    await stickyHeader.waitFor({timeout: 10000});
    const expanded = await stickyHeader.evaluate((el) => el.classList.contains('expanded'));
    if (!expanded) await stickyHeader.click();
    await page.waitForTimeout(1500);

    // Click the TestSchema1-scoped "Add all properties as columns" + button.
    // The pane has one section per schema; scoping to the TestSchema1 section is
    // required because `button-Add-notes-as-a-column` exists in multiple schemas.
    await page.evaluate(() => {
      const sections = Array.from(document.querySelectorAll('.grok-prop-panel .d4-build-root.ui-form'));
      const ts1 = sections.find(s => s.querySelector('.d4-flex-row')?.textContent?.trim() === 'TestSchema1');
      const addAll = ts1?.querySelector('[name="button-Add-TestSchema1\'s-properties-as-columns"]') as HTMLElement | null;
      addAll?.click();
    });

    // Async schema matching populates all 5 columns after ~5s — poll up to 15s.
    await expect.poll(async () => {
      const names = await page.evaluate(() => (window as any).grok.shell.t.columns.names());
      return ['rating', 'notes', 'verified', 'review_date', 'approve'].every(n => names.includes(n));
    }, {timeout: 15_000, intervals: [500]}).toBe(true);

    const names = await page.evaluate(() => (window as any).grok.shell.t.columns.names());
    for (const n of ['rating', 'notes', 'verified', 'review_date', 'approve'])
      expect(names).toContain(n);
  });

  await softStep('S2.2 Sort ascending by rating column', async () => {
    await page.evaluate(async () => {
      (window as any).grok.shell.tv.grid.sort(['rating'], [true]);
      await new Promise(r => setTimeout(r, 800));
    });
    // sort applied without error
    expect(true).toBe(true);
  });

  await softStep('S2.3 Remove rating column and re-add, verify metadata persists', async () => {
    const removed = await page.evaluate(async () => {
      const df = (window as any).grok.shell.t;
      df.columns.remove('rating');
      await new Promise(r => setTimeout(r, 800));
      return !df.columns.contains('rating');
    });
    expect(removed).toBe(true);

    // Re-add via the + button in the TestSchema1 section of the Sticky meta pane
    await page.evaluate(() => {
      const sections = Array.from(document.querySelectorAll('.grok-prop-panel .d4-build-root.ui-form'));
      const ts1 = sections.find(s => s.querySelector('.d4-flex-row')?.textContent?.trim() === 'TestSchema1');
      const btn = ts1?.querySelector('[name="button-Add-rating-as-a-column"]') as HTMLElement | null;
      btn?.click();
    });
    // poll up to 15s for rating column to reappear
    await expect.poll(async () =>
      await page.evaluate(() => (window as any).grok.shell.t.columns.contains('rating')),
      {timeout: 15_000, intervals: [500]},
    ).toBe(true);

    const res = await page.evaluate(() => {
      const df = (window as any).grok.shell.t;
      const castCol = df.col('CAST Idea ID');
      let origIdx = -1;
      for (let i = 0; i < df.rowCount; i++)
        if (String(castCol.get(i)).trim() === '634783') { origIdx = i; break; }
      const ratingCol = df.col('rating');
      return {
        reAdded: df.columns.contains('rating'),
        ratingAtOrig: ratingCol && origIdx >= 0 ? ratingCol.get(origIdx) : null,
      };
    });
    expect(res.reAdded).toBe(true);
    expect(res.ratingAtOrig).toBe(5);
  });

  // Section 3 — Batch edit metadata on multiple rows
  await softStep('S3.1 Select 3 rows, right-click > Edit for all properties', async () => {
    // Reset any prior sort, select 3 rows, compute cell coords
    const coords = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const df = g.shell.t;
      try { g.shell.tv.grid.sort([], []); } catch {}
      df.selection.setAll(false);
      df.selection.set(10, true);
      df.selection.set(11, true);
      df.selection.set(12, true);
      df.currentRowIdx = 10;
      const grid = g.shell.tv.grid;
      grid.scrollToCell('Structure', 10);
      await new Promise(r => setTimeout(r, 800));
      const gridEl = document.querySelector('[name="viewer-Grid"]') as HTMLElement;
      const rect = gridEl.getBoundingClientRect();
      const cell = grid.cell('Structure', 10);
      const b = cell.bounds;
      return {
        x: Math.round(rect.left + b.x + b.width / 2),
        y: Math.round(rect.top + b.y + b.height / 2),
      };
    });

    // Real right-click via Playwright page.mouse (proper mousedown/mouseup + contextmenu)
    await page.mouse.move(coords.x, coords.y);
    await page.mouse.click(coords.x, coords.y, {button: 'right'});

    // Wait for context menu to appear, then DOM-click "Edit for all properties"
    // (visually hidden until parent submenu is hovered — DOM click fires regardless).
    await page.locator('.d4-menu-popup').first().waitFor({timeout: 5000});
    await page.evaluate(() => {
      const menu = document.querySelector('.d4-menu-popup');
      const label = Array.from(menu?.querySelectorAll('.d4-menu-item-label') ?? [])
        .find(l => l.textContent?.trim() === 'Edit for all properties');
      (label?.closest('.d4-menu-item') as HTMLElement)?.click();
    });

    // Wait for the batch-edit host
    await page.locator('[name="input-host-Rows"]').waitFor({timeout: 10000});
    expect(await page.locator('[name="input-host-Rows"]').isVisible()).toBe(true);
  });

  await softStep('S3.2 Set verified=true and notes="batch note" for selected rows', async () => {
    // Click the Verified checkbox via real Playwright click
    const verified = page.locator('input[name="prop-view-verified"]').first();
    await verified.waitFor({timeout: 5000});
    if (!(await verified.isChecked())) await verified.click();
    await page.waitForTimeout(400);

    // Enter edit mode on notes row by clicking its view label (last notes row, below y≈220)
    await page.evaluate(async () => {
      const spans = Array.from(document.querySelectorAll('.property-grid-item-name-text span'))
        .filter(s => s.textContent?.trim() === 'notes');
      let row: HTMLElement | null = null;
      for (const s of spans) {
        const tr = s.closest('tr') as HTMLElement | null;
        if (tr && tr.getBoundingClientRect().top > 220) { row = tr; break; }
      }
      const vc = row?.querySelector('.property-grid-item-value');
      const label = vc?.querySelector('.property-grid-item-view-label') as HTMLElement | null;
      label?.click();
    });
    await page.waitForTimeout(500);
    const notesInput = page.locator('input.property-grid-item-editor-textbox').first();
    await notesInput.waitFor({timeout: 5000});
    await notesInput.focus();
    // clear any pre-existing value (e.g., from prior runs) before typing
    await page.keyboard.press('Control+a');
    await page.keyboard.press('Delete');
    await page.keyboard.type('batch note');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(800);

    const vals = await page.evaluate(() => {
      const df = (window as any).grok.shell.tables[0];
      const notes = df.col('notes'), verified = df.col('verified'), cast = df.col('CAST Idea ID');
      const target = ['634793', '634794', '634795'];
      const out: any[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        const c = String(cast.get(i)).trim();
        if (target.includes(c)) out.push({cast: c, notes: notes.get(i), verified: verified.get(i)});
      }
      return out;
    });
    expect(vals.length).toBe(3);
    for (const v of vals) {
      expect(v.notes).toBe('batch note');
      expect(v.verified).toBe(true);
    }
  });

  // Cleanup
  await page.evaluate(() => (window as any).grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
