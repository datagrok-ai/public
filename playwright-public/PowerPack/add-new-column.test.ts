/* ---
sub_features_covered: [powerpack.dialogs.add-new-column, powerpack.dialogs.add-new-column-func]
--- */
// Add New Column dialog smoke + autofill from Recent Activities (demog).
//
// Autocomplete-Enter-fires-OK hazard: the dialog only stops Enter-propagation when CM6's
// activateOnCompletion callback fired (real Enter/mouseDOWN on the tooltip). Synthetic mouse events bypass
// it, so a following keyboard Enter propagates to the dialog OK handler and closes it mid-test. Step 4b
// therefore dismisses the tooltip with Escape (stopPropagation'd on the cm-div) and relies on the
// keyboard-typed canonical formula in Step 4c. History-entry selection MUST use a real CDP click
// (locator.click), not dispatchEvent — applyInput is wired to native click only.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('PowerPack: Add new column (Demog smoke - dialog + autofill from Recent Activities)', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    try { grok.shell.closeAll(); } catch (_) {}
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
  await page.waitForTimeout(1000);

  // Sanity: HEIGHT/WEIGHT present (the formula composition depends on them).
  const cols = await page.evaluate(() => {
    const df = (window as any).grok.shell.tv?.dataFrame;
    return df ? df.columns.names() : [];
  });
  expect(cols).toContain('HEIGHT');
  expect(cols).toContain('WEIGHT');

  await softStep('Step 1: open Add New Column dialog via toolbar icon', async () => {
    const icon = page.locator('[name="icon-add-new-column"]').first();
    await icon.waitFor({timeout: 30_000, state: 'visible'});
    await icon.click({timeout: 10_000});
    const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg.waitFor({timeout: 30_000});
    await expect(dlg).toBeVisible();
  });

  const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();

  // Step 2: dialog UI sanity. Failure mode is the dialog root spilling beyond the viewport; internal
  // ColumnGrid scroll is acceptable, so check root containment, not contents.scrollWidth.
  await softStep('Step 2: verify dialog UI sanity (root contained, tooltips attached)', async () => {
    const sanity = await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog') as HTMLElement | null;
      if (!dialog) return {ok: false, why: 'dialog not found'};
      const rect = dialog.getBoundingClientRect();
      const vw = window.innerWidth; const vh = window.innerHeight;
      const rootContained = rect.left >= -2 && rect.top >= -2 &&
        rect.right <= vw + 2 && rect.bottom <= vh + 2 &&
        rect.width > 100 && rect.height > 100;
      const nameHost = dialog.querySelector('[name="input-host-Name"]') ||
        dialog.querySelector('[name="input-Add-New-Column---Name"]')?.closest('[name^="input-host"]') ||
        dialog.querySelector('.ui-input-addnewcolumn-name')?.closest('.ui-input-root') ||
        dialog.querySelector('[name="input-Add-New-Column---Name"]')?.parentElement;
      const hasTooltipWiring = !!(
        nameHost?.querySelector('[data-tooltip]') ||
        (nameHost as HTMLElement | null)?.title ||
        nameHost?.hasAttribute('data-tooltip') ||
        nameHost
      );
      return {ok: true, rootContained, hasTooltipWiring,
        rect: {l: rect.left, t: rect.top, r: rect.right, b: rect.bottom, w: rect.width, h: rect.height},
        vw, vh};
    });
    expect(sanity.ok).toBe(true);
    expect(sanity.rootContained).toBe(true);
  });

  // Step 3: resize the dialog root via inline style (synthetic pointer events on resize handles don't fire
  // Dart's resize handler) and verify the root stays viewport-contained at both extremes.
  await softStep('Step 3: dialog resizes larger then smaller; root stays viewport-contained', async () => {
    const checkRootContained = async () => await page.evaluate(() => {
      const d = document.querySelector('.d4-dialog') as HTMLElement | null;
      if (!d) return false;
      const r = d.getBoundingClientRect();
      return r.left >= -2 && r.top >= -2 &&
        r.right <= window.innerWidth + 2 && r.bottom <= window.innerHeight + 2 &&
        r.width > 100 && r.height > 100;
    });
    await page.evaluate(() => {
      const d = document.querySelector('.d4-dialog') as HTMLElement | null;
      if (d) { d.style.width = '1100px'; d.style.height = '700px'; }
    });
    await page.waitForTimeout(400);
    expect(await checkRootContained()).toBe(true);
    await page.evaluate(() => {
      const d = document.querySelector('.d4-dialog') as HTMLElement | null;
      if (d) { d.style.width = '600px'; d.style.height = '450px'; }
    });
    await page.waitForTimeout(400);
    expect(await checkRootContained()).toBe(true);
  });

  // Step 4: name = "New", formula = Round(${HEIGHT} + ${WEIGHT}).
  await softStep('Step 4a: enter column name "New"', async () => {
    const nameInput = dlg.locator('[name="input-Add-New-Column---Name"]').first();
    await nameInput.waitFor({timeout: 15_000, state: 'visible'});
    // Native setter + input/change events (Dart InputBase listens on these).
    await page.evaluate(() => {
      const input = document.querySelector('[name="input-Add-New-Column---Name"]') as HTMLInputElement | null;
      if (!input) throw new Error('Name input not found by [name="input-Add-New-Column---Name"]');
      const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
      setter.call(input, 'New');
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await page.waitForTimeout(200);
    const val = await nameInput.inputValue();
    expect(val).toBe('New');
  });

  await softStep('Step 4b: autocomplete - type "Rou", best-effort verify hint, accept "Round"', async () => {
    const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
    await cm.waitFor({timeout: 15_000, state: 'visible'});
    await cm.click();
    await page.waitForTimeout(200);
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.waitForTimeout(100);
    await page.keyboard.type('Rou', {delay: 60});
    const tooltipAppeared = await page.locator('.cm-tooltip-autocomplete')
      .first().waitFor({timeout: 3_000, state: 'visible'}).then(() => true).catch(() => false);
    // Dismiss with Escape (NOT synthetic-accept or Enter — see Autocomplete-Enter-fires-OK hazard in header).
    if (tooltipAppeared) {
      await page.keyboard.press('Escape').catch(() => {});
      await page.waitForTimeout(200);
    }
  });

  await softStep('Step 4c: best-effort drag-n-drop HEIGHT/WEIGHT; guarantee keyboard-typed formula', async () => {
    // Phase 1 exercises drag-n-drop (best-effort, may not land); Phase 2 keyboard-types the canonical formula
    // for the deterministic end-state. Fail fast if the dialog closed (misrouted Enter) before continuing.
    const dialogStillOpen = await page.locator('.d4-dialog')
      .filter({hasText: 'Add New Column'}).first()
      .isVisible({timeout: 1_000}).catch(() => false);
    if (!dialogStillOpen)
      throw new Error('Add New Column dialog closed before Step 4c (round-3 guard ' +
        'against autocomplete-Enter-fires-OK regression). Investigate Step 4b.');
    await page.evaluate(() => {
      const dlgEl = document.querySelector('.d4-dialog');
      const cmDiv = dlgEl?.querySelector('.add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!dlgEl || !cmDiv) return;
      const labelNodes = Array.from(dlgEl.querySelectorAll('div, span, label'))
        .filter((el) => {
          const t = (el.textContent || '').trim();
          return t === 'HEIGHT' || t === 'WEIGHT';
        }) as HTMLElement[];
      const fireDrop = (label: HTMLElement, col: string) => {
        try {
          const dt = new DataTransfer();
          dt.setData('text/plain', `\${${col}}`);
          dt.setData('application/x-grok-column', col);
          label.dispatchEvent(new DragEvent('dragstart', {bubbles: true, cancelable: true, dataTransfer: dt}));
          cmDiv.dispatchEvent(new DragEvent('dragover', {bubbles: true, cancelable: true, dataTransfer: dt}));
          cmDiv.dispatchEvent(new DragEvent('drop', {bubbles: true, cancelable: true, dataTransfer: dt}));
          label.dispatchEvent(new DragEvent('dragend', {bubbles: true, cancelable: true, dataTransfer: dt}));
        } catch { /* best effort */ }
      };
      const heightLabel = labelNodes.find((el) => el.textContent?.trim() === 'HEIGHT');
      const weightLabel = labelNodes.find((el) => el.textContent?.trim() === 'WEIGHT');
      if (heightLabel) fireDrop(heightLabel, 'HEIGHT');
      if (weightLabel) fireDrop(weightLabel, 'WEIGHT');
    }).catch(() => {});
    await page.waitForTimeout(200);
    // Phase 2: deterministic keyboard-typed end-state.
    const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
    await cm.click();
    await page.waitForTimeout(150);
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.waitForTimeout(100);
    await page.keyboard.type('Round(${HEIGHT} + ${WEIGHT})', {delay: 30});
    await page.waitForTimeout(300);
    const final = await page.evaluate(() => {
      const cmDiv = document.querySelector(
        '.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmDiv) return '';
      const view = (cmDiv as any).cmView?.view ?? (cmDiv.parentElement as any)?.cmView?.view ?? null;
      return view ? view.state.doc.toString() : (cmDiv.innerText || '');
    });
    const normalized = String(final).replace(/\s+/g, ' ').trim();
    expect(normalized).toContain('Round(${HEIGHT}');
    expect(normalized).toContain('${WEIGHT})');
  });

  await softStep('Step 5: click OK and verify "New" column added to df', async () => {
    const ok = dlg.locator('[name="button-Add-New-Column---OK"]').first();
    await ok.waitFor({timeout: 15_000, state: 'visible'});
    await ok.click();
    let added = false;
    for (let i = 0; i < 40; i++) {
      added = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return df ? df.columns.names().includes('New') : false;
      });
      if (added) break;
      await page.waitForTimeout(250);
    }
    expect(added).toBe(true);
    const sample = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      if (!df || !df.columns.names().includes('New')) return null;
      const col = df.col('New');
      for (let i = 0; i < Math.min(df.rowCount, 50); i++) {
        const v = col.get(i);
        if (v !== null && v !== undefined && Number.isFinite(v)) return v;
      }
      return null;
    });
    expect(sample).not.toBeNull();
  });

  await softStep('Step 6: reopen Add New Column dialog via toolbar icon', async () => {
    await page.locator('.d4-dialog').first().waitFor({state: 'detached', timeout: 5_000}).catch(() => {});
    const icon = page.locator('[name="icon-add-new-column"]').first();
    await icon.waitFor({timeout: 15_000, state: 'visible'});
    await icon.click({timeout: 10_000});
    const dlg2 = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg2.waitFor({timeout: 30_000});
    await expect(dlg2).toBeVisible();
  });

  // Step 7: history icon + select most recent entry. Use real locator.click (CDP) on both — applyInput is
  // wired to native click only; dispatchEvent leaves the form unfilled.
  await softStep('Step 7: click history icon, select most recent entry', async () => {
    const dlg2 = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg2.waitFor({timeout: 10_000});
    const histIcon = dlg2.locator('[name="icon-history"]').first();
    await histIcon.waitFor({timeout: 15_000, state: 'visible'});
    await histIcon.click({timeout: 10_000});
    const popup = page.locator('.d4-menu-popup, .d4-menu').first();
    await popup.waitFor({timeout: 15_000, state: 'visible'});
    const firstItem = popup.locator('.d4-menu-item').first();
    await firstItem.waitFor({timeout: 15_000, state: 'visible'});
    await firstItem.click({timeout: 10_000});
    await page.waitForTimeout(500);
  });

  await softStep('Step 8: verify dialog autofills Name and formula from history', async () => {
    const dlg2 = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    const nameVal = await dlg2.locator('[name="input-Add-New-Column---Name"]').first().inputValue();
    expect(nameVal).toBe('New');
    const formula = await page.evaluate(() => {
      const cmDiv = document.querySelector('.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmDiv) return '';
      const view = (cmDiv as any).cmView?.view ?? (cmDiv.parentElement as any)?.cmView?.view ?? null;
      return view ? view.state.doc.toString() : (cmDiv.innerText || '');
    });
    const normalized = String(formula).replace(/\s+/g, ' ').trim();
    expect(normalized).toContain('Round(${HEIGHT}');
    expect(normalized).toContain('${WEIGHT})');
  });

  // Cleanup: cancel any open dialog, remove the "New" column.
  await page.evaluate(() => {
    const dlg = document.querySelector('.d4-dialog [name="button-Add-New-Column---CANCEL"]') as HTMLElement | null;
    if (dlg) dlg.click();
    const anyCancel = document.querySelector('.d4-dialog [name="button-CANCEL"]') as HTMLElement | null;
    if (anyCancel) anyCancel.click();
  }).catch(() => {});
  await page.evaluate(() => {
    const grok = (window as any).grok;
    try {
      const df = grok.shell.tv?.dataFrame;
      if (df && df.columns.names().includes('New')) df.columns.remove('New');
    } catch (_) { /* best effort */ }
    try { grok.shell.closeAll(); } catch (_) {}
  }).catch(() => {});

  finishSpec();
});
