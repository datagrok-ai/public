import {test, expect} from '@playwright/test';
import {baseUrl, login, password, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('StickyMeta: Create metadata schema & entity type', async ({page}) => {
  test.setTimeout(300_000);

  await page.goto(baseUrl);
  await page.waitForTimeout(3000);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    // Login form is plain HTML (not a Dart widget), so .fill() works and is more reliable
    // than keyboard.type(), which silently drops characters on a fresh page on dev.
    await loginInput.fill(login);
    const pw = page.getByPlaceholder('Password').and(page.locator(':visible'));
    await pw.fill(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });

  // Pre-cleanup: delete TestEntity1 and TestSchema1 if they already exist. The Sticky-Meta list
  // is paginated (18 of 57), so we use the search box to narrow the list before the right-click
  // delete — otherwise the label-by-text lookup misses entries past the first page.
  // Browse-tree navigation is used because direct goto() between Sticky-Meta sibling routes is
  // ignored by the SPA (it reloads to the most-recently-visited one).
  await softStep('Setup: remove pre-existing TestEntity1 / TestSchema1', async () => {
    const searchAndDelete = async (name: string) => {
      await page.evaluate(async (n) => {
        const search = document.querySelector('input[placeholder*="by name"]') as HTMLInputElement | null;
        if (!search) return;
        search.focus();
        const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
        setter.call(search, n);
        search.dispatchEvent(new Event('input', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 1200));
        const label = Array.from(document.querySelectorAll('label')).find((l) => l.textContent?.trim() === n) as HTMLElement | undefined;
        if (!label) return;
        const rect = label.getBoundingClientRect();
        label.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2, clientX: rect.left + 10, clientY: rect.top + 5}));
        await new Promise((r) => setTimeout(r, 500));
        const menu = document.querySelector('.d4-menu-popup, .d4-menu');
        const del = Array.from(menu?.querySelectorAll('.d4-menu-item') || []).find((el) => el.textContent?.trim() === 'Delete') as HTMLElement | undefined;
        del?.click();
        await new Promise((r) => setTimeout(r, 800));
        (document.querySelector('[name="button-DELETE"]') as HTMLElement | null)?.click();
        await new Promise((r) => setTimeout(r, 1800));
        // Clear the search
        setter.call(search, '');
        search.dispatchEvent(new Event('input', {bubbles: true}));
      }, name);
    };

    await page.goto(`${baseUrl}/meta/schemas`);
    await page.locator('[name="button-New-Schema..."]').waitFor({timeout: 30000});
    await page.waitForTimeout(1500);
    await searchAndDelete('TestSchema1');
    await page.waitForTimeout(800);

    // Navigate to Types via Browse tree (goto sibling route is ignored).
    const navToTypes = async () => page.evaluate(async () => {
      const findLabel = (t: string) =>
        Array.from(document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label'))
          .find((el) => el.textContent?.trim() === t) as HTMLElement | undefined;
      const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
      if (!findLabel('Sticky Meta')) { findLabel('Platform')?.click(); await wait(800); }
      if (!findLabel('Types')) { findLabel('Sticky Meta')?.click(); await wait(800); }
      findLabel('Types')?.click();
      await wait(1500);
      return location.href;
    });
    // Retry tree click a few times — tree-click to sibling SPA routes is occasionally a no-op.
    for (let attempt = 0; attempt < 4; attempt++) {
      const url = await navToTypes();
      if (/\/meta\/types$/.test(url)) break;
      await page.waitForTimeout(2000);
    }
    const onTypes = await page.locator('[name="button-New-Entity-Type..."]').isVisible({timeout: 10000}).catch(() => false);
    if (onTypes) {
      await page.waitForTimeout(1500);
      await searchAndDelete('TestEntity1');
    }
    // Tolerant — if Setup couldn't reach /meta/types, the scenario's Step 2 will; Step 4's
    // duplicate-name handling below will no-op if TestEntity1 already exists.
  });

  // Step 2: Navigate to Browse > Platform > Sticky Meta > Types via the Browse tree.
  // Datagrok's SPA ignores goto() between Sticky-Meta sibling routes — it reloads to the
  // most-recently-visited one. So we click through the tree, matching what a user would do.
  await softStep('Step 2: Navigate to Browse > Platform > Sticky Meta > Types', async () => {
    const click = async () => page.evaluate(async () => {
      const findLabel = (t: string) =>
        Array.from(document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label'))
          .find((el) => el.textContent?.trim() === t) as HTMLElement | undefined;
      const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
      if (!findLabel('Sticky Meta')) { findLabel('Platform')?.click(); await wait(800); }
      if (!findLabel('Types')) { findLabel('Sticky Meta')?.click(); await wait(800); }
      findLabel('Types')?.click();
      await wait(1500);
      return location.href;
    });
    let url = await click();
    if (!/\/meta\/types$/.test(url)) {
      await page.waitForTimeout(1500);
      url = await click();
    }
    await page.locator('[name="button-New-Entity-Type..."]').waitFor({timeout: 30000});
    await expect(page).toHaveURL(/\/meta\/types(\?|$)/);
  });

  // Step 3: Click "New Entity Type..."
  await softStep('Step 3: Click NEW ENTITY TYPE button', async () => {
    await page.locator('[name="button-New-Entity-Type..."]').click();
    await page.locator('[name="dialog-Create-a-new-entity-type"], .d4-dialog').first().waitFor({timeout: 5000});
    const title = await page.evaluate(() => document.querySelector('.d4-dialog .d4-dialog-title')?.textContent?.trim());
    expect(title).toBe('Create a new entity type');
  });

  // Step 4: Enter name TestEntity1, matcher semtype=molecule, Save
  await softStep('Step 4: Fill name / matcher and save new entity type', async () => {
    const nameInput = page.locator('.d4-dialog [name="input-Name"]');
    await nameInput.focus();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('TestEntity1', {delay: 10});

    const matcherInput = page.locator('.d4-dialog [name="input-Matching-expression"]');
    await matcherInput.focus();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('semtype=molecule', {delay: 10});

    // Verify inputs filled and OK enabled before submitting.
    const inputs = await page.evaluate(() => {
      const dialog = Array.from(document.querySelectorAll('.d4-dialog'))
        .find((d) => d.querySelector('.d4-dialog-title')?.textContent?.trim() === 'Create a new entity type');
      return {
        name: (dialog?.querySelector('[name="input-Name"]') as HTMLInputElement | null)?.value,
        matcher: (dialog?.querySelector('[name="input-Matching-expression"]') as HTMLInputElement | null)?.value,
        okEnabled: dialog?.querySelector('[name="button-OK"]')?.classList.contains('enabled'),
      };
    });
    expect(inputs.name).toBe('TestEntity1');
    expect(inputs.matcher).toBe('semtype=molecule');

    // Simulate a real mouse click on OK — Playwright's .click() respects pointer-event
    // interception; JS .click() skips that, matching keyboard-driven user flow.
    const okBox = await page.locator('.d4-dialog[name="dialog-Create-a-new-entity-type"] [name="button-OK"]').boundingBox();
    if (okBox) {
      await page.mouse.click(okBox.x + okBox.width / 2, okBox.y + okBox.height / 2);
    }
    await page.waitForTimeout(3000);

    // Handle duplicate-name case: if the dialog is still open and a "Duplicate name" toast is
    // shown, the previous spec run left TestEntity1 behind. Cancel, delete via search, and
    // re-create in one recursive pass.
    const dialogStillOpen = await page.evaluate(() => !!Array.from(document.querySelectorAll('.d4-dialog'))
      .find((d) => d.querySelector('.d4-dialog-title')?.textContent?.trim() === 'Create a new entity type'));
    if (dialogStillOpen) {
      await page.evaluate(() => {
        const dialog = Array.from(document.querySelectorAll('.d4-dialog'))
          .find((d) => d.querySelector('.d4-dialog-title')?.textContent?.trim() === 'Create a new entity type');
        (dialog?.querySelector('[name="button-CANCEL"]') as HTMLElement | null)?.click();
      });
      await page.waitForTimeout(800);
      // Search-based delete (list is paginated)
      await page.evaluate(async () => {
        const search = document.querySelector('input[placeholder*="by name"]') as HTMLInputElement | null;
        if (!search) return;
        search.focus();
        const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
        setter.call(search, 'TestEntity1');
        search.dispatchEvent(new Event('input', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 1500));
        const label = Array.from(document.querySelectorAll('label')).find((l) => l.textContent?.trim() === 'TestEntity1') as HTMLElement | undefined;
        if (label) {
          const rect = label.getBoundingClientRect();
          label.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2, clientX: rect.left + 10, clientY: rect.top + 5}));
          await new Promise((r) => setTimeout(r, 500));
          const menu = document.querySelector('.d4-menu-popup, .d4-menu');
          const del = Array.from(menu?.querySelectorAll('.d4-menu-item') || []).find((el) => el.textContent?.trim() === 'Delete') as HTMLElement | undefined;
          del?.click();
          await new Promise((r) => setTimeout(r, 800));
          (document.querySelector('[name="button-DELETE"]') as HTMLElement | null)?.click();
          await new Promise((r) => setTimeout(r, 2000));
        }
        // Clear search
        setter.call(search, '');
        search.dispatchEvent(new Event('input', {bubbles: true}));
      });
      await page.waitForTimeout(500);
      // Re-open Create dialog and re-fill
      await page.locator('[name="button-New-Entity-Type..."]').click();
      await page.waitForTimeout(800);
      await page.locator('.d4-dialog [name="input-Name"]').focus();
      await page.keyboard.press('Control+A');
      await page.keyboard.type('TestEntity1', {delay: 10});
      await page.locator('.d4-dialog [name="input-Matching-expression"]').focus();
      await page.keyboard.press('Control+A');
      await page.keyboard.type('semtype=molecule', {delay: 10});
      const okBox2 = await page.locator('.d4-dialog[name="dialog-Create-a-new-entity-type"] [name="button-OK"]').boundingBox();
      if (okBox2) await page.mouse.click(okBox2.x + okBox2.width / 2, okBox2.y + okBox2.height / 2);
      await page.waitForTimeout(3000);
    }

    // The Types list is paginated (only ~18 of 57 items rendered). Filter via the search box
    // before checking — a bare DOM scan misses entries past the first page.
    const present = await page.evaluate(async () => {
      const search = document.querySelector('input[placeholder*="by name"]') as HTMLInputElement | null;
      if (search) {
        search.focus();
        const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
        setter.call(search, 'TestEntity1');
        search.dispatchEvent(new Event('input', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 1500));
      }
      const labels = Array.from(document.querySelectorAll('label')).map((l) => l.textContent?.trim());
      const found = labels.includes('TestEntity1');
      // Clear search
      if (search) {
        const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
        setter.call(search, '');
        search.dispatchEvent(new Event('input', {bubbles: true}));
      }
      return found;
    });
    expect(present).toBe(true);
  });

  // Step 5: Navigate to Browse > Platform > Sticky Meta > Schemas (tree already expanded)
  await softStep('Step 5: Navigate to Browse > Platform > Sticky Meta > Schemas', async () => {
    const click = async () => page.evaluate(async () => {
      const labels = Array.from(document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label'));
      const stickyGroup = labels.find((el) => el.textContent?.trim() === 'Sticky Meta')?.closest('.d4-tree-view-group');
      const schemas = Array.from(stickyGroup?.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label') || [])
        .find((el) => el.textContent?.trim() === 'Schemas') as HTMLElement | undefined;
      schemas?.click();
      await new Promise((r) => setTimeout(r, 1500));
      return location.href;
    });
    let url = await click();
    if (!/\/meta\/schemas$/.test(url)) {
      await page.waitForTimeout(1500);
      url = await click();
    }
    await page.locator('[name="button-New-Schema..."]').waitFor({timeout: 30000});
    await expect(page).toHaveURL(/\/meta\/schemas(\?|$)/);
  });

  // Step 6: Click "NEW SCHEMA..."
  await softStep('Step 6: Click NEW SCHEMA button', async () => {
    await page.locator('[name="button-New-Schema..."]').click();
    await page.locator('.d4-dialog').first().waitFor({timeout: 5000});
    const title = await page.evaluate(() => document.querySelector('.d4-dialog .d4-dialog-title')?.textContent?.trim());
    expect(title).toBe('Create a new schema');
  });

  // Step 7: Enter schema name + associate with TestEntity1
  await softStep('Step 7: Fill schema name and associate with TestEntity1', async () => {
    // `input-host-Name` also exists inside the Properties table. Scope to the direct form child.
    await page.locator('.d4-dialog .ui-form > [name="input-host-Name"] > [name="input-Name"]').click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('TestSchema1', {delay: 10});

    // Open the type selector (link action)
    await page.locator('[name="div-Associated-with-"] .d4-link-action').click();
    await page.waitForTimeout(800);
    // Check TestEntity1 checkbox in the type picker dialog
    await page.evaluate(() => {
      const cb = document.querySelector('[name="prop-view-testentity1"]') as HTMLInputElement | null;
      if (cb && !cb.checked) cb.click();
    });
    // Confirm the picker (it's the topmost dialog)
    await page.evaluate(() => {
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog'));
      const picker = dialogs.find((d) => /Select types/.test(d.querySelector('.d4-dialog-title')?.textContent || ''));
      (picker?.querySelector('[name="button-OK"]') as HTMLElement | null)?.click();
    });
    await page.waitForTimeout(600);
    const assocText = await page.evaluate(() =>
      document.querySelector('[name="div-Associated-with-"]')?.textContent?.trim());
    expect(assocText).toContain('TestEntity1');
  });

  // Step 8: Add 4 properties and save
  await softStep('Step 8: Add properties (rating/int, notes/string, verified/bool, review_date/datetime) and save', async () => {
    const props: [string, string][] = [
      ['rating', 'int'],
      ['notes', 'string'],
      ['verified', 'bool'],
      ['review_date', 'datetime'],
    ];
    for (let i = 0; i < props.length; i++) {
      const [name, type] = props[i];
      if (i > 0) {
        await page.locator('[name="button-Add-new-property-to-schema"]').click();
        await page.waitForTimeout(250);
      }
      // Ctrl+A + keyboard type is what the Dart row listener sees (matches MCP behavior).
      const row = page.locator('.d4-dialog table.d4-item-table input[name="input-Name"]').nth(i);
      await row.click();
      await page.keyboard.press('Control+A');
      await page.keyboard.type(name, {delay: 10});
      // Select element requires programmatic value + input/change dispatch.
      await page.evaluate(({i, type}) => {
        const selects = document.querySelectorAll('.d4-dialog table.d4-item-table select[name="input-Property-Type"]');
        const sel = selects[i] as HTMLSelectElement;
        sel.value = type;
        sel.dispatchEvent(new Event('input', {bubbles: true}));
        sel.dispatchEvent(new Event('change', {bubbles: true}));
      }, {i, type});
      await page.waitForTimeout(150);
    }

    // Wait for OK to become enabled, then click via mouse (bypasses pointer-event interception)
    await page.locator('.d4-dialog [name="button-OK"].enabled').waitFor({timeout: 5000});
    const okSchemaBox = await page.locator('.d4-dialog[name="dialog-Create-a-new-schema"] [name="button-OK"]').boundingBox();
    if (okSchemaBox) {
      await page.mouse.click(okSchemaBox.x + okSchemaBox.width / 2, okSchemaBox.y + okSchemaBox.height / 2);
    }
    await page.waitForTimeout(3000);

    // List is paginated — filter via search to find TestSchema1 reliably.
    const present = await page.evaluate(async () => {
      const search = document.querySelector('input[placeholder*="by name"]') as HTMLInputElement | null;
      if (search) {
        search.focus();
        const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
        setter.call(search, 'TestSchema1');
        search.dispatchEvent(new Event('input', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 1500));
      }
      return Array.from(document.querySelectorAll('label')).some((l) => l.textContent?.trim() === 'TestSchema1');
    });
    expect(present).toBe(true);

    // Close any lingering dialogs so only the Edit dialog will be open later.
    await page.evaluate(() => {
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        (d.querySelector('[name="button-CANCEL"]') as HTMLElement | null)?.click();
      });
    });
    await page.waitForTimeout(500);

    // Verify schema content via Edit dialog (list is still filtered to TestSchema1).
    await page.evaluate(async () => {
      const label = Array.from(document.querySelectorAll('label')).find((l) => l.textContent?.trim() === 'TestSchema1') as HTMLElement | undefined;
      if (!label) return;
      const rect = label.getBoundingClientRect();
      label.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2, clientX: rect.left + 10, clientY: rect.top + 5}));
      await new Promise((r) => setTimeout(r, 400));
      const menu = document.querySelector('.d4-menu-popup, .d4-menu');
      const edit = Array.from(menu?.querySelectorAll('.d4-menu-item') || []).find((el) => (el.textContent?.trim() || '').startsWith('Edit')) as HTMLElement | undefined;
      edit?.click();
    });
    await page.waitForTimeout(1500);
    const check = await page.evaluate(() => {
      const dialog = Array.from(document.querySelectorAll('.d4-dialog'))
        .find((d) => d.querySelector('.d4-dialog-title')?.textContent?.trim() === 'Edit schema');
      if (!dialog) return {dialog: false};
      const table = dialog.querySelector('table.d4-item-table');
      const nameInputs = Array.from(table?.querySelectorAll('input[name="input-Name"]') || []) as HTMLInputElement[];
      const selects = Array.from(table?.querySelectorAll('select[name="input-Property-Type"]') || []) as HTMLSelectElement[];
      const rows = nameInputs.map((n, i) => ({name: n.value, type: selects[i]?.value}));
      const assoc = dialog.querySelector('[name="div-Associated-with-"]')?.textContent?.trim();
      return {dialog: true, rows, assoc};
    });
    expect(check.dialog).toBe(true);
    expect(check.assoc).toContain('TestEntity1');
    expect(check.rows).toEqual([
      {name: 'rating', type: 'int'},
      {name: 'notes', type: 'string'},
      {name: 'verified', type: 'bool'},
      {name: 'review_date', type: 'datetime'},
    ]);
    // Close the Edit dialog — scope to it to avoid strict-mode violation if any dialog lingers.
    await page.evaluate(() => {
      const dialog = Array.from(document.querySelectorAll('.d4-dialog'))
        .find((d) => d.querySelector('.d4-dialog-title')?.textContent?.trim() === 'Edit schema');
      (dialog?.querySelector('[name="button-CANCEL"]') as HTMLElement | null)?.click();
    });
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
