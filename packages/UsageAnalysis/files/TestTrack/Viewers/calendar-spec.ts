import {test, expect, chromium} from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/demog.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Calendar viewer', async () => {
  test.setTimeout(600_000);

  // Reuse the existing Chrome session (user is already logged in)
  const browser = await chromium.connectOverCDP('http://127.0.0.1:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find(p => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
  }
  await page.waitForFunction(() => {
    try {
      return typeof grok !== 'undefined'
        && grok.shell
        && typeof grok.shell.closeAll === 'function'
        && grok.dapi
        && grok.dapi.files;
    } catch { return false; }
  }, {timeout: 60000});
  // Wait until Dart bindings are fully registered
  await page.waitForFunction(() => {
    try { grok.shell.closeAll(); return true; }
    catch { return false; }
  }, {timeout: 60000});

  // Phase 2: Open demog
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = false; } catch {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // #### Step 2: Add Calendar via Add Viewer dialog
  await softStep('Add Calendar via Add Viewer dialog', async () => {
    const result = await page.evaluate(async () => {
      // Click Add viewer icon
      const addBtn = document.querySelector('[name="Add-viewer"], [tooltip="Add viewer"]') as HTMLElement
        ?? Array.from(document.querySelectorAll('i, div'))
          .find(el => (el as HTMLElement).getAttribute?.('tooltip') === 'Add viewer') as HTMLElement;
      if (addBtn) addBtn.click();
      else {
        // Fallback: use JS API addViewer
        grok.shell.tv.addViewer('Calendar');
      }
      await new Promise(r => setTimeout(r, 500));
      // If dialog opened, click Calendar item inside it
      const dialog = document.querySelector('.d4-dialog');
      if (dialog) {
        const items = Array.from(dialog.querySelectorAll('*')).filter(
          el => el.textContent?.trim() === 'Calendar' && (el as HTMLElement).offsetParent !== null
        );
        (items[0] as HTMLElement)?.click();
        await new Promise(r => setTimeout(r, 1500));
      }
      return {
        hasCalendar: !!grok.shell.tv.viewers.find((v: any) => v.type === 'Calendar'),
        domHas: !!document.querySelector('[name="viewer-Calendar"]'),
      };
    });
    expect(result.hasCalendar).toBe(true);
    expect(result.domHas).toBe(true);
  });

  // #### Step 3: Close Calendar and reopen via Viewers toolbox icon
  await softStep('Close and reopen via toolbox icon', async () => {
    const result = await page.evaluate(async () => {
      const container = document.querySelector('[name="viewer-Calendar"]');
      const panelBase = container?.closest('.panel-base');
      const closeBtn = panelBase?.querySelector('[name="Close"]') as HTMLElement;
      closeBtn?.click();
      await new Promise(r => setTimeout(r, 500));
      const closed = !document.querySelector('[name="viewer-Calendar"]');
      const icon = document.querySelector('[name="icon-calendar"]') as HTMLElement;
      icon?.click();
      await new Promise(r => setTimeout(r, 1500));
      return {
        closed,
        reopened: !!document.querySelector('[name="viewer-Calendar"]'),
        hasCalendar: !!grok.shell.tv.viewers.find((v: any) => v.type === 'Calendar'),
      };
    });
    expect(result.closed).toBe(true);
    expect(result.reopened).toBe(true);
    expect(result.hasCalendar).toBe(true);
  });

  // #### Step 4a: Hover day cell shows tooltip
  await softStep('Hover day cell shows tooltip', async () => {
    const result = await page.evaluate(async () => {
      const cal = grok.shell.tv.viewers.find((v: any) => v.type === 'Calendar') as any;
      if (!cal.props.dateColumnName) {
        const dc = grok.shell.tv.dataFrame.columns.toList().find((c: any) => c.type === 'datetime');
        cal.props.dateColumnName = dc.name;
      }
      await new Promise(r => setTimeout(r, 500));
      const canvas = document.querySelector('[name="viewer-Calendar"] [name="canvas"]') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const dayX = rect.x + 300, dayY = rect.y + 400;
      canvas.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, cancelable: true, clientX: dayX, clientY: dayY,
      }));
      await new Promise(r => setTimeout(r, 500));
      const tt = document.querySelector('.d4-tooltip, [class*="tooltip"]') as HTMLElement;
      return {tooltip: tt?.innerText ?? ''};
    });
    expect(result.tooltip).toContain('Click to select');
    expect(result.tooltip).toMatch(/\d{4}-\d{2}-\d{2}|\d+ rows/);
  });

  // #### Step 4b: Click day cell selects rows
  await softStep('Click day cell selects rows', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      const canvas = document.querySelector('[name="viewer-Calendar"] [name="canvas"]') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const dayX = rect.x + 300, dayY = rect.y + 400;
      const ev = (type: string, opts: any = {}) => canvas.dispatchEvent(new MouseEvent(type, {
        bubbles: true, cancelable: true, clientX: dayX, clientY: dayY, button: 0, ...opts,
      }));
      ev('mousemove');
      await new Promise(r => setTimeout(r, 300));
      ev('mousedown');
      ev('mouseup');
      ev('click');
      await new Promise(r => setTimeout(r, 300));
      return {selected: df.selection.trueCount};
    });
    expect(result.selected).toBeGreaterThan(0);
  });

  // #### Step 4c: Click month label selects all rows in that month
  await softStep('Click month label selects month rows', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      const canvas = document.querySelector('[name="viewer-Calendar"] [name="canvas"]') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const mX = rect.x + 25, mY = rect.y + 400;
      const ev = (type: string) => canvas.dispatchEvent(new MouseEvent(type, {
        bubbles: true, cancelable: true, clientX: mX, clientY: mY, button: 0,
      }));
      ev('mousemove');
      await new Promise(r => setTimeout(r, 400));
      const tt = (document.querySelector('.d4-tooltip, [class*="tooltip"]') as HTMLElement)?.innerText ?? '';
      ev('mousedown'); ev('mouseup'); ev('click');
      await new Promise(r => setTimeout(r, 300));
      return {tooltip: tt, selected: df.selection.trueCount};
    });
    expect(result.tooltip).toMatch(/\w+ \d{4}/);
    expect(result.selected).toBeGreaterThan(0);
  });

  // #### Step 4d: Click day-of-week header selects all weekday rows
  await softStep('Click weekday header selects weekday rows', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      const canvas = document.querySelector('[name="viewer-Calendar"] [name="canvas"]') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const wX = rect.x + 350, wY = rect.y + 110;
      const ev = (type: string) => canvas.dispatchEvent(new MouseEvent(type, {
        bubbles: true, cancelable: true, clientX: wX, clientY: wY, button: 0,
      }));
      ev('mousemove');
      await new Promise(r => setTimeout(r, 400));
      const tt = (document.querySelector('.d4-tooltip, [class*="tooltip"]') as HTMLElement)?.innerText ?? '';
      ev('mousedown'); ev('mouseup'); ev('click');
      await new Promise(r => setTimeout(r, 300));
      return {tooltip: tt, selected: df.selection.trueCount};
    });
    expect(result.tooltip).toMatch(/Monday|Tuesday|Wednesday|Thursday|Friday|Saturday|Sunday/);
    expect(result.selected).toBeGreaterThan(0);
  });

  // #### Step 4e/f: Shift+click extends, Ctrl+click toggles
  await softStep('Shift+click extends, Ctrl+click toggles selection', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const canvas = document.querySelector('[name="viewer-Calendar"] [name="canvas"]') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      // Start: select a weekday
      df.selection.setAll(false);
      const wX = rect.x + 350, wY = rect.y + 110;
      const evW = (type: string, opts: any = {}) => canvas.dispatchEvent(new MouseEvent(type, {
        bubbles: true, cancelable: true, clientX: wX, clientY: wY, button: 0, ...opts,
      }));
      evW('mousemove'); evW('mousedown'); evW('mouseup'); evW('click');
      await new Promise(r => setTimeout(r, 300));
      const base = df.selection.trueCount;

      // Shift-click a day cell (different region)
      const dayX = rect.x + 450, dayY = rect.y + 500;
      const evD = (type: string, opts: any = {}) => canvas.dispatchEvent(new MouseEvent(type, {
        bubbles: true, cancelable: true, clientX: dayX, clientY: dayY, button: 0, ...opts,
      }));
      evD('mousemove');
      await new Promise(r => setTimeout(r, 300));
      evD('mousedown', {shiftKey: true});
      evD('mouseup', {shiftKey: true});
      evD('click', {shiftKey: true});
      await new Promise(r => setTimeout(r, 300));
      const afterShift = df.selection.trueCount;

      // Ctrl-click same cell to toggle off
      evD('mousedown', {ctrlKey: true});
      evD('mouseup', {ctrlKey: true});
      evD('click', {ctrlKey: true});
      await new Promise(r => setTimeout(r, 300));
      const afterCtrl = df.selection.trueCount;

      return {base, afterShift, afterCtrl};
    });
    expect(result.base).toBeGreaterThan(0);
    expect(result.afterShift).toBeGreaterThan(result.base);
    expect(result.afterCtrl).toBe(result.base);
  });

  // #### Step 5: Gear icon opens Property Pane
  await softStep('Gear icon opens Property Pane', async () => {
    const result = await page.evaluate(async () => {
      const container = document.querySelector('[name="viewer-Calendar"]');
      const panelBase = container?.closest('.panel-base');
      const gear = panelBase?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      gear?.click();
      await new Promise(r => setTimeout(r, 500));
      const pp = document.querySelector('.grok-prop-panel, .d4-property-grid');
      const hasDateProp = !!pp?.querySelector('[name="prop-date"]');
      const hasOnClick = !!pp?.querySelector('[name="prop-on-click"]');
      return {hasPP: !!pp, hasDateProp, hasOnClick};
    });
    expect(result.hasPP).toBe(true);
    expect(result.hasDateProp).toBe(true);
    expect(result.hasOnClick).toBe(true);
  });

  // #### Step 6: Modify properties
  await softStep('Modify properties', async () => {
    const result = await page.evaluate(async () => {
      const cal = grok.shell.tv.viewers.find((v: any) => v.type === 'Calendar') as any;
      const errors: string[] = [];
      const onErr = (e: ErrorEvent) => errors.push(e.message);
      window.addEventListener('error', onErr);

      const before = {
        showHeader: cal.props.showHeader,
        redWeekends: cal.props.redWeekends,
        showFilteredOnly: cal.props.showFilteredOnly,
        onClick: cal.props.onClick,
      };

      cal.props.showHeader = false;
      await new Promise(r => setTimeout(r, 200));
      cal.props.redWeekends = false;
      await new Promise(r => setTimeout(r, 200));
      cal.props.showFilteredOnly = false;
      await new Promise(r => setTimeout(r, 200));
      cal.props.onClick = 'Filter';
      await new Promise(r => setTimeout(r, 300));

      const after = {
        showHeader: cal.props.showHeader,
        redWeekends: cal.props.redWeekends,
        showFilteredOnly: cal.props.showFilteredOnly,
        onClick: cal.props.onClick,
      };

      // reset
      cal.props.showHeader = true;
      cal.props.redWeekends = true;
      cal.props.showFilteredOnly = true;
      cal.props.onClick = 'Select';

      window.removeEventListener('error', onErr);
      return {before, after, errors};
    });
    expect(result.before.showHeader).toBe(true);
    expect(result.after.showHeader).toBe(false);
    expect(result.after.redWeekends).toBe(false);
    expect(result.after.showFilteredOnly).toBe(false);
    expect(result.after.onClick).toBe('Filter');
    expect(result.errors).toEqual([]);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
