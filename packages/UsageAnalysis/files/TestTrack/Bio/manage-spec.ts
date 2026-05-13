import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Bio Manage Monomer Libraries', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:AppData/Bio/samples/HELM.csv');
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 3000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const hasBioChem = cols.some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Open HELM.csv and detect Macromolecule semType', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      return {
        rows: df.rowCount,
        helm: cols.find((c: any) => c.name === 'HELM')?.semType,
      };
    });
    expect(info.rows).toBe(540);
    expect(info.helm).toBe('Macromolecule');
  });

  await softStep('Open Bio > Manage > Monomer Libraries', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      const manage = document.querySelector('[name="div-Bio---Manage"]')!;
      manage.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      manage.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 500));
      (document.querySelector('[name="div-Bio---Manage---Monomer-Libraries"]') as HTMLElement).click();
    });
    await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Manage Monomer Libraries',
      null, {timeout: 30000});
    const view = await page.evaluate(() => ({name: grok.shell.v?.name, type: grok.shell.v?.type}));
    expect(view.name).toBe('Manage Monomer Libraries');
    expect(view.type).toBe('view');
  });

  await softStep('Library list shows 5 checkboxes, all checked by default', async () => {
    const state = await page.evaluate(() => {
      const root = grok.shell.v.root;
      const checkboxes = Array.from(root.querySelectorAll('input[type="checkbox"]')) as HTMLInputElement[];
      const rows = checkboxes.map((cb) => {
        let el = cb.parentElement;
        let label = '';
        for (let i = 0; i < 5 && el; i++) {
          const t = (el.textContent || '').trim();
          if (t && t.length < 200) { label = t; break; }
          el = el.parentElement;
        }
        return {checked: cb.checked, label};
      });
      return {count: checkboxes.length, rows, allChecked: checkboxes.every((c) => c.checked)};
    });
    expect(state.count).toBe(5);
    expect(state.allChecked).toBe(true);
    expect(state.rows.some((r) => /HELMCoreLibrary\.json/i.test(r.label))).toBe(true);
    expect(state.rows.some((r) => /polytool-lib\.json/i.test(r.label))).toBe(true);
  });

  await softStep('Toggling a library saves settings and updates duplicates panel', async () => {
    const result = await page.evaluate(async () => {
      const root = grok.shell.v.root;
      const before = root.querySelectorAll('.d4-card, .ui-card, [class*="duplicate"]').length;
      const beforeChildren = root.querySelectorAll('*').length;
      const checkboxes = Array.from(root.querySelectorAll('input[type="checkbox"]')) as HTMLInputElement[];
      const idx = 1;
      const target = checkboxes[idx];
      const wasChecked = target.checked;
      target.click();
      await new Promise((r) => setTimeout(r, 3000));
      const afterToggle = root.querySelectorAll('.d4-card, .ui-card, [class*="duplicate"]').length;
      const afterToggleChildren = root.querySelectorAll('*').length;
      const balloons = Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon, .grok-notification, [class*="notification"]'))
        .map((n) => (n.textContent || '').trim()).filter(Boolean);
      const toggledOff = !(Array.from(root.querySelectorAll('input[type="checkbox"]')) as HTMLInputElement[])[idx].checked;
      // toggle back
      (Array.from(root.querySelectorAll('input[type="checkbox"]')) as HTMLInputElement[])[idx].click();
      await new Promise((r) => setTimeout(r, 2500));
      const restored = (Array.from(root.querySelectorAll('input[type="checkbox"]')) as HTMLInputElement[])[idx].checked;
      const afterRestore = root.querySelectorAll('.d4-card, .ui-card, [class*="duplicate"]').length;
      return {wasChecked, toggledOff, restored, before, afterToggle, afterRestore, beforeChildren, afterToggleChildren, balloons};
    });
    expect(result.wasChecked).toBe(true);
    expect(result.toggledOff).toBe(true);
    expect(result.restored).toBe(true);
    expect(result.afterToggle).toBeLessThan(result.before);
    expect(result.afterRestore).toBe(result.before);
    expect(result.balloons.some((b) => b.includes('Monomer library user settings saved'))).toBe(true);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
