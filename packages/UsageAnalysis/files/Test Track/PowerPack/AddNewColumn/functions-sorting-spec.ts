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

test('Add New Column: functions sorting', async ({page}) => {
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

  // Setup: open SPGI.csv, wait for semType and cell rendering
  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});

  await softStep('Step 2: Open Add New Column dialog', async () => {
    await page.evaluate(async () => {
      const editMenu = document.querySelector('[name="div-Edit"]') as HTMLElement;
      editMenu.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      editMenu.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
      editMenu.dispatchEvent(new MouseEvent('mouseup', {bubbles: true}));
      editMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 300));
      const addCol = document.querySelector('[name="div-Edit---Add-New-Column..."]') as HTMLElement;
      addCol.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('[name="dialog-Add-New-Column"]').waitFor({timeout: 5000});
    expect(await page.locator('[name="dialog-Add-New-Column"]').count()).toBe(1);
  });

  await softStep('Step 3a: Click Structure column — molecule-input functions on top', async () => {
    const first15 = await page.evaluate(async () => {
      const gridEl = document.querySelector('[name="dialog-Add-New-Column"] .add-new-column-columns-grid')!;
      const canvases = gridEl.querySelectorAll('canvas');
      const mainCanvas = canvases[canvases.length - 1] as HTMLCanvasElement;
      const r = mainCanvas.getBoundingClientRect();
      const rowHeight = r.height / 21;
      const xClick = r.x + r.width / 2;
      const yClick = r.y + rowHeight * 1 + rowHeight / 2; // row 1 = Structure
      const opts = (clientX: number, clientY: number): MouseEventInit => ({
        bubbles: true, cancelable: true, view: window, clientX, clientY, button: 0, buttons: 1,
      });
      mainCanvas.dispatchEvent(new MouseEvent('mousedown', opts(xClick, yClick)));
      mainCanvas.dispatchEvent(new MouseEvent('mouseup', opts(xClick, yClick)));
      mainCanvas.dispatchEvent(new MouseEvent('click', opts(xClick, yClick)));
      await new Promise((r) => setTimeout(r, 800));
      return Array.from(document.querySelectorAll('[name="dialog-Add-New-Column"] #actions span[name^="span-"]'))
        .slice(0, 15)
        .map((s) => {
          const name = s.getAttribute('name')!.replace(/^span-/, '');
          const params = s.parentElement!.querySelector('.grok-function-params')?.textContent?.trim() ?? '';
          return {name, params};
        });
    });
    const molTerms = /molecule|smiles|mol[A-Z]|molCol/;
    const molCount = first15.slice(0, 5).filter((f) => molTerms.test(f.params) || molTerms.test(f.name)).length;
    expect(molCount).toBeGreaterThanOrEqual(3);
  });

  await softStep('Step 3b: Click Chemical Space X column — numeric-input functions on top', async () => {
    const first15 = await page.evaluate(async () => {
      const gridEl = document.querySelector('[name="dialog-Add-New-Column"] .add-new-column-columns-grid')!;
      const canvases = gridEl.querySelectorAll('canvas');
      const mainCanvas = canvases[canvases.length - 1] as HTMLCanvasElement;
      const r = mainCanvas.getBoundingClientRect();
      const rowHeight = r.height / 21;
      const xClick = r.x + r.width / 2;
      const yClick = r.y + rowHeight * 18 + rowHeight / 2; // row 18 = Chemical Space X
      const opts = (clientX: number, clientY: number): MouseEventInit => ({
        bubbles: true, cancelable: true, view: window, clientX, clientY, button: 0, buttons: 1,
      });
      mainCanvas.dispatchEvent(new MouseEvent('mousedown', opts(xClick, yClick)));
      mainCanvas.dispatchEvent(new MouseEvent('mouseup', opts(xClick, yClick)));
      mainCanvas.dispatchEvent(new MouseEvent('click', opts(xClick, yClick)));
      await new Promise((r) => setTimeout(r, 800));
      return Array.from(document.querySelectorAll('[name="dialog-Add-New-Column"] #actions span[name^="span-"]'))
        .slice(0, 15)
        .map((s) => {
          const name = s.getAttribute('name')!.replace(/^span-/, '');
          const params = s.parentElement!.querySelector('.grok-function-params')?.textContent?.trim() ?? '';
          return {name, params};
        });
    });
    const numericReturn = /: (num|double|int|qnum)\b/;
    const numCount = first15.slice(0, 5).filter((f) => numericReturn.test(f.params)).length;
    expect(numCount).toBeGreaterThanOrEqual(3);
  });

  await softStep('Step 4: Click sort icon → By name — functions alphabetical', async () => {
    await page.evaluate(async () => {
      const sortIcon = document.querySelector('[name="dialog-Add-New-Column"] [name="icon-sort-alt"]') as HTMLElement;
      sortIcon.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
      sortIcon.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, button: 0}));
      sortIcon.dispatchEvent(new MouseEvent('click', {bubbles: true, button: 0}));
      await new Promise((r) => setTimeout(r, 500));
      const byName = document.querySelector('[name="div-By-name"]') as HTMLElement;
      byName.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 800));
    });
    const first5 = await page.evaluate(() =>
      Array.from(document.querySelectorAll('[name="dialog-Add-New-Column"] #actions span[name^="span-"]'))
        .slice(0, 5)
        .map((s) => s.getAttribute('name')!.replace(/^span-/, '')),
    );
    const sortedCopy = [...first5].sort((a, b) => a.localeCompare(b, undefined, {sensitivity: 'base'}));
    expect(first5).toEqual(sortedCopy);
  });

  await softStep('Step 5: Click different columns — function order does not change', async () => {
    const snapshots = await page.evaluate(async () => {
      const gridEl = document.querySelector('[name="dialog-Add-New-Column"] .add-new-column-columns-grid')!;
      const canvases = gridEl.querySelectorAll('canvas');
      const mainCanvas = canvases[canvases.length - 1] as HTMLCanvasElement;
      const r = mainCanvas.getBoundingClientRect();
      const rowHeight = r.height / 21;
      const xClick = r.x + r.width / 2;
      const opts = (clientX: number, clientY: number): MouseEventInit => ({
        bubbles: true, cancelable: true, view: window, clientX, clientY, button: 0, buttons: 1,
      });
      const clickRow = async (idx: number) => {
        const y = r.y + rowHeight * idx + rowHeight / 2;
        mainCanvas.dispatchEvent(new MouseEvent('mousedown', opts(xClick, y)));
        mainCanvas.dispatchEvent(new MouseEvent('mouseup', opts(xClick, y)));
        mainCanvas.dispatchEvent(new MouseEvent('click', opts(xClick, y)));
        await new Promise((r) => setTimeout(r, 600));
      };
      const snap = () =>
        Array.from(document.querySelectorAll('[name="dialog-Add-New-Column"] #actions span[name^="span-"]'))
          .slice(0, 15)
          .map((s) => s.getAttribute('name')!.replace(/^span-/, ''));
      const out: string[][] = [];
      await clickRow(1); out.push(snap());  // Structure (Molecule)
      await clickRow(18); out.push(snap()); // Chemical Space X (double)
      await clickRow(4); out.push(snap());  // Chemist (string)
      return out;
    });
    expect(snapshots[1]).toEqual(snapshots[0]);
    expect(snapshots[2]).toEqual(snapshots[0]);
  });

  await page.evaluate(() => {
    const cancel = document.querySelector('[name="dialog-Add-New-Column"] [name="button-CANCEL"]') as HTMLElement | null;
    cancel?.click();
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
