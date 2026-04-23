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

test('Word cloud', async ({page}) => {
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

  await page.evaluate(async () => {
    const w = window as any;
    document.body.classList.add('selenium');
    w.grok.shell.settings.showFiltersIconsConstantly = true;
    w.grok.shell.windows.simpleMode = true;
    w.grok.shell.closeAll();
    const df = await w.grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
    const tv = w.grok.shell.addTableView(df);
    await new Promise((resolve: any) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const cols: any[] = [];
    for (let i = 0; i < df.columns.length; i++) cols.push(df.columns.byIndex(i));
    const hasBioChem = cols.some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r: any) => setTimeout(r, 200));
      }
      await new Promise((r: any) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Step 2: Open Add Viewer gallery, pick Word Cloud', async () => {
    await page.locator('i[aria-label="Add viewer"]').click();
    await page.waitForTimeout(800);
    // Click Word Cloud tile in gallery
    await page.locator('.d4-dialog, .ui-dialog').locator('text="Word Cloud"').first().click();
    await page.waitForTimeout(1500);
    const added = await page.evaluate(() => {
      const w = window as any;
      return !!w.grok.shell.tv.viewers.find((v: any) => v.type === 'Word cloud');
    });
    expect(added).toBe(true);
  });

  await softStep('Step 3: Close, re-add via Toolbox Viewers > Word Cloud icon', async () => {
    await page.evaluate(async () => {
      const w = window as any;
      const root = document.querySelector('[name="viewer-Word-cloud"]')!;
      (root.closest('.panel-base')!.querySelector('[name="Close"]') as HTMLElement).click();
      await new Promise((r: any) => setTimeout(r, 500));
    });
    await page.locator('[name="icon-Word-cloud"]').click();
    await page.waitForTimeout(1500);
    const added = await page.evaluate(() => {
      const w = window as any;
      return !!w.grok.shell.tv.viewers.find((v: any) => v.type === 'Word cloud');
    });
    expect(added).toBe(true);
  });

  await softStep('Step 4: Hamburger menu opens (interaction with displayed text uses canvas — not DOM-testable)', async () => {
    const popupShown = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Word-cloud"]')!;
      const panel = root.closest('.panel-base')!;
      (panel.querySelector('[name="icon-font-icon-menu"]') as HTMLElement).click();
      await new Promise((r: any) => setTimeout(r, 600));
      const popup = Array.from(document.querySelectorAll('.d4-menu-popup'))
        .find(e => !!(e as HTMLElement).offsetParent);
      const text = popup ? (popup as HTMLElement).innerText : '';
      document.body.click();
      return text;
    });
    expect(popupShown).toContain('Properties');
  });

  await softStep('Step 5: Open Property Pane via Gear icon', async () => {
    await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Word-cloud"]')!;
      const panel = root.closest('.panel-base')!;
      (panel.querySelector('[name="icon-font-icon-settings"]') as HTMLElement).click();
    });
    await page.waitForTimeout(800);
    const hasData = await page.evaluate(() => {
      return !!Array.from(document.querySelectorAll('.property-grid-category'))
        .find(e => e.textContent?.trim() === 'Data');
    });
    expect(hasData).toBe(true);
  });

  await softStep('Step 6: Modify column, minTextSize, maxTextSize, bold, rotationStep', async () => {
    const after = await page.evaluate(async () => {
      const w = window as any;
      const wc = w.grok.shell.tv.viewers.find((v: any) => v.type === 'Word cloud');
      wc.props.columnColumnName = 'Stereo Category';
      await new Promise((r: any) => setTimeout(r, 400));
      wc.props.minTextSize = 20;
      wc.props.maxTextSize = 80;
      wc.props.bold = true;
      wc.props.rotationStep = 90;
      await new Promise((r: any) => setTimeout(r, 600));
      return {
        column: wc.props.columnColumnName,
        minTextSize: wc.props.minTextSize,
        maxTextSize: wc.props.maxTextSize,
        bold: wc.props.bold,
        rotationStep: wc.props.rotationStep,
      };
    });
    expect(after.column).toBe('Stereo Category');
    expect(after.minTextSize).toBe(20);
    expect(after.maxTextSize).toBe(80);
    expect(after.rotationStep).toBe(90);
  });

  await softStep('Close viewer cleanly', async () => {
    const gone = await page.evaluate(async () => {
      const w = window as any;
      const root = document.querySelector('[name="viewer-Word-cloud"]')!;
      (root.closest('.panel-base')!.querySelector('[name="Close"]') as HTMLElement).click();
      await new Promise((r: any) => setTimeout(r, 500));
      return !w.grok.shell.tv.viewers.find((v: any) => v.type === 'Word cloud');
    });
    expect(gone).toBe(true);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n'));
});
