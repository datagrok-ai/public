import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';

const earthquakesPath = 'System:DemoFiles/geo/earthquakes.csv';
const demogPath = 'System:DemoFiles/demog.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e?.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e?.message ?? e}`);
  }
}

test('Radar viewer (Charts package)', async ({page}) => {
  test.setTimeout(300_000);

  // Login
  await page.goto(baseUrl);
  await page.waitForTimeout(2000);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});

  // Baseline environment setup
  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    (window as any).grok.shell.closeAll();
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
  });

  // Step 1: Open earthquakes.csv and add Radar viewer
  await softStep('Step 1: Open earthquakes.csv and add Radar viewer', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      tv.addViewer('Radar');
      await new Promise((r) => setTimeout(r, 2000));
      const viewerTypes: string[] = [];
      for (const v of tv.viewers) viewerTypes.push(v.type);
      return {rowCount: df.rowCount, viewerTypes};
    }, earthquakesPath);
    expect(result.rowCount).toBe(2426);
    expect(result.viewerTypes).toContain('Radar');
  });

  // Step 2: Open demog.csv and add Radar viewer
  await softStep('Step 2: Open demog.csv and add Radar viewer', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      tv.addViewer('Radar');
      await new Promise((r) => setTimeout(r, 2000));
      const viewerTypes: string[] = [];
      for (const v of tv.viewers) viewerTypes.push(v.type);
      return {rowCount: df.rowCount, viewerTypes};
    }, demogPath);
    expect(result.rowCount).toBe(5850);
    expect(result.viewerTypes).toContain('Radar');
  });

  // Step 3: Context Panel properties — check categories present; toggle one Style color
  await softStep('Step 3: Verify property categories and toggle a Style color', async () => {
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let radar: any = null;
      for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; break; }
      if (!radar) return {categories: [] as string[], colorEchoOk: false, newColor: 0, echoed: 0};

      const props = radar.props.getProperties();
      const categories: string[] = [];
      for (const p of props) {
        const c = p.category as string;
        if (c && !categories.includes(c)) categories.push(c);
      }

      // Toggle one Style color property and verify the getter echoes the new value
      const newColor = 0xFF123456 | 0;
      radar.props.set('backgroundMinColor', newColor);
      await new Promise((r) => setTimeout(r, 500));
      const echoed = radar.props.get('backgroundMinColor');
      const colorEchoOk = echoed === newColor;

      return {categories, colorEchoOk, newColor, echoed};
    });
    expect(result.categories).toEqual(expect.arrayContaining(['Data', 'Selection', 'Value', 'Style', 'Legend']));
    expect(result.colorEchoOk).toBe(true);
  });

  // Cleanup
  await page.evaluate(() => (window as any).grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
