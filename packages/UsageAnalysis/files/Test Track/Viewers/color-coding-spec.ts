import {test, expect, chromium} from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Color Coding: linear, categorical, disable, re-enable, copy column', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find(p => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
    await page.waitForFunction(() => {
      try { return typeof grok !== 'undefined' && typeof grok.shell.closeAll === 'function'; }
      catch { return false; }
    }, {timeout: 45000});
  }

  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;

    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  });

  // Steps 2-5: Apply color coding
  await softStep('Steps 2-5: Apply color coding to columns', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      df.col('AGE').meta.colors.setLinear();
      df.col('SEX').meta.colors.setCategorical();
      df.col('CONTROL').meta.colors.setCategorical();
      df.col('STARTED').meta.colors.setLinear();
      await new Promise(r => setTimeout(r, 1000));
      return {
        age: df.col('AGE').meta.colors.getType(),
        sex: df.col('SEX').meta.colors.getType(),
        control: df.col('CONTROL').meta.colors.getType(),
        started: df.col('STARTED').meta.colors.getType(),
      };
    });
    expect(result.age).toBe('Linear');
    expect(result.sex).toBe('Categorical');
    expect(result.control).toBe('Categorical');
    expect(result.started).toBe('Linear');
  });

  // Steps 10-12: Disable and re-enable
  await softStep('Steps 10-12: Disable and re-enable color coding', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      df.col('AGE').meta.colors.setDisabled();
      df.col('SEX').meta.colors.setDisabled();
      await new Promise(r => setTimeout(r, 500));
      const ageOff = df.col('AGE').meta.colors.getType();

      df.col('AGE').meta.colors.setLinear();
      df.col('SEX').meta.colors.setCategorical();
      await new Promise(r => setTimeout(r, 500));
      const ageOn = df.col('AGE').meta.colors.getType();

      return {ageOff, ageOn};
    });
    expect(result.ageOff).toBe('Off');
    expect(result.ageOn).toBe('Linear');
  });

  // Step 15: Copy column with color coding
  await softStep('Step 15: Copy Race column, apply categorical coloring', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      const raceCol = df.col('RACE');
      const raceCopy = df.columns.addNewString('Race_copy');
      for (let i = 0; i < df.rowCount; i++)
        raceCopy.set(i, raceCol.get(i));
      raceCopy.meta.colors.setCategorical();
      await new Promise(r => setTimeout(r, 500));
      return {type: raceCopy.meta.colors.getType(), cols: df.columns.length};
    });
    expect(result.type).toBe('Categorical');
    expect(result.cols).toBe(12);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
