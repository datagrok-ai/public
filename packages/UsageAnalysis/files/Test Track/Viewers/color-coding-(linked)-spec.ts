import {test, expect, chromium} from '@playwright/test';

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

test('Linked Color Coding: link columns, propagate changes, 5-level chain', async () => {
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

  // Steps 1-2: WEIGHT categorical, RACE linked to WEIGHT
  await softStep('Steps 1-2: WEIGHT categorical, RACE linked', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      df.col('WEIGHT').meta.colors.setCategorical();
      df.col('RACE').tags['.color-coding-type'] = 'Linked';
      df.col('RACE').tags['.color-coding-source-column'] = 'WEIGHT';
      await new Promise(r => setTimeout(r, 1000));
      return {
        weightType: df.col('WEIGHT').meta.colors.getType(),
        raceType: df.col('RACE').meta.colors.getType(),
      };
    });
    expect(result.weightType).toBe('Categorical');
    expect(result.raceType).toBe('Linked');
  });

  // Steps 3-4: HEIGHT linked to WEIGHT (text)
  await softStep('Steps 3-4: HEIGHT linked to WEIGHT', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      df.col('HEIGHT').tags['.color-coding-type'] = 'Linked';
      df.col('HEIGHT').tags['.color-coding-source-column'] = 'WEIGHT';
      await new Promise(r => setTimeout(r, 1000));
      return {heightType: df.col('HEIGHT').meta.colors.getType()};
    });
    expect(result.heightType).toBe('Linked');
  });

  // Step 5: Change WEIGHT to Linear, verify links persist
  await softStep('Step 5: Source changes propagate', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      df.col('WEIGHT').meta.colors.setLinear();
      await new Promise(r => setTimeout(r, 1000));
      return {
        weightType: df.col('WEIGHT').meta.colors.getType(),
        raceType: df.col('RACE').meta.colors.getType(),
        heightType: df.col('HEIGHT').meta.colors.getType(),
      };
    });
    expect(result.weightType).toBe('Linear');
    expect(result.raceType).toBe('Linked');
    expect(result.heightType).toBe('Linked');
  });

  // Step 7: 5-level chain
  await softStep('Step 7: 5-level linking chain', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      df.col('AGE').meta.colors.setLinear();
      df.col('SEX').tags['.color-coding-type'] = 'Linked';
      df.col('SEX').tags['.color-coding-source-column'] = 'AGE';
      df.col('DIS_POP').tags['.color-coding-type'] = 'Linked';
      df.col('DIS_POP').tags['.color-coding-source-column'] = 'SEX';
      df.col('CONTROL').tags['.color-coding-type'] = 'Linked';
      df.col('CONTROL').tags['.color-coding-source-column'] = 'DIS_POP';
      await new Promise(r => setTimeout(r, 1000));
      return {
        sexType: df.col('SEX').meta.colors.getType(),
        disPopType: df.col('DIS_POP').meta.colors.getType(),
        controlType: df.col('CONTROL').meta.colors.getType(),
      };
    });
    expect(result.sexType).toBe('Linked');
    expect(result.disPopType).toBe('Linked');
    expect(result.controlType).toBe('Linked');
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
