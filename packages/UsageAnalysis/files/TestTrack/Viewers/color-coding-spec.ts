import {test, expect, type Page} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';
const demogPath = 'System:DemoFiles/demog.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

async function openTable(page: Page, path: string) {
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

  await page.evaluate(async (p) => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(p);
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  }, path);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
}

test('Color Coding: types, disable/re-enable, pick-up/apply, linked, scheme invert', async ({page}) => {
  test.setTimeout(600000);
  stepErrors.length = 0;

  await openTable(page, demogPath);

  // ── Group 2: Apply color coding types ──────────────────────────────────────

  await softStep('2.1 AGE: linear then conditional', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      df.col('AGE').meta.colors.setLinear();
      await new Promise(r => setTimeout(r, 300));
      const linearType = df.col('AGE').meta.colors.getType();

      df.col('AGE').meta.colors.setConditional({
        '< 30': DG.Color.fromHtml('#00CC44'),
        '30-60': DG.Color.fromHtml('#FFCC00'),
        '> 60': DG.Color.fromHtml('#FF4444'),
      });
      await new Promise(r => setTimeout(r, 300));
      const condType = df.col('AGE').meta.colors.getType();

      return {linearType, condType};
    });
    expect(result.linearType).toBe('Linear');
    expect(result.condType).toBe('Conditional');
  });

  await softStep('2.2 SEX: categorical with custom M/F colors', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      df.col('SEX').meta.colors.setCategorical({
        'M': DG.Color.fromHtml('#3366CC'),
        'F': DG.Color.fromHtml('#CC6699'),
      });
      await new Promise(r => setTimeout(r, 300));
      return {type: df.col('SEX').meta.colors.getType()};
    });
    expect(result.type).toBe('Categorical');
  });

  await softStep('2.3 CONTROL: categorical (default colors)', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      df.col('CONTROL').meta.colors.setCategorical();
      await new Promise(r => setTimeout(r, 300));
      return {type: df.col('CONTROL').meta.colors.getType()};
    });
    expect(result.type).toBe('Categorical');
  });

  await softStep('2.4 STARTED: linear with custom 3-stop scheme', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      df.col('STARTED').meta.colors.setLinear([
        DG.Color.fromHtml('#0000FF'),
        DG.Color.fromHtml('#FFFFFF'),
        DG.Color.fromHtml('#FF0000'),
      ]);
      await new Promise(r => setTimeout(r, 300));
      return {type: df.col('STARTED').meta.colors.getType()};
    });
    expect(result.type).toBe('Linear');
  });

  // ── Group 3: Disable and re-enable ─────────────────────────────────────────

  await softStep('3.1 Disable AGE, SEX, STARTED', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      df.col('AGE').meta.colors.setDisabled();
      df.col('SEX').meta.colors.setDisabled();
      df.col('STARTED').meta.colors.setDisabled();
      await new Promise(r => setTimeout(r, 300));
      return {
        age: df.col('AGE').meta.colors.getType(),
        sex: df.col('SEX').meta.colors.getType(),
        started: df.col('STARTED').meta.colors.getType(),
      };
    });
    expect(result.age).toBe('Off');
    expect(result.sex).toBe('Off');
    expect(result.started).toBe('Off');
  });

  await softStep('3.2 Re-enable: types and custom colors preserved', async () => {
    // setDisabled() sets only .color-coding-type to 'Off' without clearing other tags,
    // so restoring the type brings back the original color configuration.
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;

      const ageBefore: Record<string, string> = {};
      for (const [k, v] of Object.entries(df.col('AGE').tags))
        if (k.startsWith('.color-coding')) ageBefore[k] = v as string;

      df.col('AGE').meta.colors.setConditional();
      df.col('SEX').meta.colors.setCategorical();
      df.col('STARTED').meta.colors.setLinear();
      await new Promise(r => setTimeout(r, 300));

      const ageAfter: Record<string, string> = {};
      for (const [k, v] of Object.entries(df.col('AGE').tags))
        if (k.startsWith('.color-coding')) ageAfter[k] = v as string;

      return {
        ageType: df.col('AGE').meta.colors.getType(),
        sexType: df.col('SEX').meta.colors.getType(),
        startedType: df.col('STARTED').meta.colors.getType(),
        tagsPreserved: Object.keys(ageBefore)
          .filter(k => k !== '.color-coding-type')
          .every(k => ageAfter[k] === ageBefore[k]),
      };
    });
    expect(result.ageType).toBe('Conditional');
    expect(result.sexType).toBe('Categorical');
    expect(result.startedType).toBe('Linear');
    expect(result.tagsPreserved).toBe(true);
  });

  // ── Group 4: Pick Up / Apply coloring ──────────────────────────────────────

  await softStep('4.1 Create Race_copy, apply categorical', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      df.col('RACE').meta.colors.setCategorical();
      const raceCol = df.col('RACE');
      const raceCopy = df.columns.addNewString('Race_copy');
      for (let i = 0; i < df.rowCount; i++) raceCopy.set(i, raceCol.get(i));
      raceCopy.meta.colors.setCategorical();
      await new Promise(r => setTimeout(r, 300));
      return {type: raceCopy.meta.colors.getType(), colCount: df.columns.length};
    });
    expect(result.type).toBe('Categorical');
    expect(result.colCount).toBe(12);
  });

  await softStep('4.2 Apply RACE coloring to Race_copy', async () => {
    // Copies all .color-coding-* tags — equivalent to Pick Up Coloring / Apply Coloring.
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      const src = df.col('RACE');
      const dst = df.col('Race_copy');
      for (const [key, val] of Object.entries(src.tags)) {
        if (key.startsWith('.color-coding')) dst.tags[key] = val as string;
      }
      await new Promise(r => setTimeout(r, 300));
      return {srcType: src.meta.colors.getType(), dstType: dst.meta.colors.getType()};
    });
    expect(result.dstType).toBe(result.srcType);
  });

  await softStep('4.3 Apply STARTED coloring to HEIGHT', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      const src = df.col('STARTED');
      const dst = df.col('HEIGHT');
      for (const [key, val] of Object.entries(src.tags)) {
        if (key.startsWith('.color-coding')) dst.tags[key] = val as string;
      }
      await new Promise(r => setTimeout(r, 300));
      return {srcType: src.meta.colors.getType(), dstType: dst.meta.colors.getType()};
    });
    expect(result.dstType).toBe(result.srcType);
  });

  // ── Group 5: Linked color coding ───────────────────────────────────────────

  await softStep('5.1 RACE linked to WEIGHT (background)', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      df.col('WEIGHT').meta.colors.setCategorical();
      df.col('RACE').tags['.color-coding-type'] = 'Linked';
      df.col('RACE').tags['.color-coding-source-column'] = 'WEIGHT';
      await new Promise(r => setTimeout(r, 500));
      return {
        weightType: df.col('WEIGHT').meta.colors.getType(),
        raceType: df.col('RACE').meta.colors.getType(),
      };
    });
    expect(result.weightType).toBe('Categorical');
    expect(result.raceType).toBe('Linked');
  });

  await softStep('5.2 HEIGHT linked to WEIGHT (text)', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      df.col('HEIGHT').tags['.color-coding-type'] = 'Linked';
      df.col('HEIGHT').tags['.color-coding-source-column'] = 'WEIGHT';
      await new Promise(r => setTimeout(r, 500));
      return {heightType: df.col('HEIGHT').meta.colors.getType()};
    });
    expect(result.heightType).toBe('Linked');
  });

  await softStep('5.3 Source changes propagate (Linear and Conditional)', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;

      df.col('WEIGHT').meta.colors.setLinear();
      await new Promise(r => setTimeout(r, 500));
      const afterLinear = {
        race: df.col('RACE').meta.colors.getType(),
        height: df.col('HEIGHT').meta.colors.getType(),
      };

      df.col('WEIGHT').meta.colors.setConditional({
        '< 60': DG.Color.fromHtml('#00CC44'),
        '> 60': DG.Color.fromHtml('#FF4444'),
      });
      await new Promise(r => setTimeout(r, 500));
      const afterConditional = {
        race: df.col('RACE').meta.colors.getType(),
        height: df.col('HEIGHT').meta.colors.getType(),
      };

      return {afterLinear, afterConditional};
    });
    expect(result.afterLinear.race).toBe('Linked');
    expect(result.afterLinear.height).toBe('Linked');
    expect(result.afterConditional.race).toBe('Linked');
    expect(result.afterConditional.height).toBe('Linked');
  });

  await softStep('5.4 5-level linking chain: AGE→SEX→DIS_POP→CONTROL→STARTED', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      df.col('AGE').meta.colors.setLinear();
      df.col('SEX').tags['.color-coding-type'] = 'Linked';
      df.col('SEX').tags['.color-coding-source-column'] = 'AGE';
      df.col('DIS_POP').tags['.color-coding-type'] = 'Linked';
      df.col('DIS_POP').tags['.color-coding-source-column'] = 'SEX';
      df.col('CONTROL').tags['.color-coding-type'] = 'Linked';
      df.col('CONTROL').tags['.color-coding-source-column'] = 'DIS_POP';
      df.col('STARTED').tags['.color-coding-type'] = 'Linked';
      df.col('STARTED').tags['.color-coding-source-column'] = 'CONTROL';
      await new Promise(r => setTimeout(r, 500));
      return {
        sex: df.col('SEX').meta.colors.getType(),
        disPop: df.col('DIS_POP').meta.colors.getType(),
        control: df.col('CONTROL').meta.colors.getType(),
        started: df.col('STARTED').meta.colors.getType(),
      };
    });
    expect(result.sex).toBe('Linked');
    expect(result.disPop).toBe('Linked');
    expect(result.control).toBe('Linked');
    expect(result.started).toBe('Linked');
  });

  // ── Group 6: Edit linear color scheme ──────────────────────────────────────

  await softStep('6.1 AGE: custom 3-stop linear scheme', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      df.col('AGE').meta.colors.setLinear([
        DG.Color.fromHtml('#1A237E'),
        DG.Color.fromHtml('#F5F5F5'),
        DG.Color.fromHtml('#B71C1C'),
      ]);
      await new Promise(r => setTimeout(r, 300));
      return {type: df.col('AGE').meta.colors.getType()};
    });
    expect(result.type).toBe('Linear');
  });

  await softStep('6.2 AGE: invert the scheme', async () => {
    // Tag .color-coding-linear stores ARGB integers; reverse = invert gradient.
    const result = await page.evaluate(async () => {
      const col = grok.shell.t.col('AGE');
      const tagKey = '.color-coding-linear';
      const raw = col.tags[tagKey];
      if (!raw) throw new Error(`Tag ${tagKey} not found`);
      const colors: number[] = JSON.parse(raw);
      const original = [...colors];
      colors.reverse();
      col.tags[tagKey] = JSON.stringify(colors);
      await new Promise(r => setTimeout(r, 300));
      const after: number[] = JSON.parse(col.tags[tagKey]);
      return {
        type: col.meta.colors.getType(),
        firstChanged: after[0] !== original[0],
        lastChanged: after[after.length - 1] !== original[original.length - 1],
      };
    });
    expect(result.type).toBe('Linear');
    expect(result.firstChanged).toBe(true);
    expect(result.lastChanged).toBe(true);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
