/* ---
realizes: [viewers.grid]
--- */
import {expect} from '@playwright/test';
import {localTest as test} from '../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../spec-login';
import * as v from '../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';

test('Color Coding: types, disable/re-enable, pick-up/apply, linked, scheme invert', async ({page}) => {
  test.setTimeout(120_000);

  await openDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
  // every write below is read back synchronously through getType(); the grid repaint each one
  // triggers is what the old fixed settles covered, so it is awaited with the same cap instead
  await page.evaluate(() => {
    const w = window as any;
    w.__painted = (act: () => void, capMs: number) => w.__settled('viewer:Grid.onAfterDrawContent', act, capMs);
  });

  await softStep('2.1 AGE: linear then conditional', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.t;
      await w.__painted(() => df.col('AGE').meta.colors.setLinear(), 300);
      const linearType = df.col('AGE').meta.colors.getType();

      await w.__painted(() => df.col('AGE').meta.colors.setConditional({
        '< 30': DG.Color.fromHtml('#00CC44'),
        '30-60': DG.Color.fromHtml('#FFCC00'),
        '> 60': DG.Color.fromHtml('#FF4444'),
      }), 300);
      const condType = df.col('AGE').meta.colors.getType();

      return {linearType, condType};
    });
    expect(result.linearType).toBe('Linear');
    expect(result.condType).toBe('Conditional');
  });

  await softStep('2.2 SEX: categorical with custom M/F colors', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.t;
      await w.__painted(() => df.col('SEX').meta.colors.setCategorical({
        'M': DG.Color.fromHtml('#3366CC'),
        'F': DG.Color.fromHtml('#CC6699'),
      }), 300);
      return {type: df.col('SEX').meta.colors.getType()};
    });
    expect(result.type).toBe('Categorical');
  });

  await softStep('2.3 CONTROL: categorical (default colors)', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.t;
      await w.__painted(() => df.col('CONTROL').meta.colors.setCategorical(), 300);
      return {type: df.col('CONTROL').meta.colors.getType()};
    });
    expect(result.type).toBe('Categorical');
  });

  await softStep('2.4 STARTED: linear with custom 3-stop scheme', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.t;
      await w.__painted(() => df.col('STARTED').meta.colors.setLinear([
        DG.Color.fromHtml('#0000FF'),
        DG.Color.fromHtml('#FFFFFF'),
        DG.Color.fromHtml('#FF0000'),
      ]), 300);
      return {type: df.col('STARTED').meta.colors.getType()};
    });
    expect(result.type).toBe('Linear');
  });

  await softStep('3.1 Disable AGE, SEX, STARTED', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.t;
      await w.__painted(() => {
        df.col('AGE').meta.colors.setDisabled();
        df.col('SEX').meta.colors.setDisabled();
        df.col('STARTED').meta.colors.setDisabled();
      }, 300);
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
    const result = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.t;

      const ageBefore: Record<string, string> = {};
      for (const [k, val] of Object.entries(df.col('AGE').tags))
        if (k.startsWith('.color-coding')) ageBefore[k] = val as string;

      await w.__painted(() => {
        df.col('AGE').meta.colors.setConditional();
        df.col('SEX').meta.colors.setCategorical();
        df.col('STARTED').meta.colors.setLinear();
      }, 300);

      const ageAfter: Record<string, string> = {};
      for (const [k, val] of Object.entries(df.col('AGE').tags))
        if (k.startsWith('.color-coding')) ageAfter[k] = val as string;

      return {
        ageType: df.col('AGE').meta.colors.getType(),
        sexType: df.col('SEX').meta.colors.getType(),
        startedType: df.col('STARTED').meta.colors.getType(),
        tagsPreserved: Object.keys(ageBefore)
          .filter((k) => k !== '.color-coding-type')
          .every((k) => ageAfter[k] === ageBefore[k]),
      };
    });
    expect(result.ageType).toBe('Conditional');
    expect(result.sexType).toBe('Categorical');
    expect(result.startedType).toBe('Linear');
    expect(result.tagsPreserved).toBe(true);
  });

  await softStep('4.1 Create Race_copy, apply categorical', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.t;
      df.col('RACE').meta.colors.setCategorical();
      const raceCol = df.col('RACE');
      const raceCopy = df.columns.addNewString('Race_copy');
      for (let i = 0; i < df.rowCount; i++) raceCopy.set(i, raceCol.get(i));
      await w.__painted(() => raceCopy.meta.colors.setCategorical(), 300);
      return {type: raceCopy.meta.colors.getType(), colCount: df.columns.length};
    });
    expect(result.type).toBe('Categorical');
    expect(result.colCount).toBe(12);
  });

  await softStep('4.2 Apply RACE coloring to Race_copy', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.t;
      const src = df.col('RACE');
      const dst = df.col('Race_copy');
      await w.__painted(() => {
        for (const [key, val] of Object.entries(src.tags))
          if (key.startsWith('.color-coding')) dst.tags[key] = val as string;
      }, 300);
      return {srcType: src.meta.colors.getType(), dstType: dst.meta.colors.getType()};
    });
    expect(result.dstType).toBe(result.srcType);
  });

  await softStep('4.3 Apply STARTED coloring to HEIGHT', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.t;
      const src = df.col('STARTED');
      const dst = df.col('HEIGHT');
      await w.__painted(() => {
        for (const [key, val] of Object.entries(src.tags))
          if (key.startsWith('.color-coding')) dst.tags[key] = val as string;
      }, 300);
      return {srcType: src.meta.colors.getType(), dstType: dst.meta.colors.getType()};
    });
    expect(result.dstType).toBe(result.srcType);
  });

  await softStep('5.1 RACE linked to WEIGHT (background)', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.t;
      await w.__painted(() => {
        df.col('WEIGHT').meta.colors.setCategorical();
        df.col('RACE').tags['.color-coding-type'] = 'Linked';
        df.col('RACE').tags['.color-coding-source-column'] = 'WEIGHT';
      }, 500);
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
      const w = window as any;
      const df = grok.shell.t;
      await w.__painted(() => {
        df.col('HEIGHT').tags['.color-coding-type'] = 'Linked';
        df.col('HEIGHT').tags['.color-coding-source-column'] = 'WEIGHT';
      }, 500);
      return {heightType: df.col('HEIGHT').meta.colors.getType()};
    });
    expect(result.heightType).toBe('Linked');
  });

  await softStep('5.3 Source changes propagate (Linear and Conditional)', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.t;

      await w.__painted(() => df.col('WEIGHT').meta.colors.setLinear(), 500);
      const afterLinear = {
        race: df.col('RACE').meta.colors.getType(),
        height: df.col('HEIGHT').meta.colors.getType(),
      };

      await w.__painted(() => df.col('WEIGHT').meta.colors.setConditional({
        '< 60': DG.Color.fromHtml('#00CC44'),
        '> 60': DG.Color.fromHtml('#FF4444'),
      }), 500);
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
      const w = window as any;
      const df = grok.shell.t;
      await w.__painted(() => {
        df.col('AGE').meta.colors.setLinear();
        df.col('SEX').tags['.color-coding-type'] = 'Linked';
        df.col('SEX').tags['.color-coding-source-column'] = 'AGE';
        df.col('DIS_POP').tags['.color-coding-type'] = 'Linked';
        df.col('DIS_POP').tags['.color-coding-source-column'] = 'SEX';
        df.col('CONTROL').tags['.color-coding-type'] = 'Linked';
        df.col('CONTROL').tags['.color-coding-source-column'] = 'DIS_POP';
        df.col('STARTED').tags['.color-coding-type'] = 'Linked';
        df.col('STARTED').tags['.color-coding-source-column'] = 'CONTROL';
      }, 500);
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

  await softStep('6.1 AGE: custom 3-stop linear scheme', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.t;
      await w.__painted(() => df.col('AGE').meta.colors.setLinear([
        DG.Color.fromHtml('#1A237E'),
        DG.Color.fromHtml('#F5F5F5'),
        DG.Color.fromHtml('#B71C1C'),
      ]), 300);
      return {type: df.col('AGE').meta.colors.getType()};
    });
    expect(result.type).toBe('Linear');
  });

  await softStep('6.2 AGE: invert the scheme', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const col = grok.shell.t.col('AGE');
      const tagKey = '.color-coding-linear';
      const raw = col.tags[tagKey];
      if (!raw) throw new Error(`Tag ${tagKey} not found`);
      const colors: number[] = JSON.parse(raw);
      const original = [...colors];
      colors.reverse();
      await w.__painted(() => { col.tags[tagKey] = JSON.stringify(colors); }, 300);
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

  await page.evaluate(() => { delete (window as any).__painted; });
  await v.cleanupShell(page);

  v.finishSpec();
});
