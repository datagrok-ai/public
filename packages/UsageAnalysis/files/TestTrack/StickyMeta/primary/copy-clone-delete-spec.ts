import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors, login, password} from '../../spec-login';

test.use(specTestOptions);

test('Sticky meta: copy / clone / delete persistence', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  const setup = await page.evaluate(async () => {

    document.body.classList.add('selenium');

    grok.shell.settings.showFiltersIconsConstantly = true;

    grok.shell.windows.simpleMode = true;

    grok.shell.closeAll();

    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');

    grok.shell.addTableView(df);
    await new Promise<void>(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const hasChem = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Molecule');
    if (hasChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }

    const ts = Date.now();
    const tagValue = 'SPGI_' + ts;

    const schema = await grok.dapi.stickyMeta.createSchema(
      'CLAUDE_SPGI_' + ts,
      [{name: 'SPGIType_' + ts, matchBy: 'source=' + tagValue}],
      [{name: 'rating', type: 'int'}, {name: 'notes', type: 'string'}],
    );
    const idCol = df.col('Id');
    idCol.setTag('source', tagValue);

    const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0), idCol.get(1), idCol.get(2)]);
    keyCol.setTag('source', tagValue);

    const values = DG.DataFrame.fromColumns([

      DG.Column.fromList('int', 'rating', [5, 4, 3]),

      DG.Column.fromList('string', 'notes', ['excellent', 'good', 'average']),
    ]);

    await grok.dapi.stickyMeta.setAllValues(schema, keyCol, values);
    await new Promise(r => setTimeout(r, 1000));

    const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
    const rows: any[] = [];
    for (let i = 0; i < read.rowCount; i++)
      rows.push({rating: read.col('rating').get(i), notes: read.col('notes').get(i)});
    (window as any)._state = { schemaName: schema.name, tagValue };
    return { rows, tagValue, schemaName: schema.name };
  });
  expect(setup.rows).toEqual([
    {rating: 5, notes: 'excellent'},
    {rating: 4, notes: 'good'},
    {rating: 3, notes: 'average'},
  ]);

  await softStep('Section 1 step 1: Open SPGI with sticky meta (setup)', async () => {
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 30000});
    expect(setup.tagValue).toBeTruthy();
  });

  await softStep('Section 1 step 2a: Cloned table preserves metadata', async () => {
    const res = await page.evaluate(async () => {

      const df = grok.shell.t;
      const cloned = df.clone();

      grok.shell.addTableView(cloned);
      await new Promise(r => setTimeout(r, 500));
      const idCol = cloned.col('Id');
      const tagPreserved = idCol.tags.get('source');

      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);

      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0), idCol.get(1), idCol.get(2)]);
      keyCol.setTag('source', (window as any)._state.tagValue);

      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      const rows: any[] = [];
      for (let i = 0; i < read.rowCount; i++)
        rows.push({rating: read.col('rating').get(i), notes: read.col('notes').get(i)});
      return { tagPreserved, rows };
    });
    expect(res.tagPreserved).toBe(setup.tagValue);
    expect(res.rows[0]).toEqual({rating: 5, notes: 'excellent'});
  });

  await softStep('Section 1 step 2b: New view on same DataFrame shows metadata', async () => {
    const res = await page.evaluate(async () => {

      const df = Array.from(grok.shell.tables).find((t: any) => t.name === 'Table');

      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 500));
      const idCol = df.col('Id');
      const tagPreserved = idCol.tags.get('source');

      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);

      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0), idCol.get(1), idCol.get(2)]);
      keyCol.setTag('source', (window as any)._state.tagValue);

      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { tagPreserved, first: { rating: read.col('rating').get(0), notes: read.col('notes').get(0) } };
    });
    expect(res.tagPreserved).toBe(setup.tagValue);
    expect(res.first).toEqual({rating: 5, notes: 'excellent'});
  });

  await softStep('Prepare small SPGI subset for project flows', async () => {
    await page.evaluate(() => {

      const orig = Array.from(grok.shell.tables).find((t: any) => t.name === 'Table');
      const idCol = orig.col('Id');
      const structCol = orig.col('Structure');
      const ids = [0, 1, 2, 3, 4].map(i => idCol.get(i));
      const structs = [0, 1, 2, 3, 4].map(i => structCol.get(i));

      const small = DG.DataFrame.fromColumns([

        DG.Column.fromList('string', 'Id', ids),

        DG.Column.fromList('string', 'Structure', structs),
      ]);
      small.name = 'SPGI_small_' + Date.now();
      small.col('Id').setTag('source', (window as any)._state.tagValue);
      small.col('Structure').semType = 'Molecule';
      small.col('Structure').setTag('units', 'molblock');
      (window as any)._smallDfName = small.name;

      grok.shell.addTableView(small);
    });
  });

  await softStep('Section 1 step 2c: Save as project and reopen', async () => {
    await page.evaluate(() => {

      const smallDf = Array.from(grok.shell.tables).find((t: any) => t.name === (window as any)._smallDfName);
      const projName = 'CLAUDE_SPGI_Proj_' + Date.now();
      (window as any)._projName = projName;
      (window as any)._saveStatus = 'starting';
      (async () => {
        try {

          const project = DG.Project.create();
          project.name = projName;
          const ti = smallDf.getTableInfo();
          project.addChild(ti);

          await grok.dapi.tables.uploadDataFrame(smallDf);

          await grok.dapi.tables.save(ti);

          await grok.dapi.projects.save(project);
          (window as any)._projId = project.id;
          (window as any)._saveStatus = 'done';
        } catch (e) { (window as any)._saveStatus = 'err'; }
      })();
    });
    await page.waitForFunction(() => (window as any)._saveStatus === 'done', {timeout: 30000});

    const res = await page.evaluate(async () => {

      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));

      const loaded = await grok.dapi.projects.find((window as any)._projId);
      await loaded.open();
      await new Promise(r => setTimeout(r, 3000));

      const df = grok.shell.t;
      const idCol = df.col('Id');
      const tagPreserved = idCol.tags.get('source');

      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);

      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0), idCol.get(1), idCol.get(2)]);
      keyCol.setTag('source', (window as any)._state.tagValue);

      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { tagPreserved, first: { rating: read.col('rating').get(0), notes: read.col('notes').get(0) } };
    });
    expect(res.tagPreserved).toBe(setup.tagValue);
    expect(res.first).toEqual({rating: 5, notes: 'excellent'});
  });

  await softStep('Section 1 step 2d: Move project to Space and reopen', async () => {
    await page.evaluate(() => {
      (window as any)._spaceStatus = 'starting';
      (async () => {
        try {
          const ts = Date.now();

          const space = await grok.dapi.spaces.createRootSpace('CLAUDE_Space_' + ts);

          const sc = grok.dapi.spaces.id(space.id);
          await sc.addEntity((window as any)._projId, false);

          grok.shell.closeAll();
          await new Promise(r => setTimeout(r, 500));

          const loaded = await grok.dapi.projects.find((window as any)._projId);
          await loaded.open();
          await new Promise(r => setTimeout(r, 3000));
          (window as any)._spaceStatus = 'done';
        } catch (e) { (window as any)._spaceStatus = 'err'; }
      })();
    });
    await page.waitForFunction(() => (window as any)._spaceStatus === 'done', {timeout: 30000});

    const res = await page.evaluate(async () => {

      const df = grok.shell.t;
      const idCol = df.col('Id');
      const tagPreserved = idCol.tags.get('source');

      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);

      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0), idCol.get(1), idCol.get(2)]);
      keyCol.setTag('source', (window as any)._state.tagValue);

      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { tagPreserved, first: { rating: read.col('rating').get(0), notes: read.col('notes').get(0) } };
    });
    expect(res.tagPreserved).toBe(setup.tagValue);
    expect(res.first).toEqual({rating: 5, notes: 'excellent'});
  });

  await softStep('Section 1 step 2e: Export and re-import (d42 binary) preserves metadata', async () => {
    const res = await page.evaluate(async () => {

      const df = grok.shell.t;
      const bytes = df.toByteArray();

      const imported = DG.DataFrame.fromByteArray(bytes);
      imported.name = (df.name || 'df') + '_imported';

      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));

      grok.shell.addTableView(imported);
      await new Promise(r => setTimeout(r, 800));
      const idCol = imported.col('Id');
      const tagPreserved = idCol.tags.get('source');

      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);

      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0), idCol.get(1), idCol.get(2)]);
      keyCol.setTag('source', (window as any)._state.tagValue);

      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { tagPreserved, bytes: bytes.length, first: { rating: read.col('rating').get(0), notes: read.col('notes').get(0) } };
    });
    expect(res.tagPreserved).toBe(setup.tagValue);
    expect(res.bytes).toBeGreaterThan(0);
    expect(res.first).toEqual({rating: 5, notes: 'excellent'});
  });

  await softStep('Section 2 steps 1-2: Clear rating (sentinel 0) and notes, save', async () => {

    const res = await page.evaluate(async () => {

      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);

      const df = grok.shell.t;
      const idCol = df.col('Id');

      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0)]);
      keyCol.setTag('source', (window as any)._state.tagValue);

      const values = DG.DataFrame.fromColumns([

        DG.Column.fromList('int', 'rating', [0]),

        DG.Column.fromList('string', 'notes', ['']),
      ]);

      await grok.dapi.stickyMeta.setAllValues(schema, keyCol, values);
      await new Promise(r => setTimeout(r, 800));

      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { rating: read.col('rating').get(0), notes: read.col('notes').get(0) };
    });
    expect(res).toEqual({rating: 0, notes: ''});
  });

  await softStep('Section 2 step 3: Deletion persists across view reopen', async () => {
    const res = await page.evaluate(async () => {

      const tables = Array.from(grok.shell.tables);

      const orig = tables.find((t: any) => t.name.includes('imported')) ?? tables[0];

      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 400));

      grok.shell.addTableView(orig);
      await new Promise(r => setTimeout(r, 800));

      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);
      const idCol = (orig as any).col('Id');

      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0)]);
      keyCol.setTag('source', (window as any)._state.tagValue);

      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { rating: read.col('rating').get(0), notes: read.col('notes').get(0) };
    });
    expect(res).toEqual({rating: 0, notes: ''});
  });

  await softStep('Section 3 step 2a: Metadata persists across page refresh', async () => {

    await page.evaluate(async () => {

      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);

      const df = grok.shell.t;
      const idCol = df.col('Id');

      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(4)]);
      keyCol.setTag('source', (window as any)._state.tagValue);

      const values = DG.DataFrame.fromColumns([

        DG.Column.fromList('int', 'rating', [2]),

        DG.Column.fromList('string', 'notes', ['refresh-test']),
      ]);

      await grok.dapi.stickyMeta.setAllValues(schema, keyCol, values);
      await new Promise(r => setTimeout(r, 500));
    });

    await page.reload({waitUntil: 'domcontentloaded'});
    await page.locator('[name="Browse"]').waitFor({timeout: 60_000});

    const res = await page.evaluate(async (state: any) => {

      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === state.schemaName);

      const keyCol = DG.Column.fromList('string', 'Id', ['CAST-634787']);
      keyCol.setTag('source', state.tagValue);

      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { rating: read.col('rating').get(0), notes: read.col('notes').get(0) };
    }, setup);
    expect(res).toEqual({rating: 2, notes: 'refresh-test'});
  });

  await softStep('Section 3 step 2b: Metadata persists across logout + login', async () => {

    const logoutResp = await page.evaluate(async () => {
      const r = await fetch('/api/users/logout', {method: 'POST', credentials: 'include'});
      return r.status;
    });
    expect(logoutResp).toBe(200);

    await page.reload({waitUntil: 'domcontentloaded'});
    const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
    await loginInput.waitFor({timeout: 30000});
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
    await page.locator('[name="Browse"]').waitFor({timeout: 60_000});

    const res = await page.evaluate(async (state: any) => {

      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === state.schemaName);

      const keyCol = DG.Column.fromList('string', 'Id', ['CAST-634787', 'CAST-634784', 'CAST-634785']);
      keyCol.setTag('source', state.tagValue);

      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      const rows: any[] = [];
      for (let i = 0; i < read.rowCount; i++)
        rows.push({rating: read.col('rating').get(i), notes: read.col('notes').get(i)});
      return rows;
    }, setup);
    expect(res).toEqual([
      {rating: 2, notes: 'refresh-test'},
      {rating: 4, notes: 'good'},
      {rating: 3, notes: 'average'},
    ]);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  • ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
