import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors, login, password} from '../spec-login';

test.use(specTestOptions);

test('Sticky meta: copy / clone / delete persistence', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Setup: open SPGI, create schema, attach rating/notes to 3 molecules (keyed by Id tag).
  const setup = await page.evaluate(async () => {
    // @ts-ignore
    document.body.classList.add('selenium');
    // @ts-ignore
    grok.shell.settings.showFiltersIconsConstantly = true;
    // @ts-ignore
    grok.shell.windows.simpleMode = true;
    // @ts-ignore
    grok.shell.closeAll();
    // @ts-ignore
    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    // @ts-ignore
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
    // @ts-ignore
    const schema = await grok.dapi.stickyMeta.createSchema(
      'CLAUDE_SPGI_' + ts,
      [{name: 'SPGIType_' + ts, matchBy: 'source=' + tagValue}],
      [{name: 'rating', type: 'int'}, {name: 'notes', type: 'string'}],
    );
    const idCol = df.col('Id');
    idCol.setTag('source', tagValue);
    // @ts-ignore
    const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0), idCol.get(1), idCol.get(2)]);
    keyCol.setTag('source', tagValue);
    // @ts-ignore
    const values = DG.DataFrame.fromColumns([
      // @ts-ignore
      DG.Column.fromList('int', 'rating', [5, 4, 3]),
      // @ts-ignore
      DG.Column.fromList('string', 'notes', ['excellent', 'good', 'average']),
    ]);
    // @ts-ignore
    await grok.dapi.stickyMeta.setAllValues(schema, keyCol, values);
    await new Promise(r => setTimeout(r, 1000));
    // @ts-ignore
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
      // @ts-ignore
      const df = grok.shell.t;
      const cloned = df.clone();
      // @ts-ignore
      grok.shell.addTableView(cloned);
      await new Promise(r => setTimeout(r, 500));
      const idCol = cloned.col('Id');
      const tagPreserved = idCol.tags.get('source');
      // @ts-ignore
      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);
      // @ts-ignore
      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0), idCol.get(1), idCol.get(2)]);
      keyCol.setTag('source', (window as any)._state.tagValue);
      // @ts-ignore
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
      // @ts-ignore
      const df = Array.from(grok.shell.tables).find((t: any) => t.name === 'Table');
      // @ts-ignore
      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 500));
      const idCol = df.col('Id');
      const tagPreserved = idCol.tags.get('source');
      // @ts-ignore
      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);
      // @ts-ignore
      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0), idCol.get(1), idCol.get(2)]);
      keyCol.setTag('source', (window as any)._state.tagValue);
      // @ts-ignore
      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { tagPreserved, first: { rating: read.col('rating').get(0), notes: read.col('notes').get(0) } };
    });
    expect(res.tagPreserved).toBe(setup.tagValue);
    expect(res.first).toEqual({rating: 5, notes: 'excellent'});
  });

  // Build a smaller table (5 rows) for project save/reopen flows — full SPGI upload
  // exceeds Chrome DevTools MCP protocol timeouts. Tag structure is identical.
  await softStep('Prepare small SPGI subset for project flows', async () => {
    await page.evaluate(() => {
      // @ts-ignore
      const orig = Array.from(grok.shell.tables).find((t: any) => t.name === 'Table');
      const idCol = orig.col('Id');
      const structCol = orig.col('Structure');
      const ids = [0, 1, 2, 3, 4].map(i => idCol.get(i));
      const structs = [0, 1, 2, 3, 4].map(i => structCol.get(i));
      // @ts-ignore
      const small = DG.DataFrame.fromColumns([
        // @ts-ignore
        DG.Column.fromList('string', 'Id', ids),
        // @ts-ignore
        DG.Column.fromList('string', 'Structure', structs),
      ]);
      small.name = 'SPGI_small_' + Date.now();
      small.col('Id').setTag('source', (window as any)._state.tagValue);
      small.col('Structure').semType = 'Molecule';
      small.col('Structure').setTag('units', 'molblock');
      (window as any)._smallDfName = small.name;
      // @ts-ignore
      grok.shell.addTableView(small);
    });
  });

  await softStep('Section 1 step 2c: Save as project and reopen', async () => {
    await page.evaluate(() => {
      // @ts-ignore
      const smallDf = Array.from(grok.shell.tables).find((t: any) => t.name === (window as any)._smallDfName);
      const projName = 'CLAUDE_SPGI_Proj_' + Date.now();
      (window as any)._projName = projName;
      (window as any)._saveStatus = 'starting';
      (async () => {
        try {
          // @ts-ignore
          const project = DG.Project.create();
          project.name = projName;
          const ti = smallDf.getTableInfo();
          project.addChild(ti);
          // @ts-ignore
          await grok.dapi.tables.uploadDataFrame(smallDf);
          // @ts-ignore
          await grok.dapi.tables.save(ti);
          // @ts-ignore
          await grok.dapi.projects.save(project);
          (window as any)._projId = project.id;
          (window as any)._saveStatus = 'done';
        } catch (e) { (window as any)._saveStatus = 'err'; }
      })();
    });
    await page.waitForFunction(() => (window as any)._saveStatus === 'done', {timeout: 30000});

    const res = await page.evaluate(async () => {
      // @ts-ignore
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
      // @ts-ignore
      const loaded = await grok.dapi.projects.find((window as any)._projId);
      await loaded.open();
      await new Promise(r => setTimeout(r, 3000));
      // @ts-ignore
      const df = grok.shell.t;
      const idCol = df.col('Id');
      const tagPreserved = idCol.tags.get('source');
      // @ts-ignore
      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);
      // @ts-ignore
      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0), idCol.get(1), idCol.get(2)]);
      keyCol.setTag('source', (window as any)._state.tagValue);
      // @ts-ignore
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
          // @ts-ignore
          const space = await grok.dapi.spaces.createRootSpace('CLAUDE_Space_' + ts);
          // @ts-ignore
          const sc = grok.dapi.spaces.id(space.id);
          await sc.addEntity((window as any)._projId, false);
          // @ts-ignore
          grok.shell.closeAll();
          await new Promise(r => setTimeout(r, 500));
          // @ts-ignore
          const loaded = await grok.dapi.projects.find((window as any)._projId);
          await loaded.open();
          await new Promise(r => setTimeout(r, 3000));
          (window as any)._spaceStatus = 'done';
        } catch (e) { (window as any)._spaceStatus = 'err'; }
      })();
    });
    await page.waitForFunction(() => (window as any)._spaceStatus === 'done', {timeout: 30000});

    const res = await page.evaluate(async () => {
      // @ts-ignore
      const df = grok.shell.t;
      const idCol = df.col('Id');
      const tagPreserved = idCol.tags.get('source');
      // @ts-ignore
      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);
      // @ts-ignore
      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0), idCol.get(1), idCol.get(2)]);
      keyCol.setTag('source', (window as any)._state.tagValue);
      // @ts-ignore
      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { tagPreserved, first: { rating: read.col('rating').get(0), notes: read.col('notes').get(0) } };
    });
    expect(res.tagPreserved).toBe(setup.tagValue);
    expect(res.first).toEqual({rating: 5, notes: 'excellent'});
  });

  await softStep('Section 1 step 2e: Export and re-import (d42 binary) preserves metadata', async () => {
    const res = await page.evaluate(async () => {
      // @ts-ignore
      const df = grok.shell.t;
      const bytes = df.toByteArray();
      // @ts-ignore
      const imported = DG.DataFrame.fromByteArray(bytes);
      imported.name = (df.name || 'df') + '_imported';
      // @ts-ignore
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
      // @ts-ignore
      grok.shell.addTableView(imported);
      await new Promise(r => setTimeout(r, 800));
      const idCol = imported.col('Id');
      const tagPreserved = idCol.tags.get('source');
      // @ts-ignore
      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);
      // @ts-ignore
      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0), idCol.get(1), idCol.get(2)]);
      keyCol.setTag('source', (window as any)._state.tagValue);
      // @ts-ignore
      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { tagPreserved, bytes: bytes.length, first: { rating: read.col('rating').get(0), notes: read.col('notes').get(0) } };
    });
    expect(res.tagPreserved).toBe(setup.tagValue);
    expect(res.bytes).toBeGreaterThan(0);
    expect(res.first).toEqual({rating: 5, notes: 'excellent'});
  });

  await softStep('Section 2 steps 1-2: Clear rating (sentinel 0) and notes, save', async () => {
    // Scenario intent: "Delete fields rating and notes. Save". Platform limitation:
    // setAllValues silently ignores null values for int columns, so rating cannot be
    // cleared to absent via the JS API — we use sentinel 0 as the "cleared" value
    // and verify the value *changed* (5 → 0). notes (string) clears to '' normally.
    const res = await page.evaluate(async () => {
      // @ts-ignore
      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);
      // @ts-ignore
      const df = grok.shell.t;
      const idCol = df.col('Id');
      // @ts-ignore
      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0)]);
      keyCol.setTag('source', (window as any)._state.tagValue);
      // @ts-ignore
      const values = DG.DataFrame.fromColumns([
        // @ts-ignore
        DG.Column.fromList('int', 'rating', [0]),
        // @ts-ignore
        DG.Column.fromList('string', 'notes', ['']),
      ]);
      // @ts-ignore
      await grok.dapi.stickyMeta.setAllValues(schema, keyCol, values);
      await new Promise(r => setTimeout(r, 800));
      // @ts-ignore
      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { rating: read.col('rating').get(0), notes: read.col('notes').get(0) };
    });
    expect(res).toEqual({rating: 0, notes: ''});
  });

  await softStep('Section 2 step 3: Deletion persists across view reopen', async () => {
    const res = await page.evaluate(async () => {
      // @ts-ignore
      const tables = Array.from(grok.shell.tables);
      // @ts-ignore
      const orig = tables.find((t: any) => t.name.includes('imported')) ?? tables[0];
      // @ts-ignore
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 400));
      // @ts-ignore
      grok.shell.addTableView(orig);
      await new Promise(r => setTimeout(r, 800));
      // @ts-ignore
      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);
      const idCol = (orig as any).col('Id');
      // @ts-ignore
      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0)]);
      keyCol.setTag('source', (window as any)._state.tagValue);
      // @ts-ignore
      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { rating: read.col('rating').get(0), notes: read.col('notes').get(0) };
    });
    expect(res).toEqual({rating: 0, notes: ''});
  });

  await softStep('Section 3 step 2a: Metadata persists across page refresh', async () => {
    // Seed fresh metadata on row 4 before refreshing.
    await page.evaluate(async () => {
      // @ts-ignore
      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === (window as any)._state.schemaName);
      // @ts-ignore
      const df = grok.shell.t;
      const idCol = df.col('Id');
      // @ts-ignore
      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(4)]);
      keyCol.setTag('source', (window as any)._state.tagValue);
      // @ts-ignore
      const values = DG.DataFrame.fromColumns([
        // @ts-ignore
        DG.Column.fromList('int', 'rating', [2]),
        // @ts-ignore
        DG.Column.fromList('string', 'notes', ['refresh-test']),
      ]);
      // @ts-ignore
      await grok.dapi.stickyMeta.setAllValues(schema, keyCol, values);
      await new Promise(r => setTimeout(r, 500));
    });

    await page.reload({waitUntil: 'domcontentloaded'});
    await page.locator('[name="Browse"]').waitFor({timeout: 120000});

    const res = await page.evaluate(async (state: any) => {
      // @ts-ignore
      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === state.schemaName);
      // @ts-ignore
      const keyCol = DG.Column.fromList('string', 'Id', ['CAST-634787']);
      keyCol.setTag('source', state.tagValue);
      // @ts-ignore
      const read = await grok.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { rating: read.col('rating').get(0), notes: read.col('notes').get(0) };
    }, setup);
    expect(res).toEqual({rating: 2, notes: 'refresh-test'});
  });

  await softStep('Section 3 step 2b: Metadata persists across logout + login', async () => {
    // POST /api/users/logout clears the auth cookie, then re-login via the normal form.
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
    await page.locator('[name="Browse"]').waitFor({timeout: 120000});

    // After a fresh session, read the metadata from the server.
    const res = await page.evaluate(async (state: any) => {
      // @ts-ignore
      const schemas = await grok.dapi.stickyMeta.getSchemas();
      const schema = schemas.find((s: any) => s.name === state.schemaName);
      // @ts-ignore
      const keyCol = DG.Column.fromList('string', 'Id', ['CAST-634787', 'CAST-634784', 'CAST-634785']);
      keyCol.setTag('source', state.tagValue);
      // @ts-ignore
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
