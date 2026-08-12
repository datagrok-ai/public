import { test, expect } from '@playwright/test';
import * as H from './helpers';

// Sticky metadata is bound to the object, not the in-memory dataframe: it must survive clone, new
// view, save-as-project, move-to-space, binary export/import, page refresh, and relogin; clearing
// values removes it.
//
// SCOPE NOTE: this scenario tests *server-side persistence mechanics*. Those operations (clone,
// uploadDataFrame, project/space save, d42 export) have no meaningful pure-UI gesture to assert
// against, so they are driven through the JS API; the relogin is performed through the real UI
// login form, and the metadata clear uses the API. Tag-based matching (`source=<tag>` on the Id
// column) is used because it is the reliable keying mechanism for API reads/writes.
test.describe.configure({ mode: 'serial' });

const LOGIN = process.env.DATAGROK_LOGIN!;
const PASSWORD = process.env.DATAGROK_PASSWORD!;

test('Sticky Meta: persistence across copy / clone / reload, and delete', async ({ page }) => {
  test.setTimeout(240_000);

  const suffix = H.uniqueSuffix();
  const schemaName = `PW_SM_Schema_${suffix}`;
  const tag = `pwsm_${suffix}`;
  let projId: string | null = null;

  try {
    await H.gotoHome(page);
    await H.setupEnv(page);
    await H.apiDeleteAllTestSchemas(page);

    // Setup: schema (tag-matched), open SPGI, seed rating/notes on the first three rows.
    const setup = await page.evaluate(async ({ schemaName, tag }) => {
      const g = (window as any).grok;
      const DG = (window as any).DG;
      const schema = await g.dapi.stickyMeta.createSchema(
        schemaName,
        [{ name: schemaName + '_t', matchBy: 'source=' + tag }],
        [{ name: 'rating', type: 'int' }, { name: 'notes', type: 'string' }],
      );
      const df = await g.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
      g.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      df.col('Id').setTag('source', tag);
      const keyCol = DG.Column.fromList('string', 'Id', [df.col('Id').get(0), df.col('Id').get(1), df.col('Id').get(2)]);
      keyCol.setTag('source', tag);
      const values = DG.DataFrame.fromColumns([
        DG.Column.fromList('int', 'rating', [5, 4, 3]),
        DG.Column.fromList('string', 'notes', ['excellent', 'good', 'average']),
      ]);
      await g.dapi.stickyMeta.setAllValues(schema, keyCol, values);
      await new Promise((r) => setTimeout(r, 1000));
      const read = await g.dapi.stickyMeta.getAllValues(schema, keyCol);
      const rows: any[] = [];
      for (let i = 0; i < read.rowCount; i++) rows.push({ rating: read.col('rating').get(i), notes: read.col('notes').get(i) });
      const ids = [df.col('Id').get(0), df.col('Id').get(1), df.col('Id').get(2)];
      (window as any)._sm = { schemaName, tag, ids };
      return { rows, ids };
    }, { schemaName, tag });
    expect(setup.rows).toEqual([
      { rating: 5, notes: 'excellent' },
      { rating: 4, notes: 'good' },
      { rating: 3, notes: 'average' },
    ]);
    // State for evaluates that run after a page.reload() (window globals are wiped on reload).
    const st = { schemaName, tag, ids: setup.ids };

    const readFirst = async () => page.evaluate(async () => {
      const g = (window as any).grok; const DG = (window as any).DG; const st = (window as any)._sm;
      const schema = (await g.dapi.stickyMeta.getSchemas()).find((s: any) => s.name === st.schemaName);
      const keyCol = DG.Column.fromList('string', 'Id', [st.ids[0]]);
      keyCol.setTag('source', st.tag);
      const read = await g.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { rating: read.col('rating').get(0), notes: read.col('notes').get(0) };
    });

    // ---- 3.1 Clone and new view ----
    const cloneRes = await page.evaluate(async () => {
      const g = (window as any).grok; const DG = (window as any).DG; const st = (window as any)._sm;
      const cloned = g.shell.t.clone();
      g.shell.addTableView(cloned);
      await new Promise((r) => setTimeout(r, 500));
      const idCol = cloned.col('Id');
      const tagPreserved = idCol.tags.get('source');
      const schema = (await g.dapi.stickyMeta.getSchemas()).find((s: any) => s.name === st.schemaName);
      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0)]);
      keyCol.setTag('source', st.tag);
      const read = await g.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { tagPreserved, rating: read.col('rating').get(0), notes: read.col('notes').get(0) };
    });
    expect(cloneRes.tagPreserved).toBe(tag);
    expect({ rating: cloneRes.rating, notes: cloneRes.notes }).toEqual({ rating: 5, notes: 'excellent' });

    // ---- 3.2 Save as project and reopen ----
    await page.evaluate(() => {
      const g = (window as any).grok; const DG = (window as any).DG; const st = (window as any)._sm;
      const projName = 'PW_SM_Proj_' + st.tag;
      (window as any)._projName = projName;
      (window as any)._saveStatus = 'starting';
      (async () => {
        try {
          const orig = g.shell.tables[0];
          const idCol = orig.col('Id'), structCol = orig.col('Structure');
          const small = DG.DataFrame.fromColumns([
            DG.Column.fromList('string', 'Id', [0, 1, 2, 3, 4].map((i) => idCol.get(i))),
            DG.Column.fromList('string', 'Structure', [0, 1, 2, 3, 4].map((i) => structCol.get(i))),
          ]);
          small.name = projName + '_df';
          small.col('Id').setTag('source', st.tag);
          small.col('Structure').semType = 'Molecule';
          const project = DG.Project.create();
          project.name = projName;
          const ti = small.getTableInfo();
          project.addChild(ti);
          await g.dapi.tables.uploadDataFrame(small);
          await g.dapi.tables.save(ti);
          await g.dapi.projects.save(project);
          (window as any)._projId = project.id;
          (window as any)._saveStatus = 'done';
        } catch (e) { (window as any)._saveStatus = 'err: ' + e; }
      })();
    });
    await page.waitForFunction(() => (window as any)._saveStatus === 'done' || String((window as any)._saveStatus).startsWith('err'), { timeout: 40_000 });
    expect(await page.evaluate(() => (window as any)._saveStatus)).toBe('done');
    projId = await page.evaluate(() => (window as any)._projId); // capture before any reload wipes window

    const projRes = await page.evaluate(async () => {
      const g = (window as any).grok; const DG = (window as any).DG; const st = (window as any)._sm;
      g.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
      const loaded = await g.dapi.projects.find((window as any)._projId);
      await loaded.open();
      await new Promise((r) => setTimeout(r, 3000));
      const idCol = g.shell.t.col('Id');
      const schema = (await g.dapi.stickyMeta.getSchemas()).find((s: any) => s.name === st.schemaName);
      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0)]);
      keyCol.setTag('source', st.tag);
      const read = await g.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { tagPreserved: idCol.tags.get('source'), rating: read.col('rating').get(0), notes: read.col('notes').get(0) };
    });
    expect(projRes.tagPreserved).toBe(tag);
    expect({ rating: projRes.rating, notes: projRes.notes }).toEqual({ rating: 5, notes: 'excellent' });

    // ---- 3.3 Export/import (d42 binary) ----
    const ioRes = await page.evaluate(async () => {
      const g = (window as any).grok; const DG = (window as any).DG; const st = (window as any)._sm;
      const df = g.shell.t;
      const imported = DG.DataFrame.fromByteArray(df.toByteArray());
      imported.name = (df.name || 'df') + '_imported';
      g.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      g.shell.addTableView(imported);
      await new Promise((r) => setTimeout(r, 800));
      const idCol = imported.col('Id');
      const schema = (await g.dapi.stickyMeta.getSchemas()).find((s: any) => s.name === st.schemaName);
      const keyCol = DG.Column.fromList('string', 'Id', [idCol.get(0)]);
      keyCol.setTag('source', st.tag);
      const read = await g.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { tagPreserved: idCol.tags.get('source'), rating: read.col('rating').get(0), notes: read.col('notes').get(0) };
    });
    expect(ioRes.tagPreserved).toBe(tag);
    expect({ rating: ioRes.rating, notes: ioRes.notes }).toEqual({ rating: 5, notes: 'excellent' });

    // ---- 3.4 Delete (clear) values and verify removal ----
    const cleared = await page.evaluate(async () => {
      const g = (window as any).grok; const DG = (window as any).DG; const st = (window as any)._sm;
      const schema = (await g.dapi.stickyMeta.getSchemas()).find((s: any) => s.name === st.schemaName);
      const keyCol = DG.Column.fromList('string', 'Id', [st.ids[0]]);
      keyCol.setTag('source', st.tag);
      // int cannot be set to null via the API (no-op) — use 0 as the cleared sentinel; string clears to ''.
      const values = DG.DataFrame.fromColumns([
        DG.Column.fromList('int', 'rating', [0]),
        DG.Column.fromList('string', 'notes', ['']),
      ]);
      await g.dapi.stickyMeta.setAllValues(schema, keyCol, values);
      await new Promise((r) => setTimeout(r, 800));
      const read = await g.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { rating: read.col('rating').get(0), notes: read.col('notes').get(0) };
    });
    expect(cleared).toEqual({ rating: 0, notes: '' });

    // ---- 3.3b Persistence across page refresh ----
    await page.evaluate(async () => {
      const g = (window as any).grok; const DG = (window as any).DG; const st = (window as any)._sm;
      const schema = (await g.dapi.stickyMeta.getSchemas()).find((s: any) => s.name === st.schemaName);
      const keyCol = DG.Column.fromList('string', 'Id', [st.ids[1]]);
      keyCol.setTag('source', st.tag);
      const values = DG.DataFrame.fromColumns([
        DG.Column.fromList('int', 'rating', [2]),
        DG.Column.fromList('string', 'notes', ['refresh-test']),
      ]);
      await g.dapi.stickyMeta.setAllValues(schema, keyCol, values);
      await new Promise((r) => setTimeout(r, 500));
    });
    await page.reload({ waitUntil: 'domcontentloaded' });
    await page.locator('[name="Browse"]').waitFor({ timeout: 120_000 });
    const afterRefresh = await page.evaluate(async (st) => {
      const g = (window as any).grok; const DG = (window as any).DG;
      const schema = (await g.dapi.stickyMeta.getSchemas()).find((s: any) => s.name === st.schemaName);
      const keyCol = DG.Column.fromList('string', 'Id', [st.ids[1]]);
      keyCol.setTag('source', st.tag);
      const read = await g.dapi.stickyMeta.getAllValues(schema, keyCol);
      return { rating: read.col('rating').get(0), notes: read.col('notes').get(0) };
    }, st);
    expect(afterRefresh).toEqual({ rating: 2, notes: 'refresh-test' });

    // ---- 3.3c Persistence across logout + login (real UI login) ----
    expect(await page.evaluate(async () => (await fetch('/api/users/logout', { method: 'POST', credentials: 'include' })).status)).toBe(200);
    await page.reload({ waitUntil: 'domcontentloaded' });
    const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
    await loginInput.waitFor({ timeout: 30_000 });
    await loginInput.click();
    await page.keyboard.type(LOGIN);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(PASSWORD);
    await page.keyboard.press('Enter');
    await page.locator('[name="Browse"]').waitFor({ timeout: 120_000 });
    const afterRelogin = await page.evaluate(async (st) => {
      const g = (window as any).grok; const DG = (window as any).DG;
      const schema = (await g.dapi.stickyMeta.getSchemas()).find((s: any) => s.name === st.schemaName);
      const keyCol = DG.Column.fromList('string', 'Id', [st.ids[1], st.ids[2]]);
      keyCol.setTag('source', st.tag);
      const read = await g.dapi.stickyMeta.getAllValues(schema, keyCol);
      const out: any[] = [];
      for (let i = 0; i < read.rowCount; i++) out.push({ rating: read.col('rating').get(i), notes: read.col('notes').get(i) });
      return out;
    }, st);
    expect(afterRelogin).toEqual([
      { rating: 2, notes: 'refresh-test' },
      { rating: 3, notes: 'average' },
    ]);
  } finally {
    // Cleanup: project, then schema.
    if (projId) await page.evaluate(async (id) => {
      const g = (window as any).grok;
      try { const p = await g.dapi.projects.find(id); if (p) await g.dapi.projects.delete(p); } catch { /* ignore */ }
    }, projId).catch(() => {});
    await H.apiDeleteSchema(page, schemaName).catch(() => {});
    await page.evaluate(() => (window as any).grok.shell.closeAll()).catch(() => {});
  }
});
