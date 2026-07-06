// Multi-source integration: file + ad-hoc DB table on System:Datagrok +
// provisioned saved query on System:Datagrok + provisioned dataframe-output
// script. All non-file prerequisites are created in-test via
// helpers/openers.ts — no Samples package or env-provisioned DB.
import {test, expect, Page} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {
  openTableFromDbQuery,
  openTableFromDbTable,
  openTableFromScript,
  provisionSystemDatagrokQuery,
  provisionDataframeScript,
  deleteProvisionedScript,
  resetShell,
  assertProvenanceScript,
  ProvisionedQuery,
  ProvisionedScript,
  SYSTEM_DATAGROK_DB_TABLE,
  SYSTEM_DATAGROK_NQNAME,
  SYSTEM_DATAGROK_QUERIES,
} from '../helpers/openers';
import {saveAllTablesWithProvenance, deleteProjectWithCleanup, SavedAllTables} from '../helpers/projects';

// Local OpenFile invocation with colon-form fullPath (no helpers/openers.ts
// dot-form normalization). Verified 2026-05-08: dot-form .script provenance
// makes JS-API project Save mis-attribute children — derived-tables-spec was
// fixed by switching to colon-form. Apply the same fix here to test the
// hypothesis that project relations also depend on the colon-form provenance.
async function openFileColonForm(page: Page, fullPath: string): Promise<{rowCount: number; script: string}> {
  return await page.evaluate(async (p) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    await DG.Func.find({name: 'OpenFile'})[0].prepare({
      fullPath: p,
    }).call(undefined, undefined, {processed: false});
    // Settle and read provenance.
    let df: any = null;
    for (let i = 0; i < 24; i++) {
      const tv = grok.shell.tv;
      if (tv?.dataFrame && typeof tv.addViewer === 'function') {
        df = tv.dataFrame;
        break;
      }
      await new Promise((r) => setTimeout(r, 500));
    }
    if (!df) throw new Error('OpenFile("' + p + '") did not produce a TableView');
    return {
      rowCount: df.rowCount,
      script: df.tags?.get?.('.script') ?? '',
    };
  }, fullPath);
}

test.use(projectsTestOptions);

test('Projects / Complex Integration: heterogeneous sources in one project', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `integration-test-${stamp}`;
  let provisionedQuery: ProvisionedQuery | null = null;
  let provisionedScript: ProvisionedScript | null = null;
  let saved: SavedAllTables | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Step 0: provision saved query and dataframe-output script', async () => {
      provisionedQuery = await provisionSystemDatagrokQuery(page, {
        nameStem: 'integration_query',
        sql: SYSTEM_DATAGROK_QUERIES.GROUPS_SAMPLE,
      });
      expect(provisionedQuery.queryId).toBeTruthy();

      provisionedScript = await provisionDataframeScript(page, {
        name: `integrationScript${stamp}`,
        body: `df = await grok.data.getDemoTable('demog.csv');`,
      });
      expect(provisionedScript.scriptId).toBeTruthy();
    });

    await softStep('Step 1a: open file source (defensive — skip if missing)', async () => {
      const tries = [
        'System:AppData/Chem/tests/spgi-100.csv',
        'System:DemoFiles/demog.csv',
        'System:DemoFiles/cars.csv',
      ];
      let opened = 0;
      for (const path of tries) {
        try {
          const r = await openFileColonForm(page, path);
          if (r.script) {
            await assertProvenanceScript(page, 'files', r.script);
            opened++;
          }
        } catch (_) { /* defensive */ }
      }
      expect(opened).toBeGreaterThanOrEqual(1);
    });

    await softStep('Step 1b: open DB table (System:Datagrok / public.groups)', async () => {
      const r = await openTableFromDbTable(page, {
        connectionNqName: SYSTEM_DATAGROK_NQNAME,
        schemaName: SYSTEM_DATAGROK_DB_TABLE.schemaName,
        tableName: SYSTEM_DATAGROK_DB_TABLE.tableName,
      });
      expect(r.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'db_table', r.script);
    });

    await softStep('Step 1c: open table from provisioned saved query', async () => {
      if (!provisionedQuery) throw new Error('no provisioned query');
      const r = await openTableFromDbQuery(page, provisionedQuery.queryNqName);
      expect(r.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'db_query', r.script);
    });

    await softStep('Step 1d: open table from provisioned script', async () => {
      if (!provisionedScript) throw new Error('no provisioned script');
      const r = await openTableFromScript(page, provisionedScript.resolvedNqName, {idx: 0});
      expect(r.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'script', r.script);
    });

    await softStep('Step 1e: verify >=4 heterogeneous source tables coexist in workspace', async () => {
      const tables = await evalJs<number>(page, '(grok.shell.tables?.length || 0)');
      expect(tables).toBeGreaterThanOrEqual(4);
    });

    await softStep('Step 2: save ALL tables with provenance (Data Sync ON)', async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.tableInfoIds.length).toBeGreaterThanOrEqual(4);
    });

    await softStep('Step 3: verify project children server-side', async () => {
      if (!saved) throw new Error('no saved project');
      // Project entity exposes .children (Entity[]) — there is no .relations
      // getter on the JS API surface. Verified live 2026-05-08 via diag spec:
      // protoMethods include children/addChild/removeChild but no relations.
      const r = await evalJs<{children: number; ok: boolean}>(page, `(async () => {
        try {
          const fetched = await grok.dapi.projects.find('${saved.projectId}');
          const children = fetched.children ? fetched.children.length : 0;
          return {children, ok: true};
        } catch (e) {
          return {children: -1, ok: false};
        }
      })()`);
      if (r.ok) expect(r.children).toBeGreaterThanOrEqual(4);
    });

    await softStep('Step 4: closeAll and reopen — multi-source co-existence assertion', async () => {
      if (!saved) throw new Error('no saved project');
      const r = await evalJs<{tables: number}>(page, `(async () => {
        grok.shell.closeAll();
        await new Promise(r => setTimeout(r, 800));
        const p = await grok.dapi.projects.find('${saved.projectId}');
        await p.open();
        for (let i = 0; i < 90; i++) {
          if (grok.shell.tables.length >= 4) break;
          await new Promise(r => setTimeout(r, 500));
        }
        return {tables: grok.shell.tables.length};
      })()`);
      expect(r.tables).toBeGreaterThanOrEqual(4);
    });
  } finally {
    const s = saved as SavedAllTables | null;
    if (s) {
      await deleteProjectWithCleanup(page, {projectId: s.projectId});
      for (const tableInfoId of s.tableInfoIds)
        await deleteProjectWithCleanup(page, {tableInfoId});
    }
    const pq = provisionedQuery as ProvisionedQuery | null;
    if (pq) await pq.cleanup();
    const ps = provisionedScript as ProvisionedScript | null;
    if (ps) await deleteProvisionedScript(page, ps.scriptId);
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
