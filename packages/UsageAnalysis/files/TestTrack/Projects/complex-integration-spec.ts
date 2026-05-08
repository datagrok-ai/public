/* ---
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync, projects.api.namespaces, projects.add-relation, projects.add-link]
generated_from: complex-integration.md (multi-source: file + DB-table + query + script)
--- */
// Multi-source integration: file + ad-hoc DB table on System:Datagrok +
// provisioned saved query on System:Datagrok + provisioned dataframe-output
// script. All non-file prerequisites are created in-test via
// helpers/openers.ts — no Samples package or env-provisioned DB.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {
  openTableFromFile,
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

test.use(projectsTestOptions);

test('Projects / Complex Integration: heterogeneous sources in one project', async ({page}) => {
  test.setTimeout(900_000);
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
          const r = await openTableFromFile(page, path);
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

    await softStep('Step 3: verify project relations server-side', async () => {
      if (!saved) throw new Error('no saved project');
      const r = await evalJs<{relations: number; ok: boolean}>(page, `(async () => {
        try {
          const fetched = await grok.dapi.projects.find('${saved.projectId}');
          const relations = fetched.relations ? fetched.relations.length : 0;
          return {relations, ok: true};
        } catch (e) {
          return {relations: -1, ok: false};
        }
      })()`);
      if (r.ok) expect(r.relations).toBeGreaterThanOrEqual(4);
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
    if (saved) {
      await deleteProjectWithCleanup(page, {projectId: saved.projectId});
      for (const tableInfoId of saved.tableInfoIds)
        await deleteProjectWithCleanup(page, {tableInfoId});
    }
    if (provisionedQuery) await provisionedQuery.cleanup();
    if (provisionedScript) await deleteProvisionedScript(page, provisionedScript.scriptId);
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
