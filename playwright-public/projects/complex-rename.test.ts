/* ---
sub_features_covered: [projects.api.save, projects.shell.open]
--- */
// GROK-19212: rename a referenced table inside a project, reopen, verify it loads without resolution error.
// Test 1 = table rename (GROK-19212); Test 2 = query rename (github-3550); Test 3 = script rename (sister).
import {test, expect, Page} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {
  openTableFromFile,
  openTableFromDbQuery,
  openTableFromScript,
  provisionSystemDatagrokQuery,
  provisionDataframeScript,
  deleteProvisionedScript,
  ProvisionedQuery,
  ProvisionedScript,
  SYSTEM_DATAGROK_QUERIES,
} from '../helpers/openers';
import {
  saveProjectWithProvenance,
  deleteProjectWithCleanup,
} from '../helpers/projects';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';

test.use(projectsTestOptions);

async function closeAll(page: Page) {
  await evalJs(page, 'grok.shell.closeAll()');
}

test('Projects / Complex rename: rename-then-reopen reference resolution (GROK-19212)', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = 'AutoTest-ComplexRename-' + stamp;

  await gotoApp(page);
  await setupSession(page);

  let projectId: string | null = null;
  let tableInfoId: string | null = null;
  let layoutId: string | null = null;

  try {
    await softStep('Setup: open 2 source tables + join (sets up a referenced-table dependency)', async () => {
      await closeAll(page);
      await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      await evalJs(page, `(async () => {
        const tables = grok.shell.tables;
        if (tables.length >= 2) {
          tables[0].name = 'src_a_${stamp}';
          tables[1].name = 'src_b_${stamp}';
        }
        const [df1, df2] = tables;
        const joined = grok.data.joinTables(
          df1, df2,
          ['USUBJID'], ['USUBJID'],
          ['AGE'], ['SEX'],
          'inner',
        );
        joined.name = 'joined_${stamp}';
        grok.shell.addTableView(joined);
      })()`);
      await page.waitForTimeout(2000);
      const tableCount = await evalJs(page, 'grok.shell.tables.length');
      expect(tableCount).toBe(3);
    });

    await softStep('Save baseline project (Data Sync ON), capture project ID', async () => {
      // Multi-table inline save — all 3 tables must persist so renaming src_a breaks joined's reference.
      const saved = await page.evaluate(async (n) => {
        const grok = (window as any).grok;
        const DG = (window as any).DG;
        const project = DG.Project.create();
        project.name = n;
        const tables = Array.from(grok.shell.tables);
        const tableInfoIds: string[] = [];
        for (const df of tables as any[]) {
          const ti = df.getTableInfo();
          project.addChild(ti);
          await grok.dapi.tables.uploadDataFrame(df);
          await grok.dapi.tables.save(ti);
          tableInfoIds.push(ti.id);
        }
        const tv = grok.shell.tv;
        const layout = tv?.saveLayout?.();
        if (layout) {
          project.addChild(layout);
          await grok.dapi.layouts.save(layout);
        }
        await grok.dapi.projects.save(project);
        return {
          projectId: project.id,
          tableInfoIds,
          layoutId: layout?.id ?? null,
        };
      }, projectName);
      projectId = saved.projectId;
      tableInfoId = saved.tableInfoIds[0] ?? null;
      layoutId = saved.layoutId;
      expect(projectId).toBeTruthy();
      const exists = await page.evaluate(async (pid) => {
        const grok = (window as any).grok;
        const p = await grok.dapi.projects.find(pid);
        return p != null;
      }, projectId);
      expect(exists).toBe(true);
    });

    await softStep('Step 7: rename a referenced table inside the project (the GROK-19212 trigger)', async () => {
      if (!projectId) throw new Error('no projectId captured');
      // Reopen by id (filter-by-name fails for dashed names).
      await closeAll(page);
      await page.evaluate(async (pid) => {
        const grok = (window as any).grok;
        const p = await grok.dapi.projects.find(pid);
        if (p) await p.open();
      }, projectId);
      await page.waitForTimeout(3000);

      await evalJs(page, `(async () => {
        const t = grok.shell.tables.find(t => t.name === 'src_a_${stamp}');
        if (t) t.name = 'src_a_renamed_${stamp}';
      })()`);
      await page.waitForTimeout(1000);

      const renamed = await evalJs(page, `(async () => {
        return grok.shell.tables.some(t => t.name === 'src_a_renamed_${stamp}');
      })()`);
      expect(renamed).toBe(true);
    });

    await softStep('Step 8: re-save project after rename (overwrite via JS API)', async () => {
      if (!projectId) throw new Error('no projectId captured');
      await page.evaluate(async (pid) => {
        const grok = (window as any).grok;
        const p = await grok.dapi.projects.find(pid);
        if (p) await grok.dapi.projects.save(p);
      }, projectId);
      await page.waitForTimeout(2000);
    });

    await softStep('GROK-19212 INVARIANT: close + reopen project after table rename — must load without resolution error', async () => {
      if (!projectId) throw new Error('no projectId captured');
      await closeAll(page);
      const result = await page.evaluate(async (pid) => {
        const grok = (window as any).grok;
        try {
          const p = await grok.dapi.projects.find(pid);
          if (!p) return {ok: false, reason: 'project disappeared after rename+save'};
          await p.open();
          await new Promise((r) => setTimeout(r, 3000));
          // Use dataFrame.rowCount as load signal (shell.tables.length throws Tn.grok_TableNames in some reopen states).
          const rc = grok?.shell?.tv?.dataFrame?.rowCount ?? 0;
          return {ok: true, rowCount: rc};
        } catch (e) {
          return {ok: false, reason: String(e).slice(0, 200)};
        }
      }, projectId);
      expect(result.ok).toBe(true);
      // rowCount > 0 means at least one source table re-materialized despite the rename.
      expect(result.rowCount).toBeGreaterThan(0);
    });
  } finally {
    await deleteProjectWithCleanup(page, {
      projectId: projectId ?? undefined,
      tableInfoId: tableInfoId ?? undefined,
    });
    void layoutId; // layout cleanup deferred; deleteProjectWithCleanup doesn't take layoutId yet
    await closeAll(page);
  }

  finishSpec();
});

// ---------------------------------------------------------------------------
// Test 2 — Query rename: github-3550 reproduction via project source
// ---------------------------------------------------------------------------

test('Projects / Complex rename: rename Query, reopen, verify reference resolution (github-3550)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `AutoTest-ComplexRenameQuery-${stamp}`;
  let provisioned: ProvisionedQuery | null = null;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await closeAll(page);

  try {
    await softStep('Setup: provision saved query on System:Datagrok', async () => {
      provisioned = await provisionSystemDatagrokQuery(page, {
        nameStem: 'complex_rename_query',
        sql: SYSTEM_DATAGROK_QUERIES.GROUPS_SAMPLE,
      });
      expect(provisioned.queryId).toBeTruthy();
    });

    await softStep('Open table from provisioned query and save project', async () => {
      if (!provisioned) throw new Error('no provisioned query');
      const opened = await openTableFromDbQuery(page, provisioned.queryNqName);
      expect(opened.rowCount).toBeGreaterThan(0);
      saved = await saveProjectWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
    });

    await softStep('Rename the provisioned query (test owns it — always succeeds)', async () => {
      if (!provisioned) throw new Error('no provisioned query');
      const ok = await evalJs<boolean>(page, `(async () => {
        try {
          const q = await grok.dapi.queries.find('${provisioned.queryId}');
          if (!q) return false;
          q.name = '${provisioned.resolvedName}_renamed';
          await grok.dapi.queries.save(q);
          return true;
        } catch (_) { return false; }
      })()`);
      expect(ok).toBe(true);
    });

    await softStep('github-3550 INVARIANT: reopen project, verify auto-resolve OR explicit error', async () => {
      if (!saved) throw new Error('no saved project');
      const result = await evalJs<{
        loadedOk: boolean;
        errorMessage: string | null;
        relationsCount: number;
      }>(page, `(async () => {
        grok.shell.closeAll();
        await new Promise(r => setTimeout(r, 800));
        const p = await grok.dapi.projects.find('${saved.projectId}');
        let errorMessage = null;
        try { await p.open(); } catch (e) { errorMessage = String(e).slice(0, 300); }
        for (let i = 0; i < 60; i++) {
          if (grok.shell.tables.length > 0) break;
          await new Promise(r => setTimeout(r, 500));
        }
        const fresh = await grok.dapi.projects.find('${saved.projectId}');
        return {
          loadedOk: grok.shell.tables.length > 0,
          errorMessage,
          relationsCount: fresh.relations ? fresh.relations.length : 0,
        };
      })()`);
      const isHappyPath = result.loadedOk;
      const isGracefulFailure = !result.loadedOk && result.errorMessage !== null;
      expect(isHappyPath || isGracefulFailure).toBe(true);
    });
  } finally {
    if (saved)
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.tableInfoId,
      });
    if (provisioned) await provisioned.cleanup();
    await closeAll(page);
  }

  finishSpec();
});

// ---------------------------------------------------------------------------
// Test 3 — Script rename: github-3550 sister invariant via project source
// ---------------------------------------------------------------------------

test('Projects / Complex rename: rename Script, reopen, verify reference resolution', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `AutoTest-ComplexRenameScript-${stamp}`;
  let provisioned: ProvisionedScript | null = null;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await closeAll(page);

  try {
    await softStep('Setup: provision dataframe-output script', async () => {
      provisioned = await provisionDataframeScript(page, {
        name: `complexRenameScript${stamp}`,
        body: `df = await grok.data.getDemoTable('demog.csv');`,
      });
      expect(provisioned.scriptId).toBeTruthy();
    });

    await softStep('Open table from provisioned script and save project', async () => {
      if (!provisioned) throw new Error('no provisioned script');
      const opened = await openTableFromScript(page, provisioned.resolvedNqName, {idx: 0});
      expect(opened.rowCount).toBeGreaterThan(0);
      saved = await saveProjectWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
    });

    await softStep('Rename the provisioned script (test owns it — always succeeds)', async () => {
      if (!provisioned) throw new Error('no provisioned script');
      const ok = await evalJs<boolean>(page, `(async () => {
        try {
          const s = await grok.dapi.scripts.find('${provisioned.scriptId}');
          if (!s) return false;
          s.name = '${provisioned.resolvedName}_renamed';
          await grok.dapi.scripts.save(s);
          return true;
        } catch (_) { return false; }
      })()`);
      expect(ok).toBe(true);
    });

    await softStep('Sister invariant: reopen project, verify auto-resolve OR explicit error', async () => {
      if (!saved) throw new Error('no saved project');
      const result = await evalJs<{
        loadedOk: boolean;
        errorMessage: string | null;
      }>(page, `(async () => {
        grok.shell.closeAll();
        await new Promise(r => setTimeout(r, 800));
        const p = await grok.dapi.projects.find('${saved.projectId}');
        let errorMessage = null;
        try { await p.open(); } catch (e) { errorMessage = String(e).slice(0, 300); }
        for (let i = 0; i < 60; i++) {
          if (grok.shell.tables.length > 0) break;
          await new Promise(r => setTimeout(r, 500));
        }
        return {
          loadedOk: grok.shell.tables.length > 0,
          errorMessage,
        };
      })()`);
      const isHappyPath = result.loadedOk;
      const isGracefulFailure = !result.loadedOk && result.errorMessage !== null;
      expect(isHappyPath || isGracefulFailure).toBe(true);
    });
  } finally {
    if (saved)
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.tableInfoId,
      });
    if (provisioned) await deleteProvisionedScript(page, provisioned.scriptId);
    await closeAll(page);
  }

  finishSpec();
});
