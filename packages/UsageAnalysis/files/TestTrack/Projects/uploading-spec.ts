/* ---
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync]
generated_from: uploading.md (8 cases × 2 sync states)
--- */
// Source-matrix scenario covering uploading.md cases 1, 2, 3, 4, 5, 6, 8, 9
// in both Sync ON and Sync OFF variants. All non-file prerequisites are
// created in-test via helpers/openers.ts — no Samples package or
// env-provisioned DB connection is required.
//
// Sync ON: source `df.tags['.script']` survives the save round-trip so
// reopen re-executes the creation script.
// Sync OFF: `.script` tags are stripped before save, so reopen relies on
// persisted dataframe bytes (snapshot mode).
//
// Cases 4, 5, 6 (Spaces) provision a transient root Space via the shared
// `provisionSpaceFixture` helper, copy demog.csv into it via
// `client.files.write(...)` (so the file is openable via
// `Spaces:<spaceName>/demog.csv`), then run the same save+reopen flow as
// the file/query/derived cases. The Space is released in `finally` after
// reopen — kept alive until then so Sync-ON reopen can re-run
// `OpenFile("Spaces.<spaceName>/demog.csv")` from the persisted .script.
import {test, expect, type Page} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {
  projectsTestOptions,
  evalJs,
  gotoApp,
  setupSession,
  provisionSpaceFixture,
  releaseSpaceFixture,
  SpaceFixture,
} from './_helpers';
import {
  openTableFromFile,
  openTableFromDbQuery,
  openTableFromDbTable,
  provisionSystemDatagrokQuery,
  addAggregateToWorkspace,
  assertProvenanceScript,
  resetShell,
  PROVENANCE_PATTERNS,
  ProvisionedQuery,
  SYSTEM_DATAGROK_DB_TABLE,
  SYSTEM_DATAGROK_NQNAME,
  SYSTEM_DATAGROK_QUERIES,
} from '../helpers/openers';
import {
  saveAllTablesWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
  SavedAllTables,
} from '../helpers/projects';

test.use(projectsTestOptions);

// ---------------------------------------------------------------------------
// Local helpers
// ---------------------------------------------------------------------------

/**
 * Strip the `.script` provenance tag from every open dataframe — emulates
 * the Save Project dialog's "Data Sync OFF" toggle. After strip, reopen
 * relies on persisted dataframe bytes only (snapshot mode).
 */
async function stripProvenance(page: Page): Promise<void> {
  await evalJs(page, `(async () => {
    for (const df of grok.shell.tables) {
      if (df.tags?.has?.('.script')) df.tags.delete('.script');
    }
  })()`);
}

/**
 * Provision a transient root Space with a physical copy of demog.csv. The
 * file lands at `Spaces:<spaceName>/demog.csv` and is openable via the
 * canonical `OpenFile` opener (which writes the same `.script` provenance
 * tag as a System:DemoFiles double-click).
 *
 * Returns either a usable fixture or an env-skip blocker. Callers that hit
 * a blocker should `test.skip(true, reason)` — the Spaces createRootSpace
 * API may legitimately not exist on older platform builds.
 */
async function provisionSpaceWithDemog(
  page: Page, namePrefix: string,
): Promise<{fixture: SpaceFixture} | {blocked: true; reason: string}> {
  const probe = await provisionSpaceFixture(page, {
    namePrefix,
    fileName: 'demog.csv',
  });
  if (probe.blocked || !probe.fixture) {
    if (probe.fixture) await releaseSpaceFixture(page, probe.fixture);
    return {blocked: true, reason: probe.reason};
  }
  return {fixture: probe.fixture};
}

function throwOnStepErrors() {
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
}

// ---------------------------------------------------------------------------
// Case 1 — Files + Files (Link Tables UI delegated to projects-ui-smoke)
// ---------------------------------------------------------------------------

async function runCase1(page: Page, sync: 'on' | 'off') {
  const stamp = Date.now();
  const projectName = `Test_Case1_${sync === 'on' ? 'Sync' : 'NoSync'}_${stamp}`;
  let saved: SavedAllTables | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('open demog.csv and cars.csv', async () => {
      const t1 = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      expect(t1.rowCount).toBeGreaterThan(0);
      const t2 = await openTableFromFile(page, 'System:DemoFiles/cars.csv');
      expect(t2.rowCount).toBeGreaterThan(0);
    });

    if (sync === 'off') await stripProvenance(page);

    await softStep(`save all tables with Data Sync ${sync.toUpperCase()}`, async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.tableInfoIds.length).toBeGreaterThanOrEqual(2);
    });

    await softStep('reopen verifies both tables re-materialize', async () => {
      if (!saved) throw new Error('no saved project');
      const expectedPattern = sync === 'on' ? PROVENANCE_PATTERNS.files : undefined;
      const result = await reopenAndAssertProvenance(page, saved.projectId, expectedPattern);
      expect(result.tablesAfter).toBeGreaterThanOrEqual(2);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
    });
  } finally {
    if (saved) {
      await deleteProjectWithCleanup(page, {projectId: saved.projectId});
      for (const tableInfoId of saved.tableInfoIds)
        await deleteProjectWithCleanup(page, {tableInfoId});
    }
  }

  throwOnStepErrors();
}

test('Projects / Uploading / Case 1: Files + Files (Sync ON)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase1(page, 'on');
});

test('Projects / Uploading / Case 1: Files + Files (Sync OFF)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase1(page, 'off');
});

// ---------------------------------------------------------------------------
// Case 2 — Query + Query (provisioned on System:Datagrok)
// ---------------------------------------------------------------------------

async function runCase2(page: Page, sync: 'on' | 'off') {
  const stamp = Date.now();
  const projectName = `Test_Case2_${sync === 'on' ? 'Sync' : 'NoSync'}_${stamp}`;
  let provisioned: ProvisionedQuery | null = null;
  let saved: SavedAllTables | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('provision query and run twice (two query result tables)', async () => {
      provisioned = await provisionSystemDatagrokQuery(page, {
        nameStem: 'case2_query',
        sql: SYSTEM_DATAGROK_QUERIES.GROUPS_SAMPLE,
      });
      expect(provisioned.queryId).toBeTruthy();
      const t1 = await openTableFromDbQuery(page, provisioned.queryNqName);
      expect(t1.rowCount).toBeGreaterThan(0);
      const t2 = await openTableFromDbQuery(page, provisioned.queryNqName);
      expect(t2.rowCount).toBeGreaterThan(0);
    });

    if (sync === 'off') await stripProvenance(page);

    await softStep(`save all tables with Data Sync ${sync.toUpperCase()}`, async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.tableInfoIds.length).toBeGreaterThanOrEqual(2);
    });

    await softStep('reopen verifies both query results re-materialize', async () => {
      if (!saved) throw new Error('no saved project');
      const expectedPattern = sync === 'on' ? PROVENANCE_PATTERNS.db_query : undefined;
      const result = await reopenAndAssertProvenance(page, saved.projectId, expectedPattern);
      expect(result.tablesAfter).toBeGreaterThanOrEqual(2);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
    });
  } finally {
    if (saved) {
      await deleteProjectWithCleanup(page, {projectId: saved.projectId});
      for (const tableInfoId of saved.tableInfoIds)
        await deleteProjectWithCleanup(page, {tableInfoId});
    }
    if (provisioned) await provisioned.cleanup();
  }

  throwOnStepErrors();
}

test('Projects / Uploading / Case 2: Query + Query (Sync ON)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase2(page, 'on');
});

test('Projects / Uploading / Case 2: Query + Query (Sync OFF)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase2(page, 'off');
});

// ---------------------------------------------------------------------------
// Case 3 — Query + File
// ---------------------------------------------------------------------------

async function runCase3(page: Page, sync: 'on' | 'off') {
  const stamp = Date.now();
  const projectName = `Test_Case3_${sync === 'on' ? 'Sync' : 'NoSync'}_${stamp}`;
  let provisioned: ProvisionedQuery | null = null;
  let saved: SavedAllTables | null = null;
  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('provision query, run it, and open spgi-100 file', async () => {
      provisioned = await provisionSystemDatagrokQuery(page, {
        nameStem: 'case3_query',
        sql: SYSTEM_DATAGROK_QUERIES.GROUPS_SAMPLE,
      });
      expect(provisioned.queryId).toBeTruthy();
      const t1 = await openTableFromDbQuery(page, provisioned.queryNqName);
      expect(t1.rowCount).toBeGreaterThan(0);
      const t2 = await openTableFromFile(page, 'System:AppData/Chem/tests/spgi-100.csv');
      expect(t2.rowCount).toBeGreaterThan(0);
    });

    if (sync === 'off') await stripProvenance(page);

    await softStep(`save all tables with Data Sync ${sync.toUpperCase()}`, async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.tableInfoIds.length).toBeGreaterThanOrEqual(2);
    });

    await softStep('reopen verifies query result + file table re-materialize', async () => {
      if (!saved) throw new Error('no saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThanOrEqual(2);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
    });
  } finally {
    if (saved) {
      await deleteProjectWithCleanup(page, {projectId: saved.projectId});
      for (const tableInfoId of saved.tableInfoIds)
        await deleteProjectWithCleanup(page, {tableInfoId});
    }
    if (provisioned) await provisioned.cleanup();
  }

  throwOnStepErrors();
}

test('Projects / Uploading / Case 3: Query + File (Sync ON)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase3(page, 'on');
});

test('Projects / Uploading / Case 3: Query + File (Sync OFF)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase3(page, 'off');
});

// ---------------------------------------------------------------------------
// Case 4 — Spaces + Spaces (open demog.csv twice from the same Space)
// ---------------------------------------------------------------------------

async function runCase4(page: Page, sync: 'on' | 'off') {
  const stamp = Date.now();
  const projectName = `Test_Case4_${sync === 'on' ? 'Sync' : 'NoSync'}_${stamp}`;
  let saved: SavedAllTables | null = null;
  let fixture: SpaceFixture | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  const provisioning = await provisionSpaceWithDemog(page, `case4-${sync}-${stamp}`);
  if ('blocked' in provisioning) {
    console.warn(`Case 4 env-skip — ${provisioning.reason}`);
    test.skip(true, provisioning.reason);
    return;
  }
  fixture = provisioning.fixture;

  try {
    await softStep('open demog.csv from Space twice', async () => {
      const t1 = await openTableFromFile(page, fixture!.filePath);
      expect(t1.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', t1.script);
      const t2 = await openTableFromFile(page, fixture!.filePath);
      expect(t2.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', t2.script);
    });

    if (sync === 'off') await stripProvenance(page);

    await softStep(`save all tables with Data Sync ${sync.toUpperCase()}`, async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.tableInfoIds.length).toBeGreaterThanOrEqual(2);
    });

    await softStep('reopen verifies both Space tables re-materialize', async () => {
      if (!saved) throw new Error('no saved project');
      const expectedPattern = sync === 'on' ? PROVENANCE_PATTERNS.files : undefined;
      const result = await reopenAndAssertProvenance(page, saved.projectId, expectedPattern);
      expect(result.tablesAfter).toBeGreaterThanOrEqual(2);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
    });
  } finally {
    if (saved) {
      await deleteProjectWithCleanup(page, {projectId: saved.projectId});
      for (const tableInfoId of saved.tableInfoIds)
        await deleteProjectWithCleanup(page, {tableInfoId});
    }
    if (fixture) await releaseSpaceFixture(page, fixture);
  }

  throwOnStepErrors();
}

test('Projects / Uploading / Case 4: Spaces + Spaces (Sync ON)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase4(page, 'on');
});

test('Projects / Uploading / Case 4: Spaces + Spaces (Sync OFF)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase4(page, 'off');
});

// ---------------------------------------------------------------------------
// Case 5 — Spaces + File (demog from Space + spgi-100 from System:AppData)
// ---------------------------------------------------------------------------

async function runCase5(page: Page, sync: 'on' | 'off') {
  const stamp = Date.now();
  const projectName = `Test_Case5_${sync === 'on' ? 'Sync' : 'NoSync'}_${stamp}`;
  let saved: SavedAllTables | null = null;
  let fixture: SpaceFixture | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  const provisioning = await provisionSpaceWithDemog(page, `case5-${sync}-${stamp}`);
  if ('blocked' in provisioning) {
    console.warn(`Case 5 env-skip — ${provisioning.reason}`);
    test.skip(true, provisioning.reason);
    return;
  }
  fixture = provisioning.fixture;

  try {
    await softStep('open demog.csv from Space + spgi-100.csv from Files', async () => {
      const t1 = await openTableFromFile(page, fixture!.filePath);
      expect(t1.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', t1.script);
      const t2 = await openTableFromFile(page, 'System:AppData/Chem/tests/spgi-100.csv');
      expect(t2.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', t2.script);
    });

    if (sync === 'off') await stripProvenance(page);

    await softStep(`save all tables with Data Sync ${sync.toUpperCase()}`, async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.tableInfoIds.length).toBeGreaterThanOrEqual(2);
    });

    await softStep('reopen verifies Space + File tables re-materialize', async () => {
      if (!saved) throw new Error('no saved project');
      // Mixed sources — both happen to be `files` pattern but reopen verifies
      // the active TableView only; pattern check is loose.
      const expectedPattern = sync === 'on' ? PROVENANCE_PATTERNS.files : undefined;
      const result = await reopenAndAssertProvenance(page, saved.projectId, expectedPattern);
      expect(result.tablesAfter).toBeGreaterThanOrEqual(2);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
    });
  } finally {
    if (saved) {
      await deleteProjectWithCleanup(page, {projectId: saved.projectId});
      for (const tableInfoId of saved.tableInfoIds)
        await deleteProjectWithCleanup(page, {tableInfoId});
    }
    if (fixture) await releaseSpaceFixture(page, fixture);
  }

  throwOnStepErrors();
}

test('Projects / Uploading / Case 5: Spaces + File (Sync ON)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase5(page, 'on');
});

test('Projects / Uploading / Case 5: Spaces + File (Sync OFF)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase5(page, 'off');
});

// ---------------------------------------------------------------------------
// Case 6 — Spaces + Query (demog from Space + provisioned System:Datagrok query)
// ---------------------------------------------------------------------------

async function runCase6(page: Page, sync: 'on' | 'off') {
  const stamp = Date.now();
  const projectName = `Test_Case6_${sync === 'on' ? 'Sync' : 'NoSync'}_${stamp}`;
  let saved: SavedAllTables | null = null;
  let fixture: SpaceFixture | null = null;
  let provisioned: ProvisionedQuery | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  const provisioning = await provisionSpaceWithDemog(page, `case6-${sync}-${stamp}`);
  if ('blocked' in provisioning) {
    console.warn(`Case 6 env-skip — ${provisioning.reason}`);
    test.skip(true, provisioning.reason);
    return;
  }
  fixture = provisioning.fixture;

  try {
    await softStep('open demog.csv from Space + run provisioned query', async () => {
      const t1 = await openTableFromFile(page, fixture!.filePath);
      expect(t1.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', t1.script);
      provisioned = await provisionSystemDatagrokQuery(page, {
        nameStem: 'case6_query',
        sql: SYSTEM_DATAGROK_QUERIES.GROUPS_SAMPLE,
      });
      expect(provisioned.queryId).toBeTruthy();
      const t2 = await openTableFromDbQuery(page, provisioned.queryNqName);
      expect(t2.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'db_query', t2.script);
    });

    if (sync === 'off') await stripProvenance(page);

    await softStep(`save all tables with Data Sync ${sync.toUpperCase()}`, async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.tableInfoIds.length).toBeGreaterThanOrEqual(2);
    });

    await softStep('reopen verifies Space + Query tables re-materialize', async () => {
      if (!saved) throw new Error('no saved project');
      // Mixed sources — pattern check skipped; verify multi-table reopen.
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThanOrEqual(2);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
    });
  } finally {
    if (saved) {
      await deleteProjectWithCleanup(page, {projectId: saved.projectId});
      for (const tableInfoId of saved.tableInfoIds)
        await deleteProjectWithCleanup(page, {tableInfoId});
    }
    if (provisioned) await provisioned.cleanup();
    if (fixture) await releaseSpaceFixture(page, fixture);
  }

  throwOnStepErrors();
}

test('Projects / Uploading / Case 6: Spaces + Query (Sync ON)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase6(page, 'on');
});

test('Projects / Uploading / Case 6: Spaces + Query (Sync OFF)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase6(page, 'off');
});

// ---------------------------------------------------------------------------
// Case 8 — Files + Pivot Table (Add to workspace)
// ---------------------------------------------------------------------------

async function runCase8(page: Page, sync: 'on' | 'off') {
  const stamp = Date.now();
  const projectName = `Test_Case8_${sync === 'on' ? 'Sync' : 'NoSync'}_${stamp}`;
  let saved: SavedAllTables | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('open demog.csv as base + Pivot Table → Add', async () => {
      const base = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      expect(base.rowCount).toBeGreaterThan(0);
      const derived = await addAggregateToWorkspace(page, {via: 'pivot-viewer'});
      expect(derived.rowCount).toBeGreaterThan(0);
    });

    if (sync === 'off') await stripProvenance(page);

    await softStep(`save all tables with Data Sync ${sync.toUpperCase()}`, async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.tableInfoIds.length).toBeGreaterThanOrEqual(2);
    });

    await softStep('reopen verifies base + pivot derived re-materialize', async () => {
      if (!saved) throw new Error('no saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThanOrEqual(2);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
    });
  } finally {
    if (saved) {
      await deleteProjectWithCleanup(page, {projectId: saved.projectId});
      for (const tableInfoId of saved.tableInfoIds)
        await deleteProjectWithCleanup(page, {tableInfoId});
    }
  }

  throwOnStepErrors();
}

test('Projects / Uploading / Case 8: Files + Pivot Table > Add (Sync ON)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase8(page, 'on');
});

test('Projects / Uploading / Case 8: Files + Pivot Table > Add (Sync OFF)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase8(page, 'off');
});

// ---------------------------------------------------------------------------
// Case 9 — DB table (System:Datagrok / public.groups) + Aggregate Rows
// ---------------------------------------------------------------------------

async function runCase9(page: Page, sync: 'on' | 'off') {
  const stamp = Date.now();
  const projectName = `Test_Case9_${sync === 'on' ? 'Sync' : 'NoSync'}_${stamp}`;
  let saved: SavedAllTables | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('open public.groups via DbQuery + Aggregate Rows → Add', async () => {
      const base = await openTableFromDbTable(page, {
        connectionNqName: SYSTEM_DATAGROK_NQNAME,
        schemaName: SYSTEM_DATAGROK_DB_TABLE.schemaName,
        tableName: SYSTEM_DATAGROK_DB_TABLE.tableName,
      });
      expect(base.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'db_table', base.script);
      const derived = await addAggregateToWorkspace(page, {via: 'menu'});
      expect(derived.rowCount).toBeGreaterThan(0);
    });

    if (sync === 'off') await stripProvenance(page);

    await softStep(`save all tables with Data Sync ${sync.toUpperCase()}`, async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.tableInfoIds.length).toBeGreaterThanOrEqual(2);
    });

    await softStep('reopen verifies base + aggregate derived re-materialize', async () => {
      if (!saved) throw new Error('no saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThanOrEqual(2);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
    });
  } finally {
    if (saved) {
      await deleteProjectWithCleanup(page, {projectId: saved.projectId});
      for (const tableInfoId of saved.tableInfoIds)
        await deleteProjectWithCleanup(page, {tableInfoId});
    }
  }

  throwOnStepErrors();
}

test('Projects / Uploading / Case 9: DB + Aggregate Rows > Add (Sync ON)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase9(page, 'on');
});

test('Projects / Uploading / Case 9: DB + Aggregate Rows > Add (Sync OFF)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await runCase9(page, 'off');
});
