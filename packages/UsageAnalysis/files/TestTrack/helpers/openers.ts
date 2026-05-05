/**
 * Playwright openers for table sources (file, DB query, ad-hoc DB table, script).
 *
 * Namespace per helpers-candidates.yaml: `helpers.playwright.openers.*`.
 * Imported as: `import * as openers from '../helpers/openers';` (or named).
 *
 * Authoring cycle: Phase B 2026-05-05 — first openers module.
 *
 * ----------------------------------------------------------------------------
 * Why this module exists
 * ----------------------------------------------------------------------------
 * Prior lifecycle specs opened tables via `grok.dapi.files.readCsv(path)` /
 * `grok.functions.eval('Postgres:NorthwindTest:select_top_n(...)')` /
 * `grok.functions.eval('Samples:Cars()')`. None of those paths write the
 * `df.tags['.script']` provenance tag that Datagrok needs to re-materialize
 * the table on project reopen with Data Sync ON. As a result, save+reopen
 * silently degraded to snapshot-only (or skipped entirely when the function
 * name was wrong).
 *
 * The canonical "recorder" path for opening a source-bound table is:
 *
 *     DG.Func.find({name})[0]
 *       .prepare(args)
 *       .call(undefined, undefined, {processed: false});
 *
 * The `{processed: false}` options bag is the key: it engages Datagrok's
 * function-call recorder, which writes:
 *   - `df.tags['.script']`   — single-line creation script
 *   - `df.tags['.history']`  — full execution log
 *   - `df.tags['.VariableName']` — assigned shell variable
 *
 * On project save (canonical `uploadProject` pattern, see
 * `helpers/projects.ts:saveProjectWithProvenance`), `.script` is persisted
 * server-side. On reopen, Datagrok re-executes the script — re-reading the
 * file / re-running the query / re-running the script, depending on the
 * source class.
 *
 * Verified live against dev.datagrok.ai 2026-05-05 (qa-pw). MCP capture
 * findings + per-source `.script` regex patterns: see
 * `.claude/diagnostics/mcp-capture-{files,db,scripts}.md`.
 *
 * Reference: `public/packages/UITests/src/gui/gui-utils.ts:100-111` — the
 * platform's own `uploadProject` helper, which is the second half of the
 * lifecycle (save with provenance attached).
 */

import {Page} from '@playwright/test';

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

export interface OpenedTable {
  /** Active TableView's dataFrame name (often lowercased). */
  name: string;
  /** Row count of the opened df. */
  rowCount: number;
  /** Column count of the opened df. */
  colCount: number;
  /** TableInfo id (df.getTableInfo().id) — useful for downstream cleanup. */
  tableInfoId: string;
  /** `df.tags['.script']` — the provenance creation script. */
  script: string;
}

export type SourceClass = 'files' | 'db_query' | 'db_table' | 'script' | 'derived';

/**
 * Per-source-class regex patterns for the `df.tags['.script']` tag.
 * Used by Gate E-PROV-01 to verify that a Step 1 actually wired up
 * provenance and didn't silently degrade to snapshot.
 */
export const PROVENANCE_PATTERNS: Record<SourceClass, RegExp> = {
  files: /^\w+\s*=\s*OpenFile\("[^"]+"\)/,
  db_query: /^\w+\s*=\s*[\w:]+\(\)/,
  db_table: /^\w+\s*=\s*DbQuery\([\w:]+,\s*"\w+"/,
  script: /^\w+\s*=\s*[\w:]+\([\w\s,'"=.-]*\)/,
  // Derived (Aggregate Rows / Pivot Table → Add to workspace) — references
  // the base table by quoted name, so on reopen the base must be present
  // first. Server resolves dependency order across project children.
  derived: /^\w+\s*=\s*Aggregate\("[\w\s]+",/,
};

// ---------------------------------------------------------------------------
// 0. resetShell — call at the start of each test in a multi-test file.
// ---------------------------------------------------------------------------

/**
 * Close all open views, dialogs, toasts, menus. Use at the start of each
 * test in a multi-test file to prevent cross-test state leakage in a
 * single Playwright worker.
 *
 * Why this is needed: `gotoApp` calls `page.goto(BASE_URL)` which
 * navigates the page, but the Datagrok SPA caches some state across
 * navigations (open dialogs, transient menus). When test #1 leaves an
 * orphan dialog up, test #2's `openTableFromFile` can fail with
 * "did not produce a TableView" because click events get intercepted.
 *
 * Verified empirical pattern 2026-05-05 (Phase B re-runs): solo runs
 * pass, multi-test runs fail without this helper.
 *
 * Idempotent — safe to call when shell is already clean.
 */
export async function resetShell(page: Page): Promise<void> {
  await page.evaluate(() => {
    const grok = (window as any).grok;
    try { grok?.shell?.closeAll?.(); } catch (_) {}
    document.querySelectorAll('.d4-dialog, .d4-toast, .d4-balloon, .d4-menu-popup')
      .forEach((e) => { try { (e as any).remove?.(); } catch (_) {} });
  });
  await page.waitForTimeout(300);
}

// ---------------------------------------------------------------------------
// 1. openTableFromFile — Browse → Files → ... → double-click semantics.
// ---------------------------------------------------------------------------

/**
 * Open a CSV / Excel / etc. file from a Datagrok file share via the
 * canonical `OpenFile` function.
 *
 * Equivalent to UI `Browse → Files → <namespace> → ... → <file> → double-click`,
 * but does not depend on file-tree expansion / DOM selectors. Both paths
 * write the same `.script: 'Var = OpenFile("<fullPath>")'` provenance tag.
 *
 * Server endpoint hit:
 *   `GET /api/connectors/connections/<namespace>/file/<relPath>?format=d42`
 *   (namespace path uses `.` not `:` — `System.DemoFiles` not `System:DemoFiles`).
 *
 * @param page - Playwright Page (login complete).
 * @param fullPath - File path with `:` namespace separator,
 *   e.g. `'System:DemoFiles/demog.csv'`, `'System:AppData/Chem/tests/spgi-100.csv'`.
 * @returns metadata about the opened table including the `.script` tag.
 * @throws if `OpenFile` is not registered or did not produce a TableView.
 */
export async function openTableFromFile(
  page: Page,
  fullPath: string,
): Promise<OpenedTable> {
  return await page.evaluate(async (p) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const fns = DG.Func.find({name: 'OpenFile'});
    if (!fns?.length) throw new Error('OpenFile function not registered');
    await fns[0].prepare({fullPath: p}).call(undefined, undefined, {processed: false});
    // Settle: poll for shell.tv.dataFrame up to 8s — covers dev cold-start
    // under multi-test load. Earlier fixed 800ms / 1200ms / 2000ms waits
    // were brittle.
    let df: any = null;
    for (let i = 0; i < 16; i++) {
      df = grok.shell.tv?.dataFrame;
      if (df) break;
      await new Promise((r) => setTimeout(r, 500));
    }
    if (!df) throw new Error(`OpenFile("${p}") did not produce a TableView (8s settle)`);
    const ti = df.getTableInfo?.();
    return {
      name: df.name,
      rowCount: df.rowCount,
      colCount: df.columns.length,
      tableInfoId: ti?.id,
      script: df.tags?.get?.('.script') ?? '',
    };
  }, fullPath);
}

// ---------------------------------------------------------------------------
// 2. openTableFromDbQuery — Browse → Platform → Functions → Queries → call.
// ---------------------------------------------------------------------------

/**
 * Run a saved query as a function and open its result as a TableView.
 *
 * Equivalent to UI `Browse → Platform → Functions → Queries → <query> →
 * double-click` for queries with no required parameters (or with values
 * supplied via `args`). `.script` tag is written as
 * `Var = <Namespace>:<QueryName>([args])`.
 *
 * Result df also carries DB-binding tags:
 *   `.data-connection-id`, `.data-connection-nqName`,
 *   `.DataQuery.query.finished`.
 *
 * @param page - Playwright Page.
 * @param queryNqName - Fully qualified query name including namespace,
 *   e.g. `'Samples:PostgresProducts'` or `'Samples:PostgresCustomers'`.
 *   Bare name (no `:`) is also accepted — first `DG.Func.find({name})` match wins.
 * @param args - Query input parameters (e.g. `{country: 'Germany'}`).
 *   Empty for parameter-less queries.
 * @returns metadata about the opened table including the `.script` tag.
 */
export async function openTableFromDbQuery(
  page: Page,
  queryNqName: string,
  args: Record<string, unknown> = {},
): Promise<OpenedTable> {
  return await page.evaluate(async ({nq, a}) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const [namespace, name] = nq.includes(':') ? nq.split(':') : [null, nq];
    const fns = namespace
      ? DG.Func.find({namespace, name})
      : DG.Func.find({name});
    if (!fns?.length) throw new Error(`Query function not registered: ${nq}`);
    await fns[0].prepare(a).call(undefined, undefined, {processed: false});
    await new Promise((r) => setTimeout(r, 1000));
    const df = grok.shell.tv?.dataFrame;
    if (!df) throw new Error(`Query ${nq} did not produce a TableView`);
    const ti = df.getTableInfo?.();
    return {
      name: df.name,
      rowCount: df.rowCount,
      colCount: df.columns.length,
      tableInfoId: ti?.id,
      script: df.tags?.get?.('.script') ?? '',
    };
  }, {nq: queryNqName, a: args});
}

// ---------------------------------------------------------------------------
// 3. openTableFromDbTable — Browse → Databases → ... → table → double-click.
// ---------------------------------------------------------------------------

/**
 * Open an ad-hoc DB table via `core:DbQuery` — the visual-query-builder
 * function that the UI uses when you double-click a schema-tree node.
 *
 * Equivalent to UI `Browse → Databases → <connection> → Schemas → <schema>
 * → <table> → double-click`. `.script` tag is written as
 * `Var = DbQuery(<conn-nqName>, "<table>", schemaName = "<schema>")`.
 *
 * Pass `limit: 100` to reproduce the right-click "Get Top 100" path.
 * Omit `limit` (or pass 0) to reproduce "Get All". Both reduce to the same
 * `DbQuery` call with different `limit` arg, so this single helper covers
 * both UI menu items.
 *
 * @param page - Playwright Page.
 * @param options.connectionNqName - Fully qualified connection name (preferred),
 *   e.g. `'Samples:PostgresNorthwind'`.
 * @param options.connectionId - Connection UUID (alternative to nqName).
 * @param options.schemaName - Database schema, e.g. `'public'`.
 * @param options.tableName - Table name within the schema, e.g. `'products'`.
 * @param options.limit - Optional LIMIT clause; 100 = "Get Top 100";
 *   omitted = "Get All".
 * @returns metadata about the opened table including the `.script` tag.
 */
export async function openTableFromDbTable(
  page: Page,
  options: {
    connectionNqName?: string;
    connectionId?: string;
    schemaName: string;
    tableName: string;
    limit?: number;
  },
): Promise<OpenedTable> {
  if (!options.connectionId && !options.connectionNqName)
    throw new Error('openTableFromDbTable: connectionId or connectionNqName required');
  return await page.evaluate(async (o) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    let conn: any = null;
    if (o.connectionId)
      conn = await grok.dapi.connections.find(o.connectionId);
    else if (o.connectionNqName) {
      const parts = o.connectionNqName.split(':');
      const bare = parts.length > 1 ? parts[1] : parts[0];
      // Note: dapi.connections.filter `name = "X"` matches against the
      // fully-qualified server name and returns 0 for bare names. The
      // canonical lookup field on dev is `shortName` — verified live
      // 2026-05-05 (Phase B re-run). The same applies to dapi.queries.
      const list = await grok.dapi.connections.filter(`shortName = "${bare}"`).list();
      conn = list.find((c: any) => c.nqName === o.connectionNqName) || list[0];
    }
    if (!conn) throw new Error(`Connection not found: ${o.connectionNqName ?? o.connectionId}`);
    const fns = DG.Func.find({name: 'DbQuery'});
    if (!fns?.length) throw new Error('core:DbQuery function not registered');
    const args: Record<string, unknown> = {
      conn,
      schemaName: o.schemaName,
      tableName: o.tableName,
    };
    if (o.limit && o.limit > 0) args.limit = o.limit;
    await fns[0].prepare(args).call(undefined, undefined, {processed: false});
    await new Promise((r) => setTimeout(r, 1500));
    const df = grok.shell.tv?.dataFrame;
    if (!df)
      throw new Error(`DbQuery on ${o.schemaName}.${o.tableName} did not produce a TableView`);
    const ti = df.getTableInfo?.();
    return {
      name: df.name,
      rowCount: df.rowCount,
      colCount: df.columns.length,
      tableInfoId: ti?.id,
      script: df.tags?.get?.('.script') ?? '',
    };
  }, options);
}

// ---------------------------------------------------------------------------
// 4. openTableFromScript — run a script with output: dataframe.
// ---------------------------------------------------------------------------

/**
 * Run a script registered as a `DG.Func` with `output: dataframe`.
 *
 * Equivalent to UI `Browse → Platform → Functions → Scripts → <script> →
 * right-click → Run... → fill inputs → OK`. `.script` tag is written as
 * `Var = <Namespace>:<ScriptName>(<args>)`.
 *
 * For tests that need a deterministic table-output script, provision one
 * inline first via {@link provisionDataframeScript}, then pass its name here.
 *
 * Note: the script must produce a DataFrame as one of its outputs. Scalar-
 * output scripts (e.g. `Samples:Cars` which returns `mpgPerCylinder, acc`)
 * will throw "did not produce a TableView" — use a different script.
 *
 * @param page - Playwright Page.
 * @param scriptName - Fully qualified script name (`'QaPw:MyScript'`) or
 *   bare name (first `DG.Func.find({name})` match wins).
 * @param inputs - Script input arguments. Defaults to `{}` — provide all
 *   required inputs that lack defaults.
 * @returns metadata about the opened table including the `.script` tag.
 */
export async function openTableFromScript(
  page: Page,
  scriptName: string,
  inputs: Record<string, unknown> = {},
): Promise<OpenedTable> {
  return await page.evaluate(async ({n, i}) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const [namespace, name] = n.includes(':') ? n.split(':') : [null, n];
    const fns = namespace
      ? DG.Func.find({namespace, name})
      : DG.Func.find({name});
    if (!fns?.length) throw new Error(`Script not registered as DG.Func: ${n}`);
    await fns[0].prepare(i).call(undefined, undefined, {processed: false});
    await new Promise((r) => setTimeout(r, 1500));
    const df = grok.shell.tv?.dataFrame;
    if (!df)
      throw new Error(`Script ${n} did not produce a TableView (scalar-only output?)`);
    const ti = df.getTableInfo?.();
    return {
      name: df.name,
      rowCount: df.rowCount,
      colCount: df.columns.length,
      tableInfoId: ti?.id,
      script: df.tags?.get?.('.script') ?? '',
    };
  }, {n: scriptName, i: inputs});
}

// ---------------------------------------------------------------------------
// 5. provisionDataframeScript — create a JS script with output: dataframe.
// ---------------------------------------------------------------------------

export interface ProvisionedScript {
  /** Server-assigned id of the saved Script entity (for cleanup). */
  scriptId: string;
  /** Name as recorded on the server — may be PascalCased server-side. */
  resolvedName: string;
  /** Fully qualified name including namespace (login for personal scripts). */
  resolvedNqName: string;
}

/**
 * Provision a Datagrok JS script that returns a DataFrame, registered as
 * a `DG.Func`. Used by lifecycle specs that need a deterministic
 * df-output script (Samples has none — Cars is scalar-only, Python/R
 * samples all require an input df).
 *
 * Pattern reference: `scripts-layout.test (1).ts:24-29` — the
 * `test_Layout` fixture wraps `grok.data.getDemoTable('cars.csv')`.
 *
 * **Caller is responsible for cleanup** — call `dapi.scripts.delete(scriptId)`
 * in a `finally` block. Does NOT auto-delete prior runs to avoid clobbering
 * other concurrent tests.
 *
 * @param page - Playwright Page.
 * @param options.name - Script name (will be PascalCased server-side; spec
 *   should use a fresh `Date.now()` suffix to avoid collisions).
 * @param options.body - Script body — should assign to the output variable.
 *   Default body wraps `grok.data.getDemoTable('demog.csv')`.
 * @param options.inputs - Optional `#input:` lines (default: `int idx=0`).
 * @returns metadata for downstream `openTableFromScript` + cleanup.
 */
export async function provisionDataframeScript(
  page: Page,
  options: {
    name: string;
    body?: string;
    inputs?: string[];
  },
): Promise<ProvisionedScript> {
  const inputs = options.inputs ?? ['int idx=0'];
  const body = options.body ?? `df = await grok.data.getDemoTable('demog.csv');`;
  const inputBlock = inputs.map((i) => `//input: ${i}`).join('\n');
  const content =
    `//name: ${options.name}\n` +
    `//language: javascript\n` +
    `${inputBlock}\n` +
    `//output: dataframe df\n` +
    body;

  return await page.evaluate(async ({n, c}) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const script = DG.Script.create(c);
    const saved = await grok.dapi.scripts.save(script);
    // Wait for DG.Func registration to propagate.
    let fn: any = null;
    for (let i = 0; i < 12; i++) {
      fn = DG.Func.find({name: n})?.[0];
      if (!fn) {
        const cap = n.charAt(0).toUpperCase() + n.slice(1);
        fn = DG.Func.find({name: cap})?.[0];
      }
      if (fn) break;
      await new Promise((r) => setTimeout(r, 300));
    }
    if (!fn) throw new Error(`Provisioned script not registered as DG.Func: ${n}`);
    return {
      scriptId: saved?.id,
      resolvedName: fn.name,
      resolvedNqName: fn.nqName,
    };
  }, {n: options.name, c: content});
}

/** Delete a provisioned script by id (best-effort, swallows errors). */
export async function deleteProvisionedScript(
  page: Page,
  scriptId: string,
): Promise<void> {
  await page.evaluate(async (id) => {
    try {
      const grok = (window as any).grok;
      const s = await grok.dapi.scripts.find(id);
      if (s) await grok.dapi.scripts.delete(s);
    } catch (_) { /* best effort */ }
  }, scriptId).catch(() => {});
}

// ---------------------------------------------------------------------------
// 6. addAggregateToWorkspace — derived source class (Aggregate Rows / Pivot Table).
// ---------------------------------------------------------------------------

/**
 * Trigger an Aggregate-Rows or Pivot-Table panel and click ADD to derive
 * a new workspace table from the active TableView.
 *
 * Both UI gateways (Top menu Data → Aggregate Rows, OR toolbox →
 * Pivot Table viewer) end up calling `core:Aggregate` with auto-picked
 * defaults — verified live 2026-05-05. The resulting df carries
 * `.script: 'Var = Aggregate("<baseName>", fields=[...], pivots=[...], ...)'`.
 *
 * Caller must have an active TableView (open a base via
 * {@link openTableFromFile} / {@link openTableFromDbQuery} / etc. first).
 *
 * Lifecycle note: the derived `.script` references the base table BY NAME
 * (not by Datagrok variable), so on project reopen the base must
 * re-materialize first. Datagrok resolves order automatically when both
 * are children of the same project.
 *
 * @param page - Playwright Page (active TableView present).
 * @param options.via - `'menu'` (default) → Top menu Data → Aggregate Rows.
 *   `'pivot-viewer'` → toolbox → Pivot Table viewer. Both lead to the same
 *   ADD button and same `core:Aggregate` provenance — the option exists
 *   because UI-coverage scenarios distinguish the two paths (Cases 8 vs 9
 *   in `Projects/uploading.md`).
 * @returns metadata about the derived table including the `.script` tag.
 */
export async function addAggregateToWorkspace(
  page: Page,
  options?: {via?: 'menu' | 'pivot-viewer'},
): Promise<OpenedTable> {
  return await page.evaluate(async (via) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const tv = grok.shell.tv;
    if (!tv?.dataFrame)
      throw new Error('addAggregateToWorkspace: no active TableView (open a base table first)');

    if (via === 'pivot-viewer') {
      tv.addViewer(DG.VIEWER.PIVOT_TABLE);
    } else {
      const cmd = DG.Func.find({name: 'CmdAggregateRows'})[0];
      if (!cmd) throw new Error('CmdAggregateRows not registered');
      await cmd.prepare({}).call(undefined, undefined, {processed: false});
    }
    await new Promise((r) => setTimeout(r, 1500));

    const btn = document.querySelector(
      '[name="button-ADD"].add-workspace-btn',
    ) as HTMLElement | null;
    if (!btn)
      throw new Error(
        'addAggregateToWorkspace: ADD button not visible — aggregate panel did not render',
      );
    btn.click();
    await new Promise((r) => setTimeout(r, 1500));

    const df = grok.shell.tv?.dataFrame;
    if (!df)
      throw new Error('addAggregateToWorkspace: ADD did not produce a TableView');
    const ti = df.getTableInfo?.();
    return {
      name: df.name,
      rowCount: df.rowCount,
      colCount: df.columns.length,
      tableInfoId: ti?.id,
      script: df.tags?.get?.('.script') ?? '',
    };
  }, options?.via ?? 'menu');
}

// ---------------------------------------------------------------------------
// 7. assertProvenanceScript — Gate E-PROV-01 inline check.
// ---------------------------------------------------------------------------

/**
 * Verify the active TableView's dataFrame carries a `.script` tag matching
 * the source-class pattern. Throws on mismatch.
 *
 * Use this immediately after a `openTableFrom*` call when you want a
 * spec-time guarantee that provenance was wired up — without it, save+
 * reopen lifecycles silently degrade to snapshot-only.
 *
 * Inline implementation of the proposed Gate E-PROV-01 rule. Critic-E may
 * eventually verify the same patterns statically; this helper provides
 * runtime verification today.
 *
 * @param page - Playwright Page.
 * @param sourceClass - Expected source class.
 * @param actualScript - Optional override of `df.tags['.script']` if the
 *   caller already captured it via `OpenedTable.script` (avoids a second
 *   roundtrip into the page).
 * @throws Error if the tag is missing or doesn't match the source-class regex.
 */
export async function assertProvenanceScript(
  page: Page,
  sourceClass: SourceClass,
  actualScript?: string,
): Promise<void> {
  const pattern = PROVENANCE_PATTERNS[sourceClass];
  const script = actualScript ?? await page.evaluate(() => {
    const grok = (window as any).grok;
    const df = grok.shell.tv?.dataFrame;
    return df?.tags?.get?.('.script') ?? '';
  });
  if (!script)
    throw new Error(
      `E-PROV-01: df.tags['.script'] is empty after Step 1. ` +
      `Expected ${sourceClass} provenance pattern ${pattern}. ` +
      `Likely cause: opener path did not engage the function-call recorder ` +
      `(see helpers/openers.ts file header for the canonical pattern).`,
    );
  if (!pattern.test(script))
    throw new Error(
      `E-PROV-01: df.tags['.script'] = "${script.slice(0, 200)}" ` +
      `does not match ${sourceClass} pattern ${pattern}.`,
    );
}
