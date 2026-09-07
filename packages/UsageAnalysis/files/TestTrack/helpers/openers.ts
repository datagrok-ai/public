import {Page} from '@playwright/test';

export interface OpenedTable {

  name: string;

  rowCount: number;

  colCount: number;

  tableInfoId: string;

  script: string;
}

export type SourceClass = 'files' | 'db_query' | 'db_table' | 'script' | 'derived';

export const PROVENANCE_PATTERNS: Record<SourceClass, RegExp> = {
  files: /^\w+\s*=\s*OpenFile\("[^"]+"\)/,
  db_query: /^\w+\s*=\s*[\w:]+\(\)/,
  db_table: /^\w+\s*=\s*DbQuery\([\w:]+,\s*"\w+"/,
  script: /^\w+\s*=\s*[\w:]+\([\w\s,'"=.-]*\)/,

  derived: /^\w+\s*=\s*Aggregate\("[\w\s]+",/,
};

export async function resetShell(page: Page): Promise<void> {
  await page.evaluate(() => {
    const grok = (window as any).grok;
    try { grok?.shell?.closeAll?.(); } catch (_) {}
    document.querySelectorAll('.d4-dialog, .d4-toast, .d4-balloon, .d4-menu-popup')
      .forEach((e) => { try { (e as any).remove?.(); } catch (_) {} });
  });
  await page.waitForTimeout(300);
}

export async function openTableFromFile(
  page: Page,
  fullPath: string,
): Promise<OpenedTable> {
  return await page.evaluate(async (p) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const fns = DG.Func.find({name: 'OpenFile'});
    if (!fns?.length) throw new Error('OpenFile function not registered');

    const fullPathForOpen = p.replace(/^([^:/]+):/, '$1.');
    await fns[0].prepare({fullPath: fullPathForOpen}).call(undefined, undefined, {processed: false});

    let df: any = null;
    for (let i = 0; i < 24; i++) {
      const tv = grok.shell.tv;
      if (tv?.dataFrame && typeof tv.addViewer === 'function') {
        df = tv.dataFrame;
        break;
      }
      await new Promise((r) => setTimeout(r, 500));
    }
    if (!df) throw new Error(`OpenFile("${p}") did not produce a TableView (12s settle)`);
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

export interface ProvisionedScript {

  scriptId: string;

  resolvedName: string;

  resolvedNqName: string;
}

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

export async function deleteProvisionedScript(
  page: Page,
  scriptId: string,
): Promise<void> {
  await page.evaluate(async (id) => {
    try {
      const grok = (window as any).grok;
      const s = await grok.dapi.scripts.find(id);
      if (s) await grok.dapi.scripts.delete(s);
    } catch (_) {  }
  }, scriptId).catch(() => {});
}

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

export const SYSTEM_DATAGROK_NQNAME = 'System:Datagrok';

export const SYSTEM_DATAGROK_QUERIES = {

  GROUPS_SAMPLE: 'select id, friendly_name as name, personal from groups limit 20',

  GROUPS_RELATIONS: 'select parent_id, child_id from groups_relations limit 30',
};

export const SYSTEM_DATAGROK_DB_TABLE = {
  schemaName: 'public',
  tableName: 'groups',
};

export interface SystemDatagrokConnection {

  id: string;

  nqName: string;
}

export async function getSystemDatagrokConnection(
  page: Page,
): Promise<SystemDatagrokConnection> {
  return await page.evaluate(async () => {
    const grok = (window as any).grok;
    const conn = await grok.functions.eval('System:Datagrok');
    if (!conn?.id)
      throw new Error('System:Datagrok connection not found — broken deploy');
    return {id: conn.id, nqName: conn.nqName ?? 'System:Datagrok'};
  });
}

export interface ProvisionedQuery {

  queryId: string;

  queryNqName: string;

  resolvedName: string;

  cleanup: () => Promise<void>;
}

export async function provisionSystemDatagrokQuery(
  page: Page,
  options: {nameStem: string; sql: string},
): Promise<ProvisionedQuery> {

  const stamp = Date.now();
  const camelStem = options.nameStem.replace(/[_\s]+([a-z0-9])/gi, (_, c) => c.toUpperCase());
  const uniqueName = `${camelStem}${stamp}`;

  const provisioned = await page.evaluate(async ({n, s}) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const conn = await grok.functions.eval('System:Datagrok');
    if (!conn) throw new Error('System:Datagrok connection not found');
    const q = conn.query(n, s);
    const saved = await grok.dapi.queries.save(q);

    const cap = (x: string) => (x ? x.charAt(0).toUpperCase() + x.slice(1) : x);
    const candidates = [saved?.name, saved?.shortName, n, cap(n)].filter(Boolean);
    let fn: any = null;
    for (let i = 0; i < 20 && !fn; i++) {
      for (const cand of candidates) {
        fn = DG.Func.find({name: cand})?.[0];
        if (fn) break;
      }
      if (!fn) await new Promise((r) => setTimeout(r, 300));
    }

    return {
      queryId: saved?.id,
      resolvedName: fn?.name ?? saved?.name ?? n,
      queryNqName: fn?.nqName ?? saved?.nqName
        ?? (saved?.namespace ? `${saved.namespace}:${saved.name}` : n),
    };
  }, {n: uniqueName, s: options.sql});

  return {
    ...provisioned,
    cleanup: () => deleteProvisionedQuery(page, provisioned.queryId),
  };
}

export async function deleteProvisionedQuery(
  page: Page,
  queryId: string,
): Promise<void> {
  await page.evaluate(async (id) => {
    try {
      const grok = (window as any).grok;
      const q = await grok.dapi.queries.find(id);
      if (q) await grok.dapi.queries.delete(q);
    } catch (_) {  }
  }, queryId).catch(() => {});
}
