import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';

// What the ribbon's Save dialog keeps of a table, and how the project brings it back: with Data sync
// the table carries its creation script and is rebuilt by it on open (marked `.data-sync = success`),
// without it the uploaded snapshot is loaded. The tables are made as a menu makes them — the call
// runs unprocessed, so the shell adds its table and records the call in the creation script.

interface Saved {
  project: _DG.Project;
  tableInfo: _DG.TableInfo;
  viewInfo: _DG.ViewInfo;
}

interface Source {
  run: () => Promise<_DG.FuncCall>;
  column: string;
  rows: number;
}

const demog: Source = {
  run: () => recorded(DG.Func.find({package: 'core', name: 'OpenServerFile'})[0], {fullPath: 'System:DemoFiles/demog.csv'}),
  column: 'WEIGHT',
  rows: 5850,
};

const products: Source = {
  run: () => recorded(DG.Func.find({package: 'Dbtests', name: 'PostgresqlProducts'})[0], {}),
  column: 'unitprice',
  rows: 77,
};

function recorded(func: _DG.Func, params: {[name: string]: any}): Promise<_DG.FuncCall> {
  return func.prepare(params).call(false, undefined, {processed: false});
}

function runCore(name: string, params: {[name: string]: any}): Promise<_DG.FuncCall> {
  return recorded(DG.Func.find({package: 'core', name})[0], params);
}

async function openSource(source: Source): Promise<_DG.DataFrame> {
  grok.shell.closeAll();
  await source.run();
  await waitUntil(() => grok.shell.tables[0]?.rowCount === source.rows, 'the table the source opened', 30000);
  return grok.shell.tables[0];
}

async function saveProject(df: _DG.DataFrame, dataSync: boolean): Promise<Saved> {
  const project = DG.Project.create();
  project.name = `apitests-sync-${DG.Utils.randomString(8)}`;
  const tableInfo = df.getTableInfo();
  if (dataSync) {
    tableInfo.tags['.data-sync'] = 'sync';
    tableInfo.tags['.script'] = df.getTag('.script');
    tableInfo.tags['.VariableName'] = df.getTag('.VariableName');
  }
  const viewInfo = DG.ViewInfo.fromJson(grok.shell.tableView(df.name).saveLayout().toJson());
  project.addChild(tableInfo);
  project.addChild(viewInfo);
  await grok.dapi.tables.uploadDataFrame(df);
  await grok.dapi.tables.save(tableInfo);
  await grok.dapi.views.save(viewInfo);
  await grok.dapi.projects.save(project);
  return {project, tableInfo, viewInfo};
}

async function reopen(saved: Saved, rows: number): Promise<_DG.DataFrame> {
  grok.shell.closeAll();
  await (await grok.dapi.projects.find(saved.project.id)).open();
  await waitUntil(() => grok.shell.tables[0]?.rowCount === rows, 'the table the project opened', 30000);
  return grok.shell.tables[0];
}

async function remove(saved: Saved[]): Promise<void> {
  grok.shell.closeAll();
  for (const s of saved) {
    for (const [source, id] of [[grok.dapi.projects, s.project.id], [grok.dapi.views, s.viewInfo.id],
      [grok.dapi.tables, s.tableInfo.id]] as [_DG.HttpDataSource<any>, string][]) {
      const entity = await source.find(id).catch(() => null);
      if (entity)
        await source.delete(entity);
    }
  }
}

function expectColumnPlus(df: _DG.DataFrame, column: string, source: string, plus: number): void {
  const col = df.col(column)!;
  const src = df.col(source)!;
  for (let i = 0; i < df.rowCount; i++) {
    if (src.isNone(i))
      continue;
    expect(Math.abs(col.get(i) - (src.get(i) + plus)) < 1e-3, true,
      `row ${i}: ${column} = ${col.get(i)}, ${source} + ${plus} = ${src.get(i) + plus}`);
  }
}

async function columnChainSurvivesDataSync(source: Source): Promise<void> {
  const saved: Saved[] = [];
  try {
    const df = await openSource(source);
    const fileValue = df.get(source.column, 0);
    await runCore('AddNewColumn', {table: df, expression: `\${${source.column}} + 100`, name: 'Plus100'});
    await runCore('AddNewColumn', {table: df, expression: '${Plus100} + 100', name: 'Plus200'});
    await runCore('RenameColumn', {table: df, oldColName: source.column, newColName: 'Base'});
    expect(df.col('Plus100')!.getTag('formula'), '${Base} + 100');
    df.set('Base', 0, 500);
    await waitUntil(() => df.get('Plus200', 0) === 700, 'Plus200 of the edited row is not recalculated', 30000);

    saved.push(await saveProject(df, true));
    let reopened = await reopen(saved[0], source.rows);
    expect(reopened.getTag('.data-sync'), 'success');
    expect(reopened.col('Plus100')?.getTag('formula'), '${Base} + 100');
    expect(reopened.col('Plus200')?.getTag('formula'), '${Plus100} + 100');
    expect(reopened.get('Base', 0), fileValue, 'the edited cell is reread from the source');
    expectColumnPlus(reopened, 'Plus100', 'Base', 100);
    expectColumnPlus(reopened, 'Plus200', 'Plus100', 100);

    await runCore('RenameColumn', {table: reopened, oldColName: 'Base', newColName: 'Base2'});
    expect(reopened.col('Plus100')!.getTag('formula'), '${Base2} + 100');
    expect(reopened.col('Plus200')!.getTag('formula'), '${Plus100} + 100');
    reopened.set('Base2', 1, 400);
    await waitUntil(() => reopened.get('Plus200', 1) === 600,
      'Plus200 of the row edited after the reopen is not recalculated', 30000);
    expectColumnPlus(reopened, 'Plus100', 'Base2', 100);
    expectColumnPlus(reopened, 'Plus200', 'Plus100', 100);
  }
  finally {
    await remove(saved);
  }
}

category('Dapi: projects: data sync', () => {
  test('A table saved with data sync keeps its creation script and is rebuilt by it on open', async () => {
    const saved: Saved[] = [];
    try {
      const df = await openSource(demog);
      saved.push(await saveProject(df, true));
      const tags = (await grok.dapi.tables.find(saved[0].tableInfo.id)).tags;
      expect(tags['.data-sync'], 'sync');
      expect((tags['.script'] ?? '').includes('System:DemoFiles/demog.csv'), true, `the creation script: ${tags['.script']}`);
      const reopened = await reopen(saved[0], demog.rows);
      expect(reopened.getTag('.data-sync'), 'success');
    }
    finally {
      await remove(saved);
    }
  }, {timeout: 120000});

  test('A snapshot opens from the uploaded data, and with data sync turned on is rebuilt', async () => {
    const saved: Saved[] = [];
    try {
      const df = await openSource(demog);
      saved.push(await saveProject(df, false));
      expect((await grok.dapi.tables.find(saved[0].tableInfo.id)).tags['.data-sync'] ?? null, null);
      const snapshot = await reopen(saved[0], demog.rows);
      expect(snapshot.getTag('.data-sync') ?? null, null);
      expect((snapshot.getTag('.script') ?? '').includes('System:DemoFiles/demog.csv'), true,
        `the snapshot keeps its creation script: ${snapshot.getTag('.script')}`);

      const tableInfo = snapshot.getTableInfo();
      tableInfo.tags['.data-sync'] = 'sync';
      tableInfo.tags['.script'] = snapshot.getTag('.script');
      tableInfo.tags['.VariableName'] = snapshot.getTag('.VariableName');
      await grok.dapi.tables.save(tableInfo);
      expect((await grok.dapi.tables.find(saved[0].tableInfo.id)).tags['.data-sync'], 'sync');
      const rebuilt = await reopen(saved[0], demog.rows);
      expect(rebuilt.getTag('.data-sync'), 'success');
    }
    finally {
      await remove(saved);
    }
  }, {timeout: 120000});

  test('A calculated-column chain over a file survives a data-sync round trip and a source rename',
    () => columnChainSurvivesDataSync(demog), {timeout: 120000});

  test('A calculated-column chain over a query result survives a data-sync round trip and a source rename',
    () => columnChainSurvivesDataSync(products), {timeout: 120000});
}, {owner: 'aparamonov@datagrok.ai'});
