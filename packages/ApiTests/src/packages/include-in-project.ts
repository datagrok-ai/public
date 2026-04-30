import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {after, category, expect, test} from '@datagrok-libraries/test/src/test';

category('Projects: include-in-project', () => {
  const createdProjectIds: string[] = [];
  const createdTableIds: string[] = [];
  const createdViewIds: string[] = [];

  after(async () => {
    for (const id of createdViewIds) {
      const vi = await grok.dapi.views.find(id);
      if (vi != null)
        await grok.dapi.views.delete(vi);
    }
    for (const id of createdProjectIds) {
      const proj = await grok.dapi.projects.find(id);
      if (proj != null)
        await grok.dapi.projects.delete(proj);
    }
    for (const id of createdTableIds) {
      const ti = await grok.dapi.tables.find(id);
      if (ti != null)
        await grok.dapi.tables.delete(ti);
    }
  });

  test('TableInfo.includeInProject default true; round-trips false locally', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromInt32Array('a', new Int32Array([1, 2, 3]))]);
    df.name = 't1';
    const ti = DG.TableInfo.fromDataFrame(df);
    expect(ti.includeInProject, true, 'absent flag must default to true');
    ti.includeInProject = false;
    expect(ti.includeInProject, false, 'after setter false, getter returns false');
    expect(ti.tags['.include-in-project'], 'false', 'underlying metaParam must store the literal string "false"');
    ti.includeInProject = true;
    expect(ti.includeInProject, true, 'after setter true, getter returns true');
    expect(ti.tags['.include-in-project'] == null, true, 'underlying metaParam must be removed when toggled back to true');
  }, {owner: 'dkovalyov@datagrok.ai'});

  test('TableInfo IncludeInProject=false persists across project save/find', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromInt32Array('a', new Int32Array([1, 2, 3]))]);
    df.name = `t_intermediate_${Date.now()}`;
    const tableId = await grok.dapi.tables.uploadDataFrame(df);
    createdTableIds.push(tableId);
    const ti = await grok.dapi.tables.find(tableId);
    ti.includeInProject = false;
    await grok.dapi.tables.save(ti);

    const project = DG.Project.create();
    project.name = `IncludeInProjectRoundTrip_${Date.now()}`;
    project.addChild(ti);
    const saved = await grok.dapi.projects.save(project);
    createdProjectIds.push(saved.id);

    const reloaded = await grok.dapi.projects.find(saved.id);
    const reloadedTi = reloaded.children.find((e) => e.id === tableId) as _DG.TableInfo | undefined;
    expect(reloadedTi != null, true, `reloaded child with id ${tableId} must be present`);
    expect(reloadedTi!.includeInProject, false, 'TableInfo IncludeInProject=false must round-trip through dapi.projects.save → find');
  }, {owner: 'dkovalyov@datagrok.ai', timeout: 30000});

  test('ViewInfo IncludeInProject=false persists across project save/find', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromInt32Array('a', new Int32Array([1, 2, 3]))]);
    df.name = `t_for_view_${Date.now()}`;
    const tableId = await grok.dapi.tables.uploadDataFrame(df);
    createdTableIds.push(tableId);
    const ti = await grok.dapi.tables.find(tableId);

    const tv = grok.shell.addTableView(df);
    const vi = tv.getInfo() as unknown as _DG.ViewInfo;
    vi.includeInProject = false;
    await grok.dapi.views.save(vi);
    createdViewIds.push(vi.id);

    const project = DG.Project.create();
    project.name = `ViewInfoIncludeRoundTrip_${Date.now()}`;
    project.addChild(ti);
    project.addChild(vi);
    const saved = await grok.dapi.projects.save(project);
    createdProjectIds.push(saved.id);

    const reloaded = await grok.dapi.projects.find(saved.id);
    const reloadedVi = reloaded.children.find((e) => e instanceof DG.ViewInfo) as _DG.ViewInfo | undefined;
    expect(reloadedVi != null, true, 'reloaded project must contain the ViewInfo');
    expect(reloadedVi!.includeInProject, false, 'ViewInfo IncludeInProject=false (stored in userData) must round-trip through dapi.views.save + dapi.projects.save → find');
  }, {owner: 'dkovalyov@datagrok.ai', timeout: 30000});
}, {owner: 'dkovalyov@datagrok.ai'});
