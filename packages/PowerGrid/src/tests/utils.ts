import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {awaitCheck, expect} from '@datagrok-libraries/test/src/test';

/** Saves the layout of [view] to the server, applies what the server stored to a new view of [table]
 * that has a histogram in its place, and hands that view to [check]; the layout is deleted afterwards. */
export async function layoutRoundTrip(view: DG.TableView, table: DG.DataFrame,
  check: (view: DG.TableView) => Promise<void>): Promise<void> {
  const layout = view.saveLayout();
  await grok.dapi.layouts.save(layout);
  try {
    grok.shell.closeAll();
    const restored = grok.shell.addTableView(table);
    restored.addViewer(DG.VIEWER.HISTOGRAM);
    restored.loadLayout((await grok.dapi.layouts.find(layout.id))!);
    const histograms = Array.from(restored.viewers).filter((v) => v.type === DG.VIEWER.HISTOGRAM);
    expect(histograms.length, 0, 'the histogram outlived the layout');
    await check(restored);
  }
  finally {
    await grok.dapi.layouts.delete(layout);
  }
}

/** Saves [view] with its table as a project, closes everything, opens the project and hands its view to
 * [check]; the project, its table and its view are deleted afterwards. */
export async function projectRoundTrip(view: DG.TableView,
  check: (view: DG.TableView) => Promise<void>): Promise<void> {
  const project = DG.Project.create();
  project.name = `powergrid-tests-${Date.now()}`;
  const tableInfo = view.dataFrame.getTableInfo();
  const viewInfo = DG.ViewInfo.fromJson(view.saveLayout({saveWithData: true}).toJson());
  project.addChild(tableInfo);
  project.addChild(viewInfo);
  try {
    await grok.dapi.tables.uploadDataFrame(view.dataFrame);
    await grok.dapi.tables.save(tableInfo);
    await grok.dapi.views.save(viewInfo);
    await grok.dapi.projects.save(project);
    grok.shell.closeAll();
    await (await grok.dapi.projects.find(project.id)).open();
    await awaitCheck(() => (grok.shell.tv?.dataFrame?.rowCount ?? 0) > 0, 'the project opened no table', 15000);
    await check(grok.shell.tv);
  }
  finally {
    const saved: [DG.HttpDataSource<any>, string][] =
      [[grok.dapi.projects, project.id], [grok.dapi.views, viewInfo.id], [grok.dapi.tables, tableInfo.id]];
    for (const [source, id] of saved) {
      const entity = await source.find(id);
      if (entity)
        await source.delete(entity);
    }
  }
}
