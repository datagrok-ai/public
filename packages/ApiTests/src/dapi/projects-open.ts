import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('Dapi: projects: open', () => {
  test('A project opens by its friendly name, short name and namespace-qualified name', async () => {
    grok.shell.closeAll();
    const df = grok.data.demo.demog(10);
    const view = grok.shell.addTableView(df);
    const tableInfo = df.getTableInfo();
    const viewInfo = DG.ViewInfo.fromJson(view.saveLayout().toJson());
    const project = DG.Project.create();
    project.name = `P${DG.Utils.randomString(12)}`;
    project.friendlyName = `apitests open ${DG.Utils.randomString(8)}`;
    project.addChild(tableInfo);
    project.addChild(viewInfo);
    try {
      await grok.dapi.tables.uploadDataFrame(df);
      await grok.dapi.tables.save(tableInfo);
      await grok.dapi.views.save(viewInfo);
      await grok.dapi.projects.save(project);
      const saved = await grok.dapi.projects.find(project.id);
      expect(saved.nqName.includes(':'), true, `the project has no namespace: ${saved.nqName}`);
      const dotted = saved.nqName.replace(/:/g, '.');
      const dottedFriendly = dotted.slice(0, dotted.length - saved.name.length) + saved.friendlyName;
      for (const name of [saved.name, saved.friendlyName, saved.nqName, dotted, dottedFriendly]) {
        grok.shell.closeAll();
        await grok.dapi.projects.open(name);
        expect(grok.shell.projects.some((p) => p.id === saved.id), true, `open('${name}') did not open the project`);
      }
      grok.shell.closeAll();
      const foreign = `ApiTestsNobody${DG.Utils.randomString(8)}.${saved.name}`;
      expect(await grok.dapi.projects.open(foreign).then(() => null, (e) => e), `Project ${foreign} not found`);
    }
    finally {
      grok.shell.closeAll();
      for (const [source, id] of [[grok.dapi.projects, project.id], [grok.dapi.views, viewInfo.id],
        [grok.dapi.tables, tableInfo.id]] as [_DG.HttpDataSource<any>, string][]) {
        const entity = await source.find(id).catch(() => null);
        if (entity)
          await source.delete(entity);
      }
    }
  }, {timeout: 60000});
}, {owner: 'aparamonov@datagrok.ai'});
