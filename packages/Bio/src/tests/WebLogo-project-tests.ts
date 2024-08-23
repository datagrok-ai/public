import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';

import {category, delay, expect, test, testViewer} from '@datagrok-libraries/utils/src/test';

import {awaitGrid} from './utils';
import {WebLogoViewer} from '../viewers/web-logo-viewer';

import {_package} from '../package-test';

const PROJECT_PREFIX: string = 'Tests.Bio.WebLogo-project';

category('WebLogo-project', () => {
  test('fasta', async () => {
    const prjName = `${PROJECT_PREFIX}.fasta`;
    const df = await _package.files.readCsv('tests/filter_FASTA.csv');
    const tableName = df.name;
    const col = df.getCol('fasta');
    await grok.data.detectSemanticTypes(df);
    const view = grok.shell.addTableView(df);
    const wlViewer = await df.plot.fromType('WebLogo',
      {sequenceColumnName: col.name}) as unknown as WebLogoViewer;
    view.dockManager.dock(wlViewer);
    await wlViewer.awaitRendered();
    await awaitGrid(view.grid);

    await uploadProject(prjName, df.getTableInfo(), view, df);
    grok.shell.closeAll();
    await delay(500);

    const prj2 = await grok.dapi.projects.open(prjName);
    const view2 = grok.shell.getTableView(tableName);

    const viewersA = wu(view2.viewers).toArray();
    expect(viewersA.length, 2);
    expect(viewersA.filter((f) => f.type === 'Grid').length, 1);
    const wlViewer2 = viewersA.find((f) => f.type === 'WebLogo') as WebLogoViewer;
    expect(!!wlViewer2, true);

    await awaitGrid(view.grid);
    await wlViewer2.awaitRendered();
    // TODO: Check WebLogo viewer content
  }, {skipReason: 'depends on 1.18'});
});

export async function uploadProject(projectName: string, tableInfo: DG.TableInfo,
  view: DG.TableView, df: DG.DataFrame): Promise<void> {
  const project = DG.Project.create();
  const viewLayout = view.saveLayout();

  project.name = projectName;
  project.addChild(tableInfo);
  project.addChild(viewLayout); // cause error

  await grok.dapi.layouts.save(view.saveLayout());
  await grok.dapi.tables.uploadDataFrame(df);
  await grok.dapi.tables.save(tableInfo);
  await grok.dapi.projects.save(project);
}
