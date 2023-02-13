import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {uploadProject} from './gui-utils';

category('Projects', () => {

  test('Project: Create', async () => {
    const projectName = 'Test Project';
    const project = DG.Project.create();
    project.name = projectName;
    expect(project.name, projectName);
    expect(project instanceof DG.Project, true);
  });

  test('Project: Populate', async () => {
    const project = DG.Project.create();
    const ti1 = grok.data.demo.demog(1000).getTableInfo();
    const ti2 = grok.data.demo.randomWalk(1000).getTableInfo();
    const connection = await grok.dapi.connections.first();
    expect(project.children.length, 0);
    project.addChild(ti1);
    project.addChild(ti2);
    project.addLink(connection);
    expect(project.children.length, 2);
    expect(project.children[0], ti1);
    expect(project.links.length, 1);
    expect(project.links[0], connection);
    project.removeChild(ti1);
    project.removeLink(connection);
    expect(project.children.length, 1);
    expect(project.children[0], ti2);
    expect(project.links.length, 0);
  });

  test('Project: Upload', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    
    await uploadProject('Test upload project', demog.getTableInfo(), v, demog);
    
    grok.shell.closeAll(); await delay(500);
    
    await grok.dapi.projects.open('Test upload project');
    expect(grok.shell.v.name, 'demog 1000')
    
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test upload project').first())
    expect(await grok.dapi.projects.filter('Test upload project').first(), undefined);
  });
});
