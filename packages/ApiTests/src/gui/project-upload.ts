import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {uploadProject} from './gui-utils';

category('Projects', () => {

  test('project.upload', async () => {
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
