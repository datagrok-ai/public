import {after, before, category, expect, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


async function testTags(obj: any, apiPath: any, error: string) {
  obj.tag('apiteststag');
  await apiPath.save(obj);
  const len = await apiPath.filter('#apiteststag').list().then((els: any[])=> els.length);
  await apiPath.delete(obj);
  if (len !== 1) throw new Error(`Expected 1 ${error}, got ${len}`);
}

category('UI: Tags', () => {
  let v: DG.ViewBase;
  let v1: DG.View;

  before(async () => {
    const data: DG.TabPane = grok.shell.sidebar.getPane('Data');
    const prapi = data.content.querySelector('[data-view=projects]') as HTMLElement;
    await prapi.click();
    v = grok.shell.v;
  });

  test('filter.projects', async () => {
    const prapi = await grok.dapi.projects
      .filter('#demo')
      .list()
      .then((projects) => projects.length);
    await awaitCheck(() => {
      if (document.querySelector('.grok-items-view-counts') !== null)
        return /[0-9]+ \/ [0-9]+/g.test((document.querySelector('.grok-items-view-counts') as HTMLElement).innerText);
      return false;
    }, 'cannot load projects', 3000);
    const search = v.root.querySelector('.ui-input-editor') as HTMLInputElement;
    search.value = '#demo';
    search.dispatchEvent(new Event('input'));
    const regex = new RegExp(`[0-9]+ / ${prapi}`, 'g');
    await awaitCheck(() => {
      if (document.querySelector('.grok-items-view-counts') !== null)
        return regex.test((document.querySelector('.grok-items-view-counts') as HTMLElement).innerText);
      return false;
    }, 'number of projects does not match', 3000);
  });

  test('tag.editor', async () => {
    v1 = grok.shell.newView('Tag Editor');
    const editor = DG.TagEditor.create();
    editor.addTag('demo');
    editor.addTag('test');
    editor.addTag('1234');
    v1.append(editor.root);
    const tags = v1.root.querySelectorAll('.d4-tag').length;
    expect(tags, 3);
    v1.close();
  });

  test('tag.add', async () => {
    await testTags(DG.Project.create(), grok.dapi.projects, 'project');
    await testTags(DG.DataConnection.create('apitests', {
      dataSource: '',
      server: '',
      db: '',
    }), grok.dapi.connections, 'connection');
    await testTags(DG.Script.create('apitests'), grok.dapi.scripts, 'script');
  });

  after(async () => {
    v.close();
  });
});
