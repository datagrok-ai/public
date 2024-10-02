import { after, before, category, expect, test, awaitCheck, delay } from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


async function testTags(obj: any, apiPath: any, error: string) {
  obj.tag('apiteststag');
  const lenBefore = await apiPath.filter('#apiteststag').list().then((els: any[]) => els.length);
  await apiPath.save(obj);
  const len = await apiPath.filter('#apiteststag').list().then((els: any[]) => els.length);
  await apiPath.delete(obj);
  if (lenBefore + 1 !== len) throw new Error(`Expected ${lenBefore} ${error}, got ${len}`);
}

category('UI: Tags', () => {
  let v1: DG.View;

  before(async () => {
    const mng: DG.TabPane = grok.shell.sidebar.getPane('Browse');
    mng.header.click();
    let groups: any;
    await delay(1000);
    await awaitCheck(() => {
      groups = Array.from(document.querySelectorAll('div.d4-tree-view-item-label'))
        .find((el) => el.textContent === 'Dashboards');
      return groups !== undefined;
    }, '', 2000);
    groups.click();
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
    await delay(3000);
    const search = grok.shell.v.root.querySelector('.grok-gallery-search-bar .ui-input-root .ui-input-editor') as HTMLInputElement;
    search.value = '#demo';
    search.dispatchEvent(new Event('input'));
    console.log('regex')
    console.log(`[0-9]+ / ${prapi}`)
    const regex = new RegExp(`[0-9]+ \\/ ((${prapi})|(...))`, 'g');
    await awaitCheck(() => {
      if (document.querySelector('.grok-items-view-counts') !== null) {
        console.log(document.querySelector('.grok-items-view-counts'));
        return regex.test((document.querySelector('.grok-items-view-counts') as HTMLElement).innerText);
      }
      return false;
    }, 'number of projects does not match', 7000);
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
      dataSource: 'Files',
      server: '',
      db: '',
    }), grok.dapi.connections, 'connection');
    await testTags(DG.Script.create('apitests'), grok.dapi.scripts, 'script');
  }, { timeout: 100000 });

  after(async () => {
    grok.shell.closeAll();
  });
}, { clear: false });
