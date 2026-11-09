import { after, before, category, expect, test, waitUntil } from '@datagrok-libraries/test/src/test';
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
    grok.shell.windows.showBrowse = true;
    let groups: any;
    await waitUntil(() => {
      groups = Array.from(document.querySelectorAll('div.d4-tree-view-item-label'))
        .find((el) => el.textContent === 'Dashboards');
      return groups !== undefined;
    }, 'the Dashboards tree item did not appear');
    groups.click();
  });

  test('filter.projects', async () => {
    const prapi = await grok.dapi.projects
      .filter('#demo')
      .list()
      .then((projects) => projects.length);
    const counts = () => (grok.shell.v.root.querySelector('.grok-items-view-counts') as HTMLElement | null)?.innerText ?? '';
    await waitUntil(() => /^([0-9]+ of )?[0-9]+$/.test(counts()), 'cannot load projects', 30000);
    const search = grok.shell.v.root.querySelector('.grok-gallery-search-bar .ui-input-root .ui-input-editor') as HTMLInputElement;
    search.value = '#demo';
    search.dispatchEvent(new Event('input'));
    // one page shows just the total, several "<loaded> of <total>"
    await waitUntil(() => new RegExp(`^([0-9]+ of )?${prapi}$`).test(counts()),
      `the #demo search should list ${prapi} projects`, 30000);
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
    // Projects now require a non-empty name on save (datlas projects_service guard).
    const project = DG.Project.create();
    project.name = 'apitests';
    await testTags(project, grok.dapi.projects, 'project');
    await testTags(DG.DataConnection.create('apitests', {
      dataSource: 'Files',
      server: '',
      db: '',
    }), grok.dapi.connections, 'connection');
    await testTags(DG.Script.create('apitests'), grok.dapi.scripts, 'script');
  }, { timeout: 100000 });
}, { clear: false, owner: 'drizhinashvili@datagrok.ai' });
