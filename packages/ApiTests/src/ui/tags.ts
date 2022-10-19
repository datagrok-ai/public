import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {waitForHTMLElement} from './utils';


category('UI: Tags', () => {
  let v: DG.ViewBase;
  let v1: DG.View;

  before(async () => {
    const data: DG.TabPane = grok.shell.sidebar.getPane('Data');
    const prapi = data.content.querySelector('[data-view=projects]') as HTMLElement;
    await prapi.click();
    v = grok.shell.v;
  });


  test('projectsFilter', async () => {
		const prapi = await grok.dapi.projects
			.filter('#demo')
			.list()
			.then(projects => projects.length);

		await waitForHTMLElement('.grok-items-view-counts', /[0-9]+ \/ [0-9]+/g, 'Error: cannot load Projects!')
		const search = v.root.querySelector('.ui-input-editor') as HTMLInputElement;
		search.value = '#demo';
		search.dispatchEvent(new Event('input'));
		const regex = new RegExp(`[0-9]+ / ${prapi}`, 'g');

		await waitForHTMLElement('.grok-items-view-counts', regex, 'Number of projects does not match!');
  });


  test('tagEditor', async () => {
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


  after(async () => {
    v.close();
  });
});