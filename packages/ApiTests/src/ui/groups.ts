import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {waitForHTMLCollection, waitForHTMLElement} from './utils';


category('UI: Groups', () => {
  let v: DG.ViewBase;
  let tb: HTMLElement;

  before(async () => {
    const mng: DG.TabPane = grok.shell.sidebar.getPane('Manage');
    const usel = mng.content.querySelector('[data-view=groups]') as HTMLElement;
    await usel.click();
    v = grok.shell.v;
    tb = v.toolbox;

    const filters = Array.from(tb.querySelectorAll('div.d4-accordion-pane-header'))
      .find(el => el.textContent === 'Filters') as HTMLElement;
    if (!filters.classList.contains('expanded')) await filters.click();

    const actions = Array.from(tb.querySelectorAll('div.d4-accordion-pane-header'))
      .find(el => el.textContent === 'Actions') as HTMLElement;
    if (!actions.classList.contains('expanded')) actions.click();
  });


  test('filters.all', async () => {
    const grapi = await grok.dapi.groups
      .list()
      .then(groups => groups.filter(group => !group.personal && !group.dart.y).length);
    const regex = new RegExp(`[0-9]+ / ${grapi - 1}`, 'g');

    const all = Array.from(tb.querySelectorAll('label'))
      .find(el => el.textContent === 'All');
    if (all === undefined) throw 'Error: cannot find All!';
    await all.click();

    await waitForHTMLElement('.grok-items-view-counts', regex, 'Number of groups does not match!');
  });


  test('filters.mine', async () => {
    const grapi = await grok.dapi.groups
      .filter('children.child = @current')
      .list()
      .then(groups => groups.filter(group => !group.personal && !group.dart.y).length);
    const regex = new RegExp(`[0-9]+ / ${grapi - 1}`, 'g');

    const mine = Array.from(tb.querySelectorAll('label'))
      .find(el => el.textContent === 'Mine');
    if (mine === undefined) throw 'Error: cannot find Mine!';
    await mine.click();

    await waitForHTMLElement('.grok-items-view-counts', regex, 'Number of groups does not match!');
  });


  test('actions.createNewGroup', async () => {
    const cng = Array.from(tb.querySelectorAll('label'))
      .find(el => el.textContent === 'Create New Group...');
    if (cng === undefined) throw 'Error: cannot find Create New Group!';
    await cng.click();

    const diag = document.getElementsByClassName('d4-dialog');
    if (diag.length === 0) throw 'Error: Create New Group does not work!';
    const cancel = diag[0].querySelector('[class="ui-btn ui-btn-ok"]') as HTMLElement;
    cancel.click();
  });


  test('group.panel', async () => {
    const group = (await waitForHTMLCollection('.grok-gallery-grid'))[0] as HTMLElement;
    grok.shell.windows.showProperties = true;
    group.click();
    
    const regex = new RegExp('MANAGE', 'g');
    const groupInfo = await waitForHTMLElement('.grok-entity-prop-panel', regex, 'Error in .grok-entity-prop-panel!'); 
    const memb = groupInfo.innerText.includes('Members');
    expect(memb, true);
    const manage = document.querySelector('.d4-pane-manage-button') as HTMLElement;
    manage.click();
    
    const diag = await waitForHTMLElement('.d4-dialog', /./g, 'Error in .d4-dialog!');
    const cancel = diag.querySelectorAll('[class="ui-btn ui-btn-ok"]')[1] as HTMLElement;
    cancel.click();
  });


  after(async () => {
    v.close();
  });
});