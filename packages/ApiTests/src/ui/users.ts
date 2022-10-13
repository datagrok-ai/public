import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
// import {checkHTMLElement} from './utils';


category('UI: Users', () => {
  let v: DG.ViewBase;
  let tb: HTMLElement;

  before(async () => {
    const mng: DG.TabPane = grok.shell.sidebar.getPane('Manage');
    const usel = mng.content.querySelector('[data-view=users]') as HTMLElement;
    await usel.click();
    v = grok.shell.v;
    tb = v.toolbox;

    // const filters = Array.from(tb.querySelectorAll('div.d4-accordion-pane-header'))
    //   .find(el => el.textContent === 'Filters') as HTMLElement;
    // if (!filters.classList.contains('expanded')) await filters.click();

    const actions = Array.from(tb.querySelectorAll('div.d4-accordion-pane-header'))
      .find(el => el.textContent === 'Actions') as HTMLElement;
    if (!actions.classList.contains('expanded')) actions.click();
  });


  test('filters.all', async () => {
    await grok.dapi.users
      .list()
      .then(users => expect(users.length > 0, true));
  });


  test('filters.recentlyJoined', async () => {
    await grok.dapi.users
      .filter('joined > -1w')
      .list()
      .then(users => expect(users.length >= 0, true));
  });


  test('actions.addUser', async () => {
    let user = DG.User.create();
    user.login = 'newlogin';
    user.status = DG.USER_STATUS.STATUS_NEW;
    user.firstName = 'new';
    user.lastName = 'user';
  });


  test('actions.addServiceUser', async () => {
    const asu = Array.from(tb.querySelectorAll('label'))
      .find(el => el.textContent === 'Add Service User...');
    if (asu === undefined) throw 'Error: cannot find Add Service User!';
    await asu.click();

    const diag = document.getElementsByClassName('d4-dialog');
    if (diag.length === 0) throw 'Error: Add Service User does not work!';
    const cancel = diag[0].querySelector('[class="ui-btn ui-btn-ok"]') as HTMLElement;
    cancel.click();
  });


  test('actions.inviteFriend', async () => {
    const iaf = Array.from(tb.querySelectorAll('label'))
      .find(el => el.textContent === 'Invite a Friend...');
    if (iaf === undefined) throw 'Error: cannot find Invite a Friend!';
    await iaf.click();

    const diag = document.getElementsByClassName('d4-dialog');
    if (diag.length === 0) throw 'Error: Invite a Friend does not work!';
    const cancel = diag[0].querySelector('[class="ui-btn ui-btn-ok"]') as HTMLElement;
    cancel.click();
  });


  test('user.panel', async () => {
    const user = await grok.dapi.users
      .list()
      .then(users => users.find(obj => obj.login === 'test')) as DG.User;
    expect(user.picture === null, false);
    expect(user.group === null, false);
  });


  after(async () => {
    v.close();
  });
});