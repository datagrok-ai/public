import {after, before, category, expect, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('UI: Users', () => {
  let v: DG.ViewBase;
  let tb: HTMLElement;

  before(async () => {
    const mng: DG.TabPane = grok.shell.sidebar.getPane('Manage');
    const usel = mng.content.querySelector('[data-view=users]') as HTMLElement;
    await usel.click();
    v = grok.shell.v;
    tb = v.toolbox;
    const filters = Array.from(tb.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Filters') as HTMLElement;
    if (!filters.classList.contains('expanded')) await filters.click();
    const actions = Array.from(tb.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Actions') as HTMLElement;
    if (!actions.classList.contains('expanded')) actions.click();
  });

  test('filters.all', async () => {
    const usapi = await grok.dapi.users
      .list()
      .then((users) => users.length);
    const regex = new RegExp(`[0-9]+ / ${usapi}`, 'g');
    const all = Array.from(tb.querySelectorAll('label'))
      .find((el) => el.textContent === 'All');
    if (all === undefined) throw new Error('cannot find All');
    await all.click();
    await awaitCheck(() => {
      if (document.querySelector('.grok-items-view-counts') !== null)
        return regex.test((document.querySelector('.grok-items-view-counts') as HTMLElement).innerText);
      return false;
    }, 'number of users does not match', 3000);
  });

  test('filters.recentlyJoined', async () => {
    const usapi = await grok.dapi.users
      .filter('joined > -1w')
      .list()
      .then((users) => users.length);
    const regex = new RegExp(`[0-9]+ / ${usapi}`, 'g');
    const rj = Array.from(tb.querySelectorAll('label'))
      .find((el) => el.textContent === 'Recently joined');
    if (rj === undefined) throw new Error('cannot find Recently Joined');
    await rj.click();
    setTimeout(async () => {
      await awaitCheck(() => {
        if (document.querySelector('.grok-items-view-counts') !== null)
          return regex.test((document.querySelector('.grok-items-view-counts') as HTMLElement).innerText);
        return false;
      }, 'number of users does not match', 3000);
    }, 10);
  });

  test('actions.addUser', async () => {
    const user = DG.User.create();
    user.login = 'newlogin';
    user.status = DG.USER_STATUS.STATUS_NEW;
    user.firstName = 'new';
    user.lastName = 'user';
    // TODO: add save and delete when will work
  }, {skipReason: 'blocked by GROK-11318'});

  test('actions.addServiceUser', async () => {
    const asu = Array.from(tb.querySelectorAll('label'))
      .find((el) => el.textContent === 'Add Service User...');
    if (asu === undefined) throw new Error('cannot find Add Service User');
    await asu.click();
    const diag = document.getElementsByClassName('d4-dialog');
    if (diag.length === 0) throw new Error('Add Service User does not work');
    const cancel = diag[0].querySelector('[class="ui-btn ui-btn-ok"]') as HTMLElement;
    cancel.click();
  });

  test('actions.inviteFriend', async () => {
    const iaf = Array.from(tb.querySelectorAll('label'))
      .find((el) => el.textContent === 'Invite a Friend...');
    if (iaf === undefined) throw new Error('cannot find Invite a Friend');
    await iaf.click();
    const diag = document.getElementsByClassName('d4-dialog');
    if (diag.length === 0) throw new Error('Invite a Friend does not work');
    const cancel = diag[0].querySelector('[class="ui-btn ui-btn-ok"]') as HTMLElement;
    cancel.click();
  });

  test('user.panel', async () => {
    const all = Array.from(tb.querySelectorAll('label'))
      .find((el) => el.textContent === 'All');
    if (all === undefined) throw new Error('cannot find All');
    await all.click();
    await awaitCheck(() => {
      if (document.querySelector('.grok-gallery-grid') !== null)
        return document.querySelector('.grok-gallery-grid')!.children.length !== 0;
      return false;
    }, 'cannot load all users', 3000);
    const user = document.querySelector('.grok-gallery-grid')!.children[0] as HTMLElement;
    grok.shell.windows.showProperties = true;
    user.click();
    // await awaitCheck(() => {
    //   if (document.querySelector('.grok-entity-prop-panel') !== null)
    //     return document.querySelector('.grok-entity-prop-panel')!.querySelectorAll('*').length === 48;
    //   return false;
    // }, 'cannot find user\'s widgets', 3000);
    // const expand = document.querySelector('.grok-icon.fal.fa-chevron-square-down') as HTMLElement;
    // expand.click();
    const regex = new RegExp('Groups', 'g');
    await awaitCheck(() => {
      if (document.querySelector('.grok-entity-prop-panel') !== null)
        return regex.test((document.querySelector('.grok-entity-prop-panel') as HTMLElement).innerText);
      return false;
    }, 'cannot find user\'s groups', 3000);
    const userInfo = document.querySelector('.grok-entity-prop-panel') as HTMLElement;
    const pict = userInfo.querySelector('.grok-user-profile-picture');
    const desc = userInfo.innerText;
    const b = (pict !== null) && desc.includes('Groups') && desc.includes('Joined');
    expect(b, true);
  });

  after(async () => {
    v.close();
  });
});
