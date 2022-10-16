import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {waitForHTMLCollection} from './utils';


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
      .find(el => el.textContent === 'Filters') as HTMLElement;
    if (!filters.classList.contains('expanded')) await filters.click();

    const actions = Array.from(tb.querySelectorAll('div.d4-accordion-pane-header'))
      .find(el => el.textContent === 'Actions') as HTMLElement;
    if (!actions.classList.contains('expanded')) actions.click();
  });


  test('filters.all', async () => {
    const usapi = await grok.dapi.users
      .list()
      .then(users => users.length);

    const all = Array.from(tb.querySelectorAll('label'))
      .find(el => el.textContent === 'All');
    if (all === undefined) throw 'Error: cannot find All!';
    await all.click();

    const usui = (await waitForHTMLCollection('.grok-gallery-grid')).length;
    expect(usapi, usui);
  });


  test('filters.recentlyJoined', async () => {
    const usapi = await grok.dapi.users
      .filter('joined > -1w')
      .list()
      .then(users => users.length);

    const rj = Array.from(tb.querySelectorAll('label'))
      .find(el => el.textContent === 'Recently joined');
    if (rj === undefined) throw 'Error: cannot find Recently Joined!';
    await rj.click();

    const usui = (await waitForHTMLCollection('.grok-gallery-grid')).length;
    expect(usapi, usui);
  });


  test('actions.addUser', async () => {
    let user = DG.User.create();
    user.login = 'newlogin';
    user.status = DG.USER_STATUS.STATUS_NEW;
    user.firstName = 'new';
    user.lastName = 'user';

    // TODO: add save and delete when will work
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
    const user = (await waitForHTMLCollection('.grok-gallery-grid'))[0] as HTMLElement;
    grok.shell.windows.showProperties = true;
    user.click();

    const user_info = await new Promise<HTMLElement>((resolve, reject) => {
      const selector = '.grok-entity-prop-panel';
      if (document.querySelector(selector) !== null)
        if ((document.querySelector(selector) as HTMLElement).innerText.includes('Groups')) {
          return resolve(document.querySelector(selector) as HTMLElement);
      }
  
      const observer = new MutationObserver(() => {
        if (document.querySelector(selector) !== null)
          if ((document.querySelector(selector) as HTMLElement).innerText.includes('Groups')) {
            clearTimeout(timeout);
            observer.disconnect();
            resolve(document.querySelector(selector) as HTMLElement);
        }
      });
  
      const timeout = setTimeout(() => {
        observer.disconnect();
        reject(`Error: cannot find ${selector}!`)
      }, 3000
      );
  
      observer.observe(document.body, {
          childList: true,
          subtree: true
      });
    });
    
    const pict = user_info.querySelector('.grok-user-profile-picture');
    const desc = user_info.innerText;
    const b = (pict !== null) && desc.includes('Groups') && desc.includes('Joined')
    expect(b, true);
  });


  after(async () => {
    v.close();
  });
});