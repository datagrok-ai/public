import {after, before, category, expect, test, awaitCheck, delay} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('UI: Groups', () => {
  before(async () => {
    const mng: DG.TabPane = grok.shell.sidebar.getPane('Browse');
    mng.header.click();
    let platform: any;
    await awaitCheck(() => {
      platform = Array.from(document.querySelectorAll('div.d4-tree-view-node'))
        .find((el) => el.textContent === 'Platform');
      return platform !== undefined;
    }, '', 2000);
    if ((platform.nextElementSibling as HTMLElement).style.display === 'none')
      (platform.firstElementChild as HTMLElement).click();
    await delay(100);
    const groups = Array.from(document.querySelectorAll('div.d4-tree-view-item-label'))
      .find((el) => el.textContent === 'Groups') as HTMLElement;
    groups.click(); 
  });
 
  /*
  test('filters.all', async () => {
    const grapi = await grok.dapi.groups
      .list()
      .then((groups) => groups.filter((group) => !group.personal && !group.dart.y).length);
    const regex = new RegExp(`[0-9]+ / ${grapi - 1}`, 'g');
    const all = Array.from(tb.querySelectorAll('label'))
      .find((el) => el.textContent === 'All');
    if (all === undefined) throw new Error('cannot find All!');
    await all.click();
    await awaitCheck(() => {
      if (document.querySelector('.grok-items-view-counts') !== null)
        return regex.test((document.querySelector('.grok-items-view-counts') as HTMLElement).innerText);
      return false;
    }, 'number of groups does not match', 3000);
  });

  test('filters.mine', async () => {
    const grapi = await grok.dapi.groups
      .filter('children.child = @current')
      .list()
      .then((groups) => groups.filter((group) => !group.personal && !group.dart.y).length);
    const regex = new RegExp(`[0-9]+ / ${grapi - 1}`, 'g');
    const mine = Array.from(tb.querySelectorAll('label'))
      .find((el) => el.textContent === 'Mine');
    if (mine === undefined) throw new Error('cannot find Mine!');
    await mine.click();
    await awaitCheck(() => {
      if (document.querySelector('.grok-items-view-counts') !== null)
        return regex.test((document.querySelector('.grok-items-view-counts') as HTMLElement).innerText);
      return false;
    }, 'number of groups does not match', 3000);
  });
  */

  test('actions.createNewGroup', async () => {
    const cng = Array.from(document.querySelectorAll('.ui-btn'))
      .find((el) => el.textContent === 'New Group...');
    if (cng === undefined) throw new Error('cannot find New Group button');
    (cng as HTMLElement).click();
    await awaitCheck(() => DG.Dialog.getOpenDialogs().length === 1, 'Create New Group dialog was not shown', 1000);
    const diag = DG.Dialog.getOpenDialogs()[0];
    const cancel = diag.root.querySelector('[class="ui-btn ui-btn-ok"]') as HTMLElement;
    cancel.click();
    const mainGroup = await grok.dapi.groups.createNew('APITests Main Group');
    const subGroup = DG.Group.create('APITests Sub Group');
    let error;
    try {
      subGroup.includeTo(mainGroup);
      const testUser = await grok.dapi.users.filter('login = "System"').first();
      subGroup.addMember(testUser.group);
      await grok.dapi.groups.saveRelations(subGroup);
      expect(subGroup.members.length, 1);
      expect(subGroup.adminMembers.length, 1);
      expect(subGroup.memberships.length, 2);
      expect(subGroup.adminMemberships.length, 0);
    } catch (e: any) {
      error = e;
    } finally {
      await grok.dapi.groups.delete(subGroup);
      await grok.dapi.groups.delete(mainGroup);
    }
    if (error) throw error;
  });

  test('group.panel', async () => {
    await awaitCheck(() => {
      if (document.querySelector('.grok-gallery-grid') !== null)
        return document.querySelector('.grok-gallery-grid')!.children.length !== 0;
      return false;
    }, 'cannot load all users', 10000);
    const group = document.querySelector('.grok-gallery-grid')!.children[0] as HTMLElement;
    grok.shell.windows.showProperties = true;
    group.click();
    const regex = new RegExp('MANAGE', 'g');
    await awaitCheck(() => {
      if (document.querySelector('.grok-entity-prop-panel') !== null)
        return regex.test((document.querySelector('.grok-entity-prop-panel') as HTMLElement).innerText);
      return false;
    }, 'error in .grok-entity-prop-panel', 10000);
    const groupInfo = document.querySelector('.grok-entity-prop-panel') as HTMLElement;
    const memb = groupInfo.innerText.includes('Members');
    expect(memb, true);
    const manage = document.querySelector('.d4-pane-manage-button') as HTMLElement;
    manage.click();
    await awaitCheck(() => document.querySelector('.d4-dialog') !== null, 'cannot find manage dialog', 3000);
    const diag = document.querySelector('.d4-dialog') as HTMLElement;
    const cancel = diag.querySelectorAll('[class="ui-btn ui-btn-ok"]')[1] as HTMLElement;
    cancel.click();
  });

  after(async () => {
    grok.shell.closeAll();
  });
}, {clear: false});
