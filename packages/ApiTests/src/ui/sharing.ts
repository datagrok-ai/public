import {before, category, test, expect, awaitCheck} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('UI: Sharing', () => {
  let v: DG.ViewBase;
  let testUser: DG.User;
  const entityName = 'apitestsentityshare';

  before(async () => {
    testUser = await grok.dapi.users.filter('login = "test.user"').first();
  });

  test('scripts.ui', async () => {
    const newScript = DG.Script.create('');
    newScript.name = entityName;
    await grok.dapi.scripts.save(newScript);
    try {
      await testEntityUI('Functions', 'scripts', '.d4-gallery-card');
    } finally {
      await grok.dapi.scripts.delete(newScript);
    }
  });

  test('projects.ui', async () => {
    const newProject = DG.Project.create();
    newProject.name = entityName;
    await grok.dapi.projects.save(newProject);
    try {
      await testEntityUI('Data', 'projects', '.d4-gallery-card');
    } finally {
      await grok.dapi.projects.delete(newProject);
    }
  });

  test('connections.ui', async () => {
    const newConnection = DG.DataConnection.create(entityName, {
      dataSource: '',
      server: '',
      db: '',
    });
    await grok.dapi.connections.save(newConnection);
    try {
      await testEntityUI('Manage', 'connections', '.d4-link-label');
    } finally {
      await grok.dapi.connections.delete(newConnection);
    }
  });

  test('scripts.api', async () => {
    testEntityAPI(DG.Script.create(''), grok.dapi.scripts);
  });

  test('projects.api', async () => {
    testEntityAPI(DG.Project.create(), grok.dapi.projects);
  });

  test('connections.api', async () => {
    testEntityAPI(DG.DataConnection.create('', {
      dataSource: '',
      server: '',
      db: ''}), grok.dapi.connections);
  });

  async function testEntityAPI(entity: DG.Entity, dapi: DG.HttpDataSource<any>) {
    await dapi.save(entity);
    try {
      // @ts-ignore
      expect((await grok.dapi.permissions.get(entity)).edit, undefined);
      await grok.dapi.permissions.grant(entity, testUser.group, true);
    } finally {
      // @ts-ignore
      expect((await grok.dapi.permissions.get(entity)).edit.length, 1);
      await dapi.delete(entity);
    }
  }

  async function testEntityUI(paneName: string, elName: string, selector: string) {
    const pane = grok.shell.sidebar.getPane(paneName);
    const el = pane.content.querySelector(`[data-view=${elName}]`) as HTMLElement;
    el.click();
    await awaitCheck(() => document.querySelector('.grok-gallery-grid')!.children.length > 0,
      'cannot load gallery grid', 3000);
    v = grok.shell.v;
    const gallery = v.root;
    const search = gallery.querySelector('.ui-input-editor') as HTMLInputElement;
    search.value = entityName;
    search.dispatchEvent(new MouseEvent('input'));
    await awaitCheck(() => document.querySelector('.grok-gallery-grid')!.children.length === 1,
      'more than one testing entity present', 3000);
    let entity = gallery.querySelector('.grok-gallery-grid')!.children[0];
    if (entity.className !== selector.slice(1)) entity = entity.querySelector(selector)!;
    entity.dispatchEvent(new MouseEvent('contextmenu'));
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'cannot find context menu', 3000);
    const menu = document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup');
    const share = Array.from(menu!.querySelectorAll('.d4-menu-item.d4-menu-item-vert'))
      .find((el) => el.textContent === 'Share...') as HTMLElement;
    share.click();
    await awaitCheck(() => DG.Dialog.getOpenDialogs().find((d) => d.root.style.position === 'fixed' &&
      d.title.toLowerCase() === `share ${entityName}`) !== undefined, 'cannot find dialog', 10000);
    const shareDialog = DG.Dialog.getOpenDialogs().find((d) => d.root.style.position === 'fixed' &&
      d.title.toLowerCase() === `share ${entityName}`);
    const shareDialogRoot = shareDialog!.root;
    await awaitCheck(() => shareDialogRoot.querySelector('.user-selector-input') !== null,
      'cannot enter user for sharing', 3000);
    const inputUser = shareDialogRoot.querySelector('.user-selector-input') as HTMLInputElement;
    inputUser.value = 'test.user';
    inputUser.dispatchEvent(new MouseEvent('input'));
    await awaitCheck(() => (shareDialogRoot.querySelector('.user-selector-drop-down') as HTMLElement)
      .style.visibility === '', 'cannot find user for sharing', 3000);
    inputUser.dispatchEvent(new KeyboardEvent('keydown', {keyCode: 13}));
    shareDialogRoot.querySelector('textarea')!.value = 'Message1!';
    (shareDialogRoot.querySelector('.ui-btn.ui-btn-ok.enabled') as HTMLElement).click();
    // To Do: add balloons cleaner
    await awaitCheck(() => (document.querySelector('.d4-balloon.info') as HTMLElement).innerText === 'Shared',
      'cannot find info balloon', 3000);
    v.close();
  }
});
