import {before, category, test, expect, awaitCheck, delay} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('UI: Sharing', () => {
  let v: DG.ViewBase;
  let testUser: DG.User;
  const entityName = 'apitestsentityshare';

  before(async () => {
    testUser = await grok.dapi.users.filter('login = "admin"').first();
  });

  test('scripts.ui', async () => {
    const newScript = DG.Script.create('');
    newScript.name = entityName;
    await grok.dapi.scripts.save(newScript);
    try {
      await testEntityUI(['Platform', 'Functions', 'Scripts'], '.d4-gallery-card');
    } finally {
      await grok.dapi.scripts.delete(newScript);
    }
  });

  test('projects.ui', async () => {
    const newProject = DG.Project.create();
    newProject.name = entityName;
    await grok.dapi.projects.save(newProject);
    try {
      await testEntityUI(['Dashboards'], '.d4-gallery-card');
    } finally {
      await grok.dapi.projects.delete(newProject);
    }
  });

  test('connections.ui', async () => {
    const newConnection = DG.DataConnection.create(entityName, {
      dataSource: 'Postgres',
      server: '',
      db: '',
    });
    await grok.dapi.connections.save(newConnection);
    try {
      await testEntityUI(['Databases'], '.d4-gallery-card');
    } finally {
      await grok.dapi.connections.delete(newConnection);
    }
  });

  test('scripts.api', async () => {
    await testEntityAPI(DG.Script.create(''), grok.dapi.scripts);
  });

  test('projects.api', async () => {
    await testEntityAPI(DG.Project.create(), grok.dapi.projects);
  });

  test('connections.api', async () => {
    await testEntityAPI(DG.DataConnection.create('', {
      dataSource: 'Files',
      server: '',
      db: ''}), grok.dapi.connections);
  });

  async function testEntityAPI(entity: DG.Entity, dapi: DG.HttpDataSource<any>) {
    await dapi.save(entity);
    try {
      // @ts-ignore
      expect((await grok.dapi.permissions.get(entity)).edit == undefined);
      await grok.dapi.permissions.grant(entity, testUser.group, true);
    } finally {
      // @ts-ignore
      expect((await grok.dapi.permissions.get(entity)).edit.length, 1);
      await dapi.delete(entity);
    }
  }

  async function testEntityUI(path: string[], selector: string, friendlyName?: string) {
    let treeGroupToClick: DG.TreeViewGroup | DG.TreeViewNode = grok.shell.browsePanel.mainTree;
    for (let i = 0; i < path.length; i++) {
      const name = path[i];
      treeGroupToClick = i === path.length - 1 ? (treeGroupToClick as DG.TreeViewGroup).children.find((g) => g.text === name)! :
        (treeGroupToClick as DG.TreeViewGroup).getOrCreateGroup(name);
    }
    treeGroupToClick.root.click();
    await awaitCheck(() => (document.querySelector('.grok-gallery-grid')?.children?.length ?? 0) > 0,
      'cannot load gallery grid', 3000);
    v = grok.shell.v;
    const gallery = v.root;
    const search = gallery.querySelector('.ui-input-editor') as HTMLInputElement;
    search.value = entityName;
    search.dispatchEvent(new MouseEvent('input'));
    await awaitCheck(() => document.querySelector('.grok-gallery-grid')?.children?.length === 1,
      'more than one testing entity present', 3000);
    await delay(2000);
    let entity = gallery.querySelector('.grok-gallery-grid')!.children[0];
    if (entity.className !== selector.slice(1)) entity = entity.querySelector(selector) ?? entity;
    entity.dispatchEvent(new MouseEvent('contextmenu'));
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'cannot find context menu', 5000);
    const menu = document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup');
    const share = Array.from(menu!.querySelectorAll('.d4-menu-item.d4-menu-item-vert'))
      .find((el) => el.textContent === 'Share...') as HTMLElement;
    share.click();
    await awaitCheck(() => DG.Dialog.getOpenDialogs().find((d) => d.root.style.position === 'fixed' &&
      (d.title.toLowerCase() === `share ${entityName}` || d.title.toLowerCase() === `share ${(friendlyName ?? '').toLowerCase()}`)) !== undefined, 'cannot find dialog', 10000);
    const shareDialog = DG.Dialog.getOpenDialogs().find((d) => d.root.style.position === 'fixed' &&
      (d.title.toLowerCase() === `share ${entityName}` || d.title.toLowerCase() === `share ${(friendlyName ?? '').toLowerCase()}`));
    const shareDialogRoot = shareDialog!.root;
    await awaitCheck(() => shareDialogRoot.querySelector('.d4-user-selector-input') !== null,
      'cannot enter user for sharing', 3000);
    const inputUser = shareDialogRoot.querySelector('.d4-user-selector-input') as HTMLInputElement;
    inputUser.value = 'Admin';
    inputUser.dispatchEvent(new MouseEvent('input'));
    await awaitCheck(() => (shareDialogRoot.querySelector('.d4-user-selector-drop-down') as HTMLElement)
      ?.style?.visibility === '', 'cannot find user for sharing', 3000);
    await delay(100); // no way to handle it
    inputUser.dispatchEvent(new KeyboardEvent('keydown', {keyCode: 40}));
    await delay(100); // no way to handle it
    inputUser.dispatchEvent(new KeyboardEvent('keydown', {keyCode: 40}));
    await delay(100); // no way to handle it
    inputUser.dispatchEvent(new KeyboardEvent('keydown', {keyCode: 13}));
    await delay(300); // no way to handle it
    shareDialogRoot.querySelector('textarea')!.value = 'Message1!';
    (shareDialogRoot.querySelector('.ui-btn.ui-btn-ok.enabled') as HTMLElement).click();
    DG.Balloon.closeAll();
    await delay(100);
    await awaitCheck(() => (document.querySelector('.d4-balloon.info') as HTMLElement)?.innerText === 'Shared',
      'cannot find info balloon', 3000);
    v.close();
  }
}, { owner: 'aparamonov@datagrok.ai' });
