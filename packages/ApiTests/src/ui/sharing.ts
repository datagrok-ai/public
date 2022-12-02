import {category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('UI: Sharing', () => {
  let v: DG.ViewBase;
  const entityName = 'apitestsentityshare';
  let error;

  test('scripts.ui', async () => {
    error = null;
    const newScript = DG.Script.create('');
    newScript.name = entityName;
    await grok.dapi.scripts.save(newScript);
    try {
      await testEntity('Functions', 'scripts', '.d4-gallery-card');
    } catch (e: any) {
      error = e;
    }
    await grok.dapi.scripts.delete(newScript);
    if (error) throw error;
  });

  test('projects.ui', async () => {
    error = null;
    const newProject = DG.Project.create();
    newProject.name = entityName;
    await grok.dapi.projects.save(newProject);
    try {
      await testEntity('Data', 'projects', '.d4-gallery-card');
    } catch (e: any) {
      error = e;
    }
    await grok.dapi.projects.delete(newProject);
    if (error) throw error;
  });

  test('connections.ui', async () => {
    error = null;
    const newConnection = DG.DataConnection.create(entityName, {
      dataSource: '',
      server: '',
      db: '',
    });
    await grok.dapi.connections.save(newConnection);
    try {
      await testEntity('Manage', 'connections', '.d4-link-label');
    } catch (e: any) {
      error = e;
    }
    await grok.dapi.connections.delete(newConnection);
    if (error) throw error;
  });

  async function testEntity(paneName: string, elName: string, selector: string) {
    const pane = grok.shell.sidebar.getPane(paneName);
    const el = pane.content.querySelector(`[data-view=${elName}]`) as HTMLElement;
    el.click();
    await awaitCheck(() => document.querySelector('.grok-gallery-grid')!.children.length > 0); // 3000
    v = grok.shell.v;
    const gallery = v.root;
    const search = gallery.querySelector('.ui-input-editor') as HTMLInputElement;
    search.value = entityName;
    search.dispatchEvent(new MouseEvent('input'));
    await awaitCheck(() => document.querySelector('.grok-gallery-grid')!.children.length === 1); // 3000
    let entity = gallery.querySelector('.grok-gallery-grid')!.children[0];
    if (entity.className !== selector.slice(1)) entity = entity.querySelector(selector)!;
    entity.dispatchEvent(new MouseEvent('contextmenu'));
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null);//3000
    const menu = document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup');
    const share = Array.from(menu!.querySelectorAll('.d4-menu-item.d4-menu-item-vert'))
      .find((el) => el.textContent === 'Share...') as HTMLElement;
    share.click();
    await awaitCheck(() => DG.Dialog.getOpenDialogs().find((d) => d.root.style.position === 'fixed' &&
      d.title.toLowerCase() === `share ${entityName}`) !== undefined); // 10000
    const shareDialog = DG.Dialog.getOpenDialogs().find((d) => d.root.style.position === 'fixed' &&
      d.title.toLowerCase() === `share ${entityName}`);
    const shareDialogRoot = shareDialog!.root;
    await awaitCheck(() => shareDialogRoot.querySelector('.user-selector-input') !== null); // 3000
    const inputUser = shareDialogRoot.querySelector('.user-selector-input') as HTMLInputElement;
    inputUser.value = 'test.user';
    inputUser.dispatchEvent(new MouseEvent('input'));
    await awaitCheck(() => (shareDialogRoot.querySelector('.user-selector-drop-down') as HTMLElement)
      .style.visibility === ''); // 3000
    inputUser.dispatchEvent(new KeyboardEvent('keydown', {keyCode: 13}));
    shareDialogRoot.querySelector('textarea')!.value = 'Message1!';
    (shareDialogRoot.querySelector('.ui-btn.ui-btn-ok.enabled') as HTMLElement).click();
    // To Do: add balloons cleaner
    await awaitCheck(() => (document.querySelector('.d4-balloon.info') as HTMLElement).innerText === 'Shared'); // 3000
    v.close();
  }
});
