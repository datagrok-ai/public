import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, awaitCheck, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
// import $ from 'cash-dom';

category('Layouts', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;

  before(async () => {
    df = grok.data.demo.demog(10);
    tv = grok.shell.addTableView(df);
  });

  test('TableView.saveLayout()', async () => {
    expect(tv.saveLayout() instanceof DG.ViewLayout, true);
  });

  test('ViewLayout.toJson()', async () => {
    const layout = tv.saveLayout();
    const json = layout.toJson();
    expect(typeof json, 'string');
    expect(JSON.parse(json).id, layout.id);
  });

  test('ViewLayout.fromJson(string)', async () => {
    const layout = tv.saveLayout();
    const json = layout.toJson();
    expect(DG.ViewLayout.fromJson(json) instanceof DG.ViewLayout, true);
  });

  test('ViewLayout.fromViewState(string)', async () => {
    const layout = tv.saveLayout();
    const state = layout.viewState;
    expect(typeof state, 'string');
    expect(DG.ViewLayout.fromViewState(state) instanceof DG.ViewLayout, true);
  });

  test('delete', async () => {
    grok.shell.windows.showToolbox = true;
    const tb = tv.toolbox;
    const layouts = Array.from(tb.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Layouts') as HTMLElement;
    if (!layouts.classList.contains('expanded'))
      layouts.click();
    tv.addViewer(DG.VIEWER.SCATTER_PLOT);
    await delay(500);
    let list: HTMLElement | null;
    await awaitCheck(() => {
      list = tb.querySelector('#layouts');
      return list !== null;
    }, '', 3000);
    const save = Array.from(tb.querySelectorAll('.ui-btn'))
      .find((el) => el.textContent === 'Save') as HTMLElement;
    await delay(3000);
    let num = list!.children.length + 1;
    save.click();
    await delay(6000);
    await awaitCheck(() => list!.children.length === num, 'Layout was not saved', 3000);
    num--;
    try {
      (list!.firstElementChild as HTMLElement).focus();
      list!.firstElementChild!.dispatchEvent(new MouseEvent('contextmenu'));
      await delay(100);
      await awaitCheck(() => document.querySelector('[d4-name="Upload"]') !== null, 'Cannot find context menu', 3000);
      (document.querySelector('[d4-name="Upload"]') as HTMLElement).click();
      await delay(1000);

      (list!.firstElementChild as HTMLElement).focus();
      list!.firstElementChild!.dispatchEvent(new MouseEvent('contextmenu'));
      await delay(100);
      await awaitCheck(() => document.querySelector('[d4-name="Delete"]') !== null, 'Cannot find context menu', 3000);
      (document.querySelector('[d4-name="Delete"]') as HTMLElement).click();
      let d: DG.Dialog;
      await awaitCheck(() => {
        const a = DG.Dialog.getOpenDialogs();
        if (a.length === 0) return false;
        d = a[0];
        return true;
      });
      const yes = Array.from(d!.root.querySelectorAll('.ui-btn'))
        .find((el) => el.textContent === 'DELETE') as HTMLElement;
      yes.click();
      await awaitCheck(() => list!.children.length === num, 'Layout was not deleted', 6000);
    } catch (e) {
      const l = await grok.dapi.layouts.getApplicable(df);
      if (l.length > num) {
        l.sort((a, b) => a.createdOn > b.createdOn ? -1 : 1);
        await grok.dapi.layouts.delete(l[0]);
      }
      throw e;
    }
  }, {timeout: 100000});

  after(async () => {
    grok.shell.closeAll();
  });
}, {clear: false, owner: 'aparamonov@datagrok.ai'});

category('Layouts: Apply', () => {
  const df: DG.DataFrame = grok.data.demo.demog(100);

  test('toolbox', async () => {
    const tv = grok.shell.addTableView(df);
    grok.shell.windows.showToolbox = true;
    const tb = tv.toolbox;
    const layouts = Array.from(tb.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Layouts') as HTMLElement;
    if (!layouts.classList.contains('expanded'))
      layouts.click();
    tv.addViewer(DG.VIEWER.SCATTER_PLOT);
    await delay(500);
    let list: HTMLElement | null;
    await awaitCheck(() => {
      list = tb.querySelector('#layouts');
      return list !== null;
    }, '', 3000);
    const save = Array.from(tb.querySelectorAll('.ui-btn'))
      .find((el) => el.textContent === 'Save') as HTMLElement;
    await delay(10000);
    const num = list!.children.length + 1;
    save.click();
    await delay(10000);
    await awaitCheck(() => list!.children.length === num, 'Layout was not saved', 3000);
    try {
      tv.resetLayout();
      await delay(100);
      await awaitCheck(() => [...tv.viewers].length === 1, 'Cannot reset layout', 3000);
      (list!.childNodes[1].childNodes[0] as HTMLElement).click();
      await awaitCheck(() => [...tv.viewers].length === 2,
        `Layout was not applied, expected 1 viewer, got ${[...tv.viewers].length - 1}`, 3000);
      await delay(100);
    } finally {
      const l = await grok.dapi.layouts.getApplicable(df);
      l.sort((a, b) => a.createdOn > b.createdOn ? -1 : 1);
      await grok.dapi.layouts.delete(l[0]);
    }
  }, {timeout: 100000});

  test('gallery', async () => {
    const tv = grok.shell.addTableView(df);
    grok.shell.windows.showToolbox = false;
    grok.shell.topMenu.find('View').find('Layout').find('Open Gallery').click();
    tv.addViewer(DG.VIEWER.SCATTER_PLOT);
    await delay(500);
    let list: HTMLElement | null;
    await awaitCheck(() => {
      list = document.querySelector('.panel-content #layouts');
      return list !== null;
    }, '', 3000);
    const num = list!.children. length + 1;
    grok.shell.topMenu.find('View').find('Layout').find('Save to Gallery').click();
    await awaitCheck(() => list!.children.length === num, 'Layout was not saved', 3000);
    try {
      tv.resetLayout();
      await delay(100);
      await awaitCheck(() => [...tv.viewers].length === 1, 'Cannot reset layout', 3000);
      (list!.firstElementChild?.firstElementChild as HTMLElement).click();
      await awaitCheck(() => [...tv.viewers].length === 2,
        `Layout was not applied, expected 1 viewer, got ${[...tv.viewers].length - 1}`, 3000);
      await delay(100);
    } finally {
      const l = await grok.dapi.layouts.getApplicable(df);
      l.sort((a, b) => a.createdOn > b.createdOn ? -1 : 1);
      await grok.dapi.layouts.delete(l[0]);
    }
  }, {skipReason: 'deadlock'});

  // test('drag-and-drop', async () => {
  //   const tv = grok.shell.addTableView(df);
  //   await delay(500);
  //   $(tv.root).trigger('drop', '');
  // });

  after(async () => {
    const p1 = Array.from(document.querySelectorAll('.panel-titlebar'))
      .find((el) => el.textContent === 'Viewer suggestions');
    if (p1)
      (p1.querySelector('.panel-titlebar-button-close') as HTMLElement).click();
    const p2 = Array.from(document.querySelectorAll('.panel-titlebar'))
      .find((el) => el.textContent === 'Layout suggestions');
    if (p2)
      (p2.querySelector('.panel-titlebar-button-close') as HTMLElement).click();
  });
}, { owner: 'aparamonov@datagrok.ai' });
