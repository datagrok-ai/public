import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
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
    let list: HTMLElement | null = null;
    // the list starts with its filter box and shows a loader while it fetches the layouts
    const loaded = () => (list = tb.querySelector('#layouts')) !== null && list.querySelector('.d4-update-shadow') === null;
    const cardIds = () => Array.from(list!.querySelectorAll('.grok-suggestions-chart-card'))
      .map((c) => c.getAttribute('data-layout-id'));
    await waitUntil(loaded, 'the layouts list did not load', 30000);
    const before = new Set(cardIds());
    const save = Array.from(tb.querySelectorAll('.ui-btn'))
      .find((el) => el.textContent === 'Save') as HTMLElement;
    save.click();
    let id: string | null | undefined;
    await waitUntil(() => loaded() && (id = cardIds().find((c) => !before.has(c))) != null, 'Layout was not saved', 30000);
    const card = () => list!.querySelector(`.grok-suggestions-chart-card[data-layout-id="${id}"]`) as HTMLElement | null;
    try {
      const openContextMenu = async () => {
        document.querySelectorAll('.d4-menu-popup').forEach((m) => m.remove());
        card()!.dispatchEvent(new MouseEvent('contextmenu'));
        await waitUntil(() => document.querySelector('.d4-menu-popup .d4-menu-item') !== null, 'Cannot find context menu');
      };

      await openContextMenu();
      // "Upload" only shows while the layout is dirty or not on the server; Delete is offered once it is uploaded
      const upload = document.querySelector('[d4-name="Upload"]') as HTMLElement | null;
      if (upload !== null) {
        // Upload fires the entity-saved event once the layout is on the server
        let uploaded = false;
        const sub = grok.events.onEvent('d4-entity-saved').subscribe(() => uploaded = true);
        try {
          upload.click();
          await waitUntil(() => uploaded, 'the layout was not uploaded');
        } finally {
          sub.unsubscribe();
        }
        await openContextMenu();
      }
      await waitUntil(() => document.querySelector('[d4-name="Delete"]') !== null, 'the layout menu has no Delete');
      (document.querySelector('[d4-name="Delete"]') as HTMLElement).click();
      let d: DG.Dialog;
      await waitUntil(() => {
        const a = DG.Dialog.getOpenDialogs();
        if (a.length === 0) return false;
        d = a[0];
        return true;
      }, 'the delete confirmation did not open');
      const yes = Array.from(d!.root.querySelectorAll('.ui-btn'))
        .find((el) => el.textContent === 'DELETE') as HTMLElement;
      yes.click();
      await waitUntil(() => card() === null, 'Layout was not deleted');
    } catch (e) {
      const saved = id ? await grok.dapi.layouts.find(id) : null;
      if (saved)
        await grok.dapi.layouts.delete(saved);
      throw e;
    }
  }, {timeout: 100000});
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
    let list: HTMLElement | null = null;
    const loaded = () => (list = tb.querySelector('#layouts')) !== null && list.querySelector('.d4-update-shadow') === null;
    const cardIds = () => Array.from(list!.querySelectorAll('.grok-suggestions-chart-card'))
      .map((c) => c.getAttribute('data-layout-id'));
    await waitUntil(loaded, 'the layouts list did not load', 30000);
    const before = new Set(cardIds());
    const save = Array.from(tb.querySelectorAll('.ui-btn'))
      .find((el) => el.textContent === 'Save') as HTMLElement;
    save.click();
    let id: string | null | undefined;
    await waitUntil(() => loaded() && (id = cardIds().find((c) => !before.has(c))) != null, 'Layout was not saved', 30000);
    try {
      tv.resetLayout();
      await waitUntil(() => [...tv.viewers].length === 1, 'Cannot reset layout');
      (list!.querySelector(`.grok-suggestions-chart-card[data-layout-id="${id}"] .grok-suggestions-chart-image`) as HTMLElement).click();
      await waitUntil(() => [...tv.viewers].length === 2,
        `Layout was not applied, expected 1 viewer, got ${[...tv.viewers].length - 1}`);
    } finally {
      const saved = id ? await grok.dapi.layouts.find(id) : null;
      if (saved)
        await grok.dapi.layouts.delete(saved);
    }
  }, {timeout: 100000});

  test('gallery', async () => {
    const tv = grok.shell.addTableView(df);
    grok.shell.windows.showToolbox = false;
    grok.shell.topMenu.find('View').find('Layout').find('Open Gallery').click();
    tv.addViewer(DG.VIEWER.SCATTER_PLOT);
    let list: HTMLElement | null = null;
    const loaded = () => (list = document.querySelector('.panel-content #layouts')) !== null &&
      list.querySelector('.d4-update-shadow') === null;
    const cardIds = () => Array.from(list!.querySelectorAll('.grok-suggestions-chart-card'))
      .map((c) => c.getAttribute('data-layout-id'));
    await waitUntil(loaded, 'the layouts gallery did not load', 30000);
    const before = new Set(cardIds());
    grok.shell.topMenu.find('View').find('Layout').find('Save to Gallery').click();
    let id: string | null | undefined;
    await waitUntil(() => loaded() && (id = cardIds().find((c) => !before.has(c))) != null,
      'Layout was not saved', 30000);
    try {
      tv.resetLayout();
      await waitUntil(() => [...tv.viewers].length === 1, 'Cannot reset layout');
      (list!.querySelector(`.grok-suggestions-chart-card[data-layout-id="${id}"] .grok-suggestions-chart-image`) as
        HTMLElement).click();
      await waitUntil(() => [...tv.viewers].length === 2,
        `Layout was not applied, expected 1 viewer, got ${[...tv.viewers].length - 1}`);
    } finally {
      const saved = id ? await grok.dapi.layouts.find(id) : null;
      if (saved)
        await grok.dapi.layouts.delete(saved);
    }
  }, {timeout: 100000});

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
