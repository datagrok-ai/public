import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, awaitCheck, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';


category('Layouts', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;

  before(async () => {
    df = grok.data.demo.demog();
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

  after(async () => {
    grok.shell.closeAll();
  });
}, {clear: false});

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
    tv.addViewer(DG.VIEWER.HISTOGRAM);
    tv.addViewer(DG.VIEWER.LINE_CHART);
    await delay(2000);
    let list: HTMLElement | null;
    await awaitCheck(() => {
      list = tb.querySelector('#layouts');
      return list !== null;
    }, '', 3000);
    const save = Array.from(tb.querySelectorAll('.ui-btn'))
      .find((el) => el.textContent === 'Save') as HTMLElement;
    const num = list!.children.length + 1;
    save.click();
    await awaitCheck(() => tb.querySelector('#layouts')?.children.length === num, 'Layout was not saved', 3000);
    try {
      tv.resetLayout();
      await delay(100);
      await awaitCheck(() => [...tv.viewers].length === 1, 'Cannot reset layout', 3000);
      (list!.firstElementChild?.firstElementChild as HTMLElement).click();
      await awaitCheck(() => [...tv.viewers].length === 4,
        `Layout was not applied, expected 3 viewers, got ${[...tv.viewers].length - 1}`, 3000);
      await delay(100);
    } finally {
      const l = await grok.dapi.layouts.getApplicable(df);
      l.sort((a, b) => a.createdOn > b.createdOn ? -1 : 1);
      await grok.dapi.layouts.delete(l[0]);
    }
  });

  test('gallery', async () => {
  });
});
