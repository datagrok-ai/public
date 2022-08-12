import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';

category('View events', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;
  const subs: Subscription[] = [];

  before(async () => {
    df = grok.data.demo.demog();
    for (const v of grok.shell.tableViews)
      if (v.table?.name === df.name)
        v.close();
    tv = grok.shell.addTableView(df);
  });
  test('On view adding', async () => {
    let v0: DG.View;
    let v1: DG.TableView;
    const openViewsCount = Array.from(grok.shell.views).length;
    try {
      subs.push(grok.events.onViewAdding.subscribe((view) => {
        expect(Array.from(grok.shell.views).length, openViewsCount);
        expect(grok.shell.view(view.name), null);
        expect(view.name, view.type === DG.TYPE.TABLE_VIEW ? v1.name : v0.name);
      }));
      v0 = grok.shell.newView('test');
      v1 = grok.shell.addTableView(df);
    } finally {
      closeViews(v0!, v1!);
    }
  });

  test('On view added', async () => {
    let v0: DG.View;
    let v1: DG.TableView;
    let v2: DG.ViewBase;
    const openViewsCount = Array.from(grok.shell.views).length;
    try {
      subs.push(grok.events.onViewAdded.subscribe((view) => {
        expect(Array.from(grok.shell.views).length > openViewsCount, true);
        expect(grok.shell.view(view.name) != null, true);
        expect(view.name, view.type === DG.TYPE.TABLE_VIEW ?
          v1.name : view.type === DG.View.FUNCTIONS ?
            v2.name : v0.name);
      }));
      v0 = grok.shell.newView('test');
      v1 = grok.shell.addTableView(df);
      v2 = grok.shell.addView(DG.View.createByType(DG.View.FUNCTIONS));
    } finally {
      closeViews(v0!, v1!, v2!);
    }
  });

  test('On view renamed', () => new Promise((resolve, reject) => {
    const name = tv.name;
    subs.push(grok.events.onViewRenamed.subscribe((view) => {
      if (view.name === 'new name')
        resolve('OK');
    }));
    try {
      tv.name = 'new name';
      setTimeout(() => reject('Failed to rename a view'), 10);
    } finally {
      tv.name = name;
    }
  }));

  test('On current view changed', async () => {
    const currentView = grok.shell.v;
    const newView = grok.shell.newView('test');
    try {
      subs.push(grok.events.onCurrentViewChanged.subscribe(() => {
        expect(currentView === grok.shell.v, false);
      }));
      grok.shell.v = newView;
    } finally {
      closeViews(newView);
    }
  });

  test('On view removed', async () => {
    const newView = grok.shell.newView('test');
    subs.push(grok.events.onViewRemoved.subscribe((view) => {
      expect(newView, view);
    }));
    newView.close();
  });

  test('On view layout generated', async () => {
    let baseLayout: DG.ViewLayout;
    subs.push(grok.events.onViewLayoutGenerated.subscribe((layout) => {
      expect(baseLayout, layout);
    }));
    baseLayout = tv.saveLayout();
  });

  test('On view layout applying', async () => {
    const v = grok.shell.addTableView(df);
    try {
      subs.push(grok.events.onViewLayoutApplying.subscribe((layout) => {
        expect(layout instanceof DG.ViewLayout, true);
        const state = JSON.parse(layout.viewState);
        const viewerElement = state.children.find((c: { [key: string]: any }) =>
          c.state.element && c.state.element.type === DG.VIEWER.SCATTER_PLOT);
        expect(viewerElement != null, true);
        expect(Array.from(v.viewers).length, 1);
      }));
      tv.scatterPlot();
      v.loadLayout(tv.saveLayout());
    } finally {
      tv.resetLayout();
      closeViews(v);
    }
  });

  test('On view layout applied', async () => {
    const v = grok.shell.addTableView(df);
    try {
      subs.push(grok.events.onViewLayoutApplied.subscribe((layout) => {
        expect(layout instanceof DG.ViewLayout, true);
        const state = JSON.parse(layout.viewState);
        const viewerElement = state.children.find((c: { [key: string]: any }) =>
          c.state.element && c.state.element.type === DG.VIEWER.HISTOGRAM);
        expect(viewerElement != null, true);
        expect(Array.from(v.viewers).length, 2);
      }));
      tv.histogram();
      v.loadLayout(tv.saveLayout());
    } finally {
      tv.resetLayout();
      closeViews(v);
    }
  });

  after(async () => {
    subs.forEach((sub) => sub.unsubscribe());
  });
});

function closeViews(...views: DG.ViewBase[]) {
  views.forEach((view) => {
    view.close();
    if (view instanceof DG.TableView && view.table)
      grok.shell.closeTable(view.table);
  });
}
