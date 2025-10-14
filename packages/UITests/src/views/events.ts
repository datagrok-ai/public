import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test, testEvent} from '@datagrok-libraries/utils/src/test';

category('View: Events', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;

  before(async () => {
    df = grok.data.demo.demog();
    for (const v of grok.shell.tableViews) {
      if (v.table?.name === df.name)
        v.close();
    }
    tv = grok.shell.addTableView(df);
  });

  test('onViewAdding', async () => {
    const openViewsCount = Array.from(grok.shell.views).length;
    // @ts-ignore
    await testEvent<DG.View>(grok.events.onViewAdding, (view) => {
      expect(Array.from(grok.shell.views).length, openViewsCount);
      expect(view.name, 'test');
    }, () => {
      grok.shell.newView('test').close();
    });
  });

  test('onViewAdded 1', async () => {
    const openViewsCount = Array.from(grok.shell.views).length;
    // @ts-ignore
    await testEvent<DG.View>(grok.events.onViewAdded, (view) => {
      expect(Array.from(grok.shell.views).length > openViewsCount, true);
      expect(grok.shell.view(view.name) != null, true);
      expect(view.name, 'test');
    }, () => {
      grok.shell.newView('test').close();
    });
  });

  test('onViewAdded 2', async () => {
    // @ts-ignore
    await testEvent<DG.View>(grok.events.onViewAdded, (view) => {
      expect(view.name, 'Functions');
    }, () => {
      grok.shell.addView(DG.View.createByType(DG.View.FUNCTIONS));
    });
  });

  test('onViewAdded 3', async () => {
    // @ts-ignore
    await testEvent<DG.View>(grok.events.onViewAdded, (view) => {
      expect(view instanceof DG.TableView, true);
      expect((view as DG.TableView).dataFrame.rowCount, 10000);
    }, () => {
      grok.shell.addTableView(grok.data.demo.demog());
    });
  });

  test('onViewRenamed', async () => {
    const name = tv.name;
    try {
      // @ts-ignore
      await testEvent<DG.View>(grok.events.onViewRenamed, (view) => {
        expect(view.name, 'new name');
      }, () => {
        tv.name = 'new name';
      });
    } finally {
      tv.name = name;
    }
  });

  test('onCurrentViewChanged', async () => {
    const currentView = grok.shell.v;
    // @ts-ignore
    await testEvent<DG.View>(grok.events.onCurrentViewChanged, (view) => {
      expect(currentView === grok.shell.v, false);
      expect(grok.shell.v.name, 'test');
    }, () => {
      grok.shell.v = grok.shell.newView('test');
    }, 100);
  });

  test('onViewRemoved', async () => {
    const newView1 = grok.shell.newView('test1');
    const newView2 = grok.shell.newView('test2');
    try {
      // @ts-ignore
      await testEvent<DG.View>(grok.events.onViewRemoved, (view) => {
        expect(view.name, 'test2');
      }, () => {
        newView2.close();
      });
    } finally {
      newView1.close();
      newView2.close();
    }
  });

  test('onViewLayoutGenerated', async () => {
    // @ts-ignore
    await testEvent<DG.ViewInfo>(grok.events.onViewLayoutGenerated, (layout) => {
      expect(layout != null, true);
    }, () => {
      tv.saveLayout();
    });
  });

  test('onViewLayoutApplying', async () => {
    const v = grok.shell.addTableView(df);
    try {
      // @ts-ignore
      await testEvent<DG.ViewInfo>(grok.events.onViewLayoutApplying, (layout) => {
        expect(layout instanceof DG.ViewInfo, true);
        const state = JSON.parse(layout.viewState);
        const viewerElement = state.children.find((c: { [key: string]: any }) =>
          c.state.element && c.state.element.type === DG.VIEWER.SCATTER_PLOT);
        expect(viewerElement != null, true);
        expect(Array.from(v.viewers).length, 1);
      }, () => {
        tv.scatterPlot();
        v.loadLayout(tv.saveLayout());
      });
    } finally {
      tv.resetLayout();
      v.close();
    }
  });

  test('onViewLayoutApplied', async () => {
    const v = grok.shell.addTableView(df);
    try {
      // @ts-ignore
      await testEvent<DG.ViewInfo>(grok.events.onViewLayoutApplied, (layout) => {
        expect(layout instanceof DG.ViewInfo, true);
        const state = JSON.parse(layout.viewState);
        const viewerElement = state.children.find((c: { [key: string]: any }) =>
          c.state.element && c.state.element.type === DG.VIEWER.HISTOGRAM);
        expect(viewerElement != null, true);
        expect(Array.from(v.viewers).length, 2);
      }, () => {
        grok.shell.v = tv;
        tv.histogram();
        v.loadLayout(tv.saveLayout());
      });
    } finally {
      tv.resetLayout();
    }
  });

  after(async () => {
    grok.shell.closeAll();
  });
}, {clear: false, owner: 'aparamonov@datagrok.ai' });
