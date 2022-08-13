import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test, delay} from '@datagrok-libraries/utils/src/test';

category('Viewers', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;
  let coreViewerTypes: string[];

  before(async () => {
    coreViewerTypes = Object.values(DG.VIEWER).filter((v) =>
      v != DG.VIEWER.TIMELINES &&
      v != DG.VIEWER.GLOBE &&
      v != DG.VIEWER.SCATTER_PLOT_3D &&
      v != DG.VIEWER.GOOGLE_MAP &&
      v != DG.VIEWER.SHAPE_MAP);
    df = grok.data.demo.demog();
    tv = grok.shell.addTableView(df);
  });

  test('addViewer(ViewerType)', async () => {
    try {
      for (const viewerType of coreViewerTypes) {
        console.log(`Adding ${viewerType}`);
        const viewer = await addViewerAndWait(tv, viewerType);//tv.addViewer(viewerType);
        if (!(viewer instanceof DG.Viewer))
          throw `TableView.addViewer('${viewerType}') should add a Viewer instance`;
        await delay(100);

        viewer.removeFromView();
        console.log(`Removed ${viewerType}`);
      }
    } finally {
      closeViewers(tv);
    }
  });

  test('addViewer(Viewer)', async () => {
    try {
      for (const viewerType of coreViewerTypes) {
        console.log(`Adding ${viewerType}`);
        const viewer = await addViewerAndWait(tv, DG.Viewer.fromType(viewerType, df));
        await delay(100);
        viewer.removeFromView();
        console.log(`Removed ${viewerType}`);
      }
    } finally {
      closeViewers(tv);
    }
  });

  test('close', async () => {
    tv.scatterPlot();
    tv.barChart();
    expect(Array.from(tv.viewers).length, 3);
    closeViewers(tv);
    expect(Array.from(tv.viewers).length, 1);
  });

  test('getViewerTypes', async () => {
    const registeredViewers = DG.Viewer.getViewerTypes();
    expect(coreViewerTypes.every((t) => registeredViewers.includes(t)), true);
  });

  test('fromType', async () => {
    for (const viewerType of coreViewerTypes) {
      const viewer = DG.Viewer.fromType(viewerType, df);
      if (!(viewer instanceof DG.Viewer))
        throw `Viewer.fromType('${viewerType}', df) should add a Viewer instance`;
      expect(viewer.table.id, df.id);
    }
  });

  test('ScatterPlotViewer.zoom', async () => {
    const sp = DG.Viewer.scatterPlot(df);
    tv.addViewer(sp);
    try {
      sp.zoom(10, 10, 100, 100);
    } finally {
      sp.close();
    }
  });

  test('ScatterPlotViewer.onZoomed', async () => {
    const sp = DG.Viewer.scatterPlot(df);
    tv.addViewer(sp);
    try {
      let rectangle: any;
      sp.onZoomed.subscribe((r) => rectangle = r);
      sp.zoom(10, 10, 100, 100);
      expect(rectangle instanceof DG.Rect, true);
    } finally {
      sp.close();
    }
  });

  after(async () => {
    tv.close();
    grok.shell.closeTable(df);
  });
});

function closeViewers(view: DG.TableView) {
  Array.from(view.viewers).slice(1).forEach((v) => v.close());
}

function addViewerAndWait(tv: DG.TableView, viewerType: string | DG.Viewer): Promise<DG.Viewer> {
  return new Promise((resolve, reject) => {
    const sub = grok.events.onViewerAdded.subscribe((data) => {
      // @ts-ignore
      if ((data.args.viewer as DG.Viewer).type == viewerType || data.args.viewer == viewerType) {
        sub.unsubscribe();
        // @ts-ignore
        resolve(data.args.viewer);
      }
    });
    tv.addViewer(viewerType);
    setTimeout(() => {
      sub.unsubscribe();
      // eslint-disable-next-line prefer-promise-reject-errors
      reject('timeout');
    }, 100);
  });
}
