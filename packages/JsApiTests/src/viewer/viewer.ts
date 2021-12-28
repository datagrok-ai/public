import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { after, before, category, expect, test } from '../test';

category('Viewer', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;
  const coreViewerTypes = Object.values(DG.VIEWER).filter((v) => v !== DG.VIEWER.TIMELINES && v !== DG.VIEWER.GLOBE);
  const coreViewers: DG.Viewer[] = [];

  before(async () => {
    df = grok.data.demo.demog();
    tv = grok.shell.addTableView(df);
  });

  test('TableView.addViewer(ViewerType)', async () => {
    for (let viewerType of coreViewerTypes) {
      const viewer = tv.addViewer(viewerType);
      if (!(viewer instanceof DG.Viewer))
        throw `TableView.addViewer('${viewerType}') should add a Viewer instance`;
    }
    const attachedViewers = Array.from(tv.viewers);
    if (attachedViewers.length - 1 !== coreViewerTypes.length)
      throw 'TableView.addViewer failed to attach some viewers to the table view';
  });

  test('Viewer.close()', async () => {
    Array.from(tv.viewers).slice(1).forEach(v => v.close());
    expect(Array.from(tv.viewers).length, 1);
  });

  test('Viewer.getViewerTypes()', async () => {
    const registeredViewers = DG.Viewer.getViewerTypes();
    expect(coreViewerTypes.every((t) => registeredViewers.includes(t)), true);
  });

  test('Viewer.fromType(ViewerType, DataFrame)', async () => {
    for (let viewerType of coreViewerTypes) {
      const viewer = DG.Viewer.fromType(viewerType, df);
      if (!(viewer instanceof DG.Viewer))
        throw `Viewer.fromType('${viewerType}', df) should add a Viewer instance`;
      expect(viewer.table.id, df.id);
    }
  });

  test('TableView.addViewer(Viewer)', async () => {
    coreViewers.forEach((v) => tv.addViewer(v));
    if (Array.from(tv.viewers).length - 1 !== coreViewers.length)
      throw 'TableView.addViewer failed to attach some viewers to the table view';
  });

  after(async () => {
    tv.close();
    grok.shell.closeTable(df);
  });
});
