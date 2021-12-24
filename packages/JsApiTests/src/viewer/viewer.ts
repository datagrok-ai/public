import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, test} from '../test';

category('Viewer', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;

  before(async () => {
    df = grok.data.demo.demog();
    tv = grok.shell.addTableView(df);
  });

  test('TableView.addViewer', async () => {
    const coreViewers = Object.values(DG.VIEWER).filter((v) => v !== DG.VIEWER.TIMELINES && v !== DG.VIEWER.GLOBE);
    for (let viewerType of coreViewers) {
      const viewer = tv.addViewer(viewerType);
      if (!(viewer instanceof DG.Viewer))
        throw `TableView.addViewer('${viewerType}') should add a Viewer instance`;
    }
    const attachedViewers = Array.from(tv.viewers);
    if (attachedViewers.length - 1 !== coreViewers.length)
      throw 'TableView.addViewer failed to attach some viewers to the table view';
  });

  after(async () => {
    tv.close();
    grok.shell.closeTable(df);
  });
});
