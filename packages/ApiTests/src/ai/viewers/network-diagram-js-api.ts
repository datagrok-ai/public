import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, df, expectNoThrow, expectRoundTripPropAndLook,
  subscribeAll, withAttachedViewer} from '../helpers';

// NetworkDiagram JS API: typed factory, events, settings round-trip, edge DataFrame.
category('AI: Viewers: NetworkDiagram JS API', () => {
  test('factory DG.Viewer.network returns typed DG.NetworkDiagramViewer; round-trips type', async () => {
    const t = demog();
    const v = DG.Viewer.network(t);
    expect(v instanceof DG.NetworkDiagramViewer, true);
    expect(v.dataFrame === t, true);
    expect(v.type, DG.VIEWER.NETWORK_DIAGRAM);
  });

  test('view.addViewer attaches a typed DG.NetworkDiagramViewer', async () => {
    await withAttachedViewer<DG.NetworkDiagramViewer>(demog(), DG.VIEWER.NETWORK_DIAGRAM, {}, (v, tv) => {
      expect(v instanceof DG.NetworkDiagramViewer, true);
      const found = tv.viewers.find((x) => x.type === DG.VIEWER.NETWORK_DIAGRAM);
      expect(found instanceof DG.NetworkDiagramViewer, true);
    });
  });

  test('onNodeClicked / onEdgeClicked are rxjs Observables', async () => {
    const v = DG.Viewer.network(demog(20));
    subscribeAll([v.onNodeClicked, v.onEdgeClicked])();
  });

  test('visual look round-trip — autoLayout, mergeNodes, click selection, column selectors, edge width', async () => {
    const v = DG.Viewer.network(demog());
    expectRoundTripPropAndLook(v, {
      autoLayout: false,
      mergeNodes: false,
      selectRowsOnClick: false,
      selectEdgesOnClick: false,
      showColumnSelectors: false,
      edgeWidth: 5,
    });
  });

  test('edge-shaped DataFrame with explicit node columns does not throw', async () => {
    const t = df([
      ['from', DG.COLUMN_TYPE.STRING, ['a', 'b', 'a']],
      ['to', DG.COLUMN_TYPE.STRING, ['b', 'c', 'c']],
    ]);
    let v: DG.NetworkDiagramViewer;
    expectNoThrow(() => {
      v = DG.Viewer.network(t, {node1ColumnName: 'from', node2ColumnName: 'to'});
    });
    expectNoThrow(() => v!.getOptions(true));
  });
});
