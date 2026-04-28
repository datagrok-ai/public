import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-13205: Trellis plot inner viewers (Line / Box /
// Pie / Bar) used to duplicate data instead of partitioning when rowSource was
// set to one of the per-row variants. JS-side coverage pins state-machine
// integrity only — the property bag accepts each RowSet literal, the viewer
// instance survives without throwing, and the value persists through
// getOptions(true).look. Visual partition correctness is not asserted here and
// must be checked in manual / image-diff QA.
//
// RowSet literals taken from public/js-api/src/interfaces/d4.d.ts:
//   All, Filtered, Selected, SelectedOrCurrent, FilteredSelected,
//   MouseOverGroup, CurrentRow, MouseOverRow.
category('AI: GROK-13205: Trellis row source variants', () => {
  function makeDf(): DG.DataFrame {
    const n = 60;
    const value = new Array<number>(n);
    const split1 = new Array<string>(n);
    const split2 = new Array<string>(n);
    const idx = new Array<number>(n);
    for (var i = 0; i < n; i++) {
      value[i] = (i * 7) % 23;
      split1[i] = ['a', 'b', 'c', 'd'][i % 4];
      split2[i] = ['x', 'y', 'z'][i % 3];
      idx[i] = i;
    }
    return DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'idx', idx),
      DG.Column.fromList('double', 'value', value),
      DG.Column.fromList('string', 'cat1', split1),
      DG.Column.fromList('string', 'cat2', split2),
    ]);
  }

  function applyRowSourceState(df: DG.DataFrame, rowSource: string): void {
    // Drive the corresponding df state so the inner viewers have something to
    // partition against. We don't assert on bar geometry — only that the
    // setter does not throw under each combination.
    if (rowSource === 'Selected' || rowSource === 'SelectedOrCurrent' ||
        rowSource === 'FilteredSelected') {
      df.selection.setAll(false);
      for (var i = 0; i < df.rowCount; i += 2)
        df.selection.set(i, true);
    }
    if (rowSource === 'CurrentRow' || rowSource === 'SelectedOrCurrent')
      df.currentRowIdx = 5;
    if (rowSource === 'MouseOverRow' || rowSource === 'MouseOverGroup')
      df.mouseOverRowIdx = 7;
    if (rowSource === 'FilteredSelected' || rowSource === 'Filtered') {
      df.filter.setAll(true);
      for (var i = 0; i < df.rowCount; i += 3)
        df.filter.set(i, false);
    }
  }

  const rowSources = [
    'Selected',
    'SelectedOrCurrent',
    'CurrentRow',
    'MouseOverRow',
    'FilteredSelected',
    'MouseOverGroup',
  ];

  test('rowSource literals persist through props and look', async () => {
    for (var rs of rowSources) {
      const df = makeDf();
      applyRowSourceState(df, rs);
      const v = DG.Viewer.trellisPlot(df, {
        viewerType: 'Line chart',
        xColumnNames: ['cat1'],
        yColumnNames: ['cat2'],
      });
      expect(v instanceof DG.Viewer, true);
      expect(v.type, DG.VIEWER.TRELLIS_PLOT);
      (v.props as any)['rowSource'] = rs;
      expect(v.props['rowSource'], rs);
      expect(v.getOptions(true).look['rowSource'], rs);
      // Detached viewer (never attached to a TableView) — do not call close().
    }
  });

  test('inner viewerType swap survives across rowSource variants', async () => {
    const innerTypes = ['Bar chart', 'Box plot', 'Pie chart', 'Line chart'];
    for (var rs of rowSources) {
      const df = makeDf();
      applyRowSourceState(df, rs);
      const v = DG.Viewer.trellisPlot(df, {
        viewerType: 'Line chart',
        xColumnNames: ['cat1'],
        yColumnNames: ['cat2'],
      });
      (v.props as any)['rowSource'] = rs;
      for (var inner of innerTypes) {
        v.props['viewerType'] = inner;
        expect(v.props['viewerType'], inner);
        expect(v.getOptions(true).look['viewerType'], inner);
        // rowSource must not be clobbered by the inner-type swap.
        expect(v.props['rowSource'], rs);
      }
      // Detached viewer (never attached to a TableView) — do not call close().
    }
  });
});
