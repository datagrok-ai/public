import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {df, expectPropAndLook} from '../helpers';

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
  const ROW_SOURCES = [
    'Selected', 'SelectedOrCurrent', 'CurrentRow', 'MouseOverRow', 'FilteredSelected', 'MouseOverGroup',
  ];

  function makeDf(): DG.DataFrame {
    const n = 60;
    const value = new Array<number>(n);
    const split1 = new Array<string>(n);
    const split2 = new Array<string>(n);
    const idx = new Array<number>(n);
    for (let i = 0; i < n; i++) {
      value[i] = (i * 7) % 23;
      split1[i] = ['a', 'b', 'c', 'd'][i % 4];
      split2[i] = ['x', 'y', 'z'][i % 3];
      idx[i] = i;
    }
    return df([
      ['idx', 'int', idx], ['value', 'double', value],
      ['cat1', 'string', split1], ['cat2', 'string', split2],
    ]);
  }

  function applyRowSourceState(d: DG.DataFrame, rs: string): void {
    if (rs === 'Selected' || rs === 'SelectedOrCurrent' || rs === 'FilteredSelected') {
      d.selection.setAll(false);
      for (let i = 0; i < d.rowCount; i += 2) d.selection.set(i, true);
    }
    if (rs === 'CurrentRow' || rs === 'SelectedOrCurrent')
      d.currentRowIdx = 5;
    if (rs === 'MouseOverRow' || rs === 'MouseOverGroup')
      d.mouseOverRowIdx = 7;
    if (rs === 'FilteredSelected' || rs === 'Filtered') {
      d.filter.setAll(true);
      for (let i = 0; i < d.rowCount; i += 3) d.filter.set(i, false);
    }
  }

  function makeTrellis(d: DG.DataFrame): DG.Viewer {
    return DG.Viewer.trellisPlot(d, {viewerType: 'Line chart', xColumnNames: ['cat1'], yColumnNames: ['cat2']});
  }

  test('rowSource literals persist through props and look', async () => {
    for (const rs of ROW_SOURCES) {
      const d = makeDf();
      applyRowSourceState(d, rs);
      const v = makeTrellis(d);
      expect(v.type, DG.VIEWER.TRELLIS_PLOT);
      (v.props as any)['rowSource'] = rs;
      expectPropAndLook(v, {rowSource: rs});
    }
  });

  test('inner viewerType swap survives across rowSource variants', async () => {
    const innerTypes = ['Bar chart', 'Box plot', 'Pie chart', 'Line chart'];
    for (const rs of ROW_SOURCES) {
      const d = makeDf();
      applyRowSourceState(d, rs);
      const v = makeTrellis(d);
      (v.props as any)['rowSource'] = rs;
      for (const inner of innerTypes) {
        v.props['viewerType'] = inner;
        expectPropAndLook(v, {viewerType: inner});
        expect(v.props['rowSource'], rs);
      }
    }
  });
});
