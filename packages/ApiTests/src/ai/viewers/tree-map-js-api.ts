import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {demog, df, look, expectRoundTripPropAndLook, expectAttachesNoThrow} from '../helpers';

// DG.TreeMap

category('AI: Viewers: Tree Map', () => {
  test('treeMap factory returns typed DG.TreeMap of the right type', async () => {
    const v = DG.Viewer.treeMap(demog(), {splitByColumnNames: ['sex']});
    expect(v instanceof DG.TreeMap, true);
    expect(v.type, DG.VIEWER.TREE_MAP);
  });

  test('split/color/size column props round-trip', async () => {
    const v = DG.Viewer.treeMap(demog());
    // Scalar column props go through the prop+look round-trip helper.
    expectRoundTripPropAndLook(v, {colorColumnName: 'age', sizeColumnName: 'height'});
    // splitByColumnNames is array-valued — compare it with a deep array check
    // (the scalar expect() in expectPropAndLook can't compare arrays).
    v.setOptions({splitByColumnNames: ['sex']});
    expectArray(look(v).splitByColumnNames, ['sex']);
  });

  test('aggregation and style props round-trip', async () => {
    const v = DG.Viewer.treeMap(demog());
    expectRoundTripPropAndLook(v,
      {colorAggrType: 'min', sizeAggrType: 'max', autoLayout: false, showColumnSelectionPanel: false});
  });

  // onViewerRendered is a base-Viewer event (covered in Viewers: Lifecycle Events); the boundary
  // test below covers attach. No per-viewer-type render-event re-test here.

  test('boundary: single-row frame and a no-categorical-column frame both attach without throwing', async () => {
    const one = df([['cat', DG.COLUMN_TYPE.STRING, ['a']], ['val', DG.COLUMN_TYPE.FLOAT, [1.0]]]);
    await expectAttachesNoThrow(one, DG.VIEWER.TREE_MAP, {splitByColumnNames: ['cat'], colorColumnName: 'val'});
    const numeric = df([['a', DG.COLUMN_TYPE.FLOAT, [1, 2, 3]], ['b', DG.COLUMN_TYPE.FLOAT, [4, 5, 6]]]);
    await expectAttachesNoThrow(numeric, DG.VIEWER.TREE_MAP, {});
  });
}, {owner: 'agolovko@datagrok.ai'});
