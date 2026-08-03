import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {demog, df, look, expectRoundTripPropAndLook, expectAttachesNoThrow} from '../helpers';

// DG.TileViewer

category('AI: Viewers: Tile Viewer', () => {
  test('tile factory returns typed DG.TileViewer of the right type', async () => {
    const v = DG.Viewer.tile(demog());
    expect(v instanceof DG.TileViewer, true);
    expect(v.type, DG.VIEWER.TILE_VIEWER);
  });

  test('scalar @Prop settings round-trip', async () => {
    const v = DG.Viewer.tile(demog());
    expectRoundTripPropAndLook(v,
      {allowDragBetweenLanes: false, autoGenerate: true, tilesFont: 'normal normal 15px "Roboto"'});
  });

  test('lanesColumnName + cardMarkup round-trip', async () => {
    const v = DG.Viewer.tile(demog());
    expectRoundTripPropAndLook(v, {lanesColumnName: 'sex', cardMarkup: '<div>x</div>'});
  });

  test('lanes array prop deep-checks via look', async () => {
    const v = DG.Viewer.tile(demog());
    // lanes is array-valued — compare it with a deep array check
    // (the scalar expect() in expectRoundTripPropAndLook can't compare arrays).
    v.setOptions({lanes: ['Male', 'Female']});
    expectArray(look(v).lanes, ['Male', 'Female']);
  });

  // onViewerRendered is a base-Viewer event (covered in Viewers: Lifecycle Events); the boundary
  // test below covers attach. No per-viewer-type render-event re-test here.

  test('boundary: single-row frame attaches without throwing', async () => {
    const one = df([['cat', DG.COLUMN_TYPE.STRING, ['a']], ['val', DG.COLUMN_TYPE.FLOAT, [1.0]]]);
    await expectAttachesNoThrow(one, DG.VIEWER.TILE_VIEWER, {lanesColumnName: 'cat'});
  });
}, {owner: 'agolovko@datagrok.ai'});
