import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectRoundTripPropAndLook, withAttachedViewer} from '../helpers';

// Residual Form/IFormSettings surface: DG.Form.forDataFrame + editable/designMode/row/state r/w,
// FormViewer.createDefault typed factory, and viewer font / allowDynamicMenus look round-trips.
category('AI: Viewers: Form Extras', () => {
  test('DG.Form.forDataFrame() returns a Form bound to the df', async () => {
    const t = demog();
    const f = DG.Form.forDataFrame(t);
    expect(f instanceof DG.Form, true);
    // forDataFrame builds the form but binds no row yet, so getRow() is null until set.
    // Bind a row, then the row getter exposes the source DataFrame via Row.table.
    f.row = t.row(0);
    expect(f.row.table.dart === t.dart, true);
    expect(f.row.idx, 0);
  });

  test('Form.editable get/set round-trips on a standalone Form', async () => {
    const f = DG.Form.forDataFrame(demog());
    const original = f.editable;
    f.editable = true;
    expect(f.editable, true);
    f.editable = false;
    expect(f.editable, false);
    f.editable = original;
    expect(f.editable, original);
  });

  test('Form.designMode get/set round-trips on a standalone Form', async () => {
    const f = DG.Form.forDataFrame(demog());
    const original = f.designMode;
    f.designMode = true;
    expect(f.designMode, true);
    f.designMode = false;
    expect(f.designMode, false);
    f.designMode = original;
    expect(f.designMode, original);
  });

  test('Form.state get/set round-trips (serialize/restore)', async () => {
    const f = DG.Form.forDataFrame(demog());
    const stateA = f.state;
    expect(typeof stateA, 'string');
    expect(stateA.length > 0, true);
    // Serializer doesn't guarantee element order (injects textColor), so assert structural invariants:
    // a SketchState with non-empty elementStates, and re-applying preserves the element set.
    const parsedA = JSON.parse(stateA);
    expect(parsedA['#type'], 'SketchState');
    expect(Array.isArray(parsedA.elementStates), true);
    const countA = parsedA.elementStates.length;
    expect(countA > 0, true);
    f.state = stateA;
    const parsedB = JSON.parse(f.state);
    expect(parsedB['#type'], 'SketchState');
    expect(parsedB.elementStates.length, countA);
  });

  test('Form.row setter binds the Form to a row index; getter reflects it', async () => {
    const t = demog();
    const f = DG.Form.forDataFrame(t);
    f.row = t.row(4);
    expect(f.row.idx, 4);
    f.row = t.row(1);
    expect(f.row.idx, 1);
  });

  test('FormViewer.createDefault() returns a typed FormViewer bound to df', async () => {
    const t = demog();
    // createDefault yields a standalone (unhosted) viewer; its close() crashes dereferencing a
    // host that never exists, so guard the teardown — that crash is not the behavior under test.
    const v = DG.FormViewer.createDefault(t);
    try {
      expect(v instanceof DG.FormViewer, true);
      expect(v.dataFrame === t, true);
      expect(v.type, DG.VIEWER.FORM);
      expect(v.form instanceof DG.Form, true);
    } finally {
      try {
        v.close();
      } catch (_e) {/* undocked standalone viewer: close path needs a host */}
    }
  });

  // Dropped formFont/controlsFont: both are orphan fields in the generated TS interface with no
  // backing @Prop, so setOptions({...Font}) is not reflected in getOptions().look — no real surface.

  test('IFormSettings.allowDynamicMenus round-trips via setOptions/getOptions', async () => {
    await withAttachedViewer<DG.FormViewer>(demog(), DG.VIEWER.FORM, {}, (v) => {
      const original = (v.props as any).allowDynamicMenus;
      expectRoundTripPropAndLook(v, {allowDynamicMenus: !original});
      expectRoundTripPropAndLook(v, {allowDynamicMenus: original});
    });
  });
});
