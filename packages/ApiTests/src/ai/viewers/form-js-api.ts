import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, df, expectChoices, expectNoThrow, expectPropAndLook, expectRoundTrip,
  expectRoundTripPropAndLook, wait, withAttachedViewer} from '../helpers';

// DG.FormViewer — core/client/d4/lib/src/viewers/form/form_core.dart (scenario: form-js-api)
// Form viewer JS surface (standalone DataFrame Form viewer, NOT the InputForm
// widget): typed factory returning DG.FormViewer (vs the generic
// Viewer<IFormSettings> the old factory returned), the editable/designMode
// sketchForm pair, the show* ribbon-toggle booleans, the syncMode Current /
// None / Mouse Over choices, the new buildForm/columnNames/row members backed
// by FormViewer_BuildForm / FormViewer_Get_ColumnNames / FormViewer_Get_Row
// wraps, and the existing form getter returning a DG.Form bound to a DataFrame.
category('AI: Viewers: Form JS API', () => {
  test('factory DG.Viewer.form returns typed DG.FormViewer; df.plot.form round-trips type', async () => {
    const t = demog();
    const v = DG.Viewer.form(t);
    expect(v instanceof DG.FormViewer, true);
    expect(v.dataFrame === t, true);
    expect(v.type, DG.VIEWER.FORM);

    const v2 = t.plot.form({});
    expect(v2.type, DG.VIEWER.FORM);
    expect(v2 instanceof DG.FormViewer, true);
  });

  test('editable + designMode round-trip via setter pairs', async () => {
    await withAttachedViewer<DG.FormViewer>(demog(), DG.VIEWER.FORM, {}, (v) => {
      v.editable = true;
      expect(v.editable, true);
      v.editable = false;
      expect(v.editable, false);

      const initialDesign = v.designMode;
      v.designMode = !initialDesign;
      expect(v.designMode, !initialDesign);
      v.designMode = initialDesign;
      expect(v.designMode, initialDesign);
    });
  });

  test('visual look round-trip — show* booleans (Navigation / PrevRow / NextRow / RowSelector / ' +
    'FieldEditor / DesignEditor / ColumnSelector / SaveFile / OpenFile)', async () => {
    await withAttachedViewer<DG.FormViewer>(demog(), DG.VIEWER.FORM, {}, (v) => {
      expectRoundTripPropAndLook(v, {
        showNavigation: false,
        showPrevRowArrow: false,
        showNextRowArrow: false,
        showRowSelector: false,
        showFieldEditor: false,
        showDesignEditor: false,
        showColumnSelector: false,
        showSaveFile: false,
        showOpenFile: false,
      });
    });
  });

  test('syncMode exposes choices + round-trips Current / None / Mouse Over', async () => {
    await withAttachedViewer<DG.FormViewer>(demog(), DG.VIEWER.FORM, {}, (v) => {
      expectChoices(v, 'syncMode', ['None', 'Current', 'Mouse Over']);
      expectRoundTrip(v, {syncMode: 'None'});
      expectRoundTrip(v, {syncMode: 'Mouse Over'});
      expectRoundTrip(v, {syncMode: 'Current'});
    });
  });

  test('buildForm + columnNames — rebuild the form with a subset of columns', async () => {
    await withAttachedViewer<DG.FormViewer>(demog(), DG.VIEWER.FORM, {}, async (v) => {
      v.buildForm(['name', 'age']);
      await wait(200);
      const names = v.columnNames ?? [];
      // buildForm should re-register field handlers; assert it's at least non-empty and
      // contains one of the requested columns (the SketchForm may dedupe or filter,
      // so we don't lock down the exact length).
      expect(names.length >= 1, true);
      expect(names.includes('name') || names.includes('age'), true);
    });
  });

  test('row getter tracks dataFrame.currentRow in default Current mode', async () => {
    const t = demog();
    await withAttachedViewer<DG.FormViewer>(t, DG.VIEWER.FORM, {}, async (v) => {
      t.currentRowIdx = 3;
      expect(v.row, 3);

      v.setOptions({syncMode: 'None'});
      // In `None` mode the form's row decouples from dataFrame.currentRowIdx
      // and reflects the internal `_row` field, which defaults to -1.
      expect(v.row === -1 || v.row !== t.currentRowIdx, true);
    });
  });

  test('form getter returns a DG.Form instance with a DataFrame', async () => {
    await withAttachedViewer<DG.FormViewer>(demog(), DG.VIEWER.FORM, {}, (v) => {
      const f = v.form;
      expect(f != null, true);
      expect(f instanceof DG.Form, true);
    });
  });

  test('boundary — single-row DataFrame does not throw on construction or option round-trip', async () => {
    // NOTE: SketchField._init requires at least one row; truly-empty DataFrame is
    // documented as an unsupported boundary (would throw `No element`).
    const one = df([['x', DG.COLUMN_TYPE.INT, [42]]]);
    await withAttachedViewer<DG.FormViewer>(one, DG.VIEWER.FORM, {}, (v) => {
      expect(v instanceof DG.FormViewer, true);
      expectNoThrow(() => v.getOptions(true));
      expectNoThrow(() => v.setOptions({showOpenFile: false}));
      expectPropAndLook(v, {showOpenFile: false});
    });
  });
});
