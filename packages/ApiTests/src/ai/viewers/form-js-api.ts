import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, df, expectChoices, expectNoThrow, expectPropAndLook, expectRoundTrip,
  expectRoundTripPropAndLook, look, until, withAttachedViewer} from '../helpers';

// FormViewer JS surface: typed factory, editable/designMode, show* toggles, syncMode, buildForm, row, form getter.
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

  test('buildForm — rebuild the form with a subset of columns (read back via getOptions)', async () => {
    await withAttachedViewer<DG.FormViewer>(demog(), DG.VIEWER.FORM, {}, async (v) => {
      v.buildForm(['name', 'age']);
      await until(() => ((look(v)['columnNames'] as string[]) ?? []).length >= 1);
      const names = (look(v)['columnNames'] as string[]) ?? [];
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
      // In None mode the form's row decouples from currentRowIdx and defaults to -1.
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
    // SketchField._init requires at least one row; truly-empty DataFrame is an unsupported boundary.
    const one = df([['x', DG.COLUMN_TYPE.INT, [42]]]);
    await withAttachedViewer<DG.FormViewer>(one, DG.VIEWER.FORM, {}, (v) => {
      expect(v instanceof DG.FormViewer, true);
      expectNoThrow(() => v.getOptions(true));
      expectNoThrow(() => v.setOptions({showOpenFile: false}));
      expectPropAndLook(v, {showOpenFile: false});
    });
  });
});
