import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-2892: Pie chart `onClick='Filter'` semantics.
// Clicking a slice should filter all rows that fall into that category. The
// actual click-to-filter pipeline is UI-only and not reachable from the JS
// API, so this file pins the *property side*: the `onClick` enum round-trips
// cleanly through the factory, `setOptions`, the `props` bag, and the `look`
// snapshot, and `getProperties()` actually exposes `onClick`.
//
// `df.plot.pie` does not exist — we use the `DG.Viewer.pieChart` factory.
// `viewer.dataFrame === df` is unreliable across the Dart→JS wrapper boundary,
// so we never assert on it. Every viewer is created inside try/finally and
// closed via `tv.close()` at the end (close() on a never-attached viewer
// throws, but here every viewer is attached via the table view).

category('AI: GROK-2892: Pie chart onClick filter semantics', () => {
  test('factory + setOptions: onClick=Filter round-trips through props and look', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = DG.Viewer.pieChart(df, {category: 'race', onClick: 'Filter'});
      tv.addViewer(v);
      expect(v instanceof DG.Viewer, true);
      expect(v.type, DG.VIEWER.PIE_CHART);

      // Factory option should be reflected on the property bag.
      expect(v.props.onClick, 'Filter');

      // setOptions() with the same value must be a no-op for the property.
      v.setOptions({onClick: 'Filter'});
      expect(v.props.onClick, 'Filter');

      // The serialized `look` snapshot — what gets persisted in layouts —
      // must agree with the live property bag.
      const look = v.getOptions(true).look;
      expect(look.onClick, 'Filter');
    }
    finally {
      tv.close();
    }
  });

  test('RowGroupAction: all three values round-trip via setOptions + props + look', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = DG.Viewer.pieChart(df, {category: 'race'});
      tv.addViewer(v);

      // RowGroupAction enum (d4.ts:987-991): Select | Filter | None.
      // Each value must round-trip cleanly — set it, then read it back via
      // both the live property bag and the persisted `look` snapshot.
      for (const action of ['Select', 'Filter', 'None']) {
        v.setOptions({onClick: action});
        expect(v.props.onClick, action);
        const look = v.getOptions(true).look;
        expect(look.onClick, action);
      }
    }
    finally {
      tv.close();
    }
  });

  test('getProperties(): onClick is exposed (best-effort description check)', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = DG.Viewer.pieChart(df, {category: 'race', onClick: 'Filter'});
      tv.addViewer(v);

      const props = v.getProperties();
      expect(Array.isArray(props), true);
      const onClickProp = props.find((p) => p.name === 'onClick');
      expect(onClickProp != null, true);

      // Description is best-effort — pin only that it's a non-empty string
      // when present, without failing the test if the platform omits it.
      if (onClickProp != null && onClickProp.description != null)
        expect(typeof onClickProp.description === 'string' && onClickProp.description.length > 0, true);
    }
    finally {
      tv.close();
    }
  });
});
