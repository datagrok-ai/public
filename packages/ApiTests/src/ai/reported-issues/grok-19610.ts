import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-19610: Box plot 'Zoom Values by Filter'.
// The user-facing fix was a tooltip clarification on the property; the test
// pins the property key (`zoomValuesByFilter`), confirms it round-trips
// through getOptions(true).look, that toggling under an active filter on
// the value column does not throw, and that the property carries a
// non-empty description (the actual remediation in the ticket).
category('AI: GROK-19610: Box plot zoom-values-by-filter toggle', () => {
  test('zoomValuesByFilter round-trips through props and look', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.box({value: 'height', category: 'race'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.BOX_PLOT);
    expect(v.dataFrame === df, true);

    v.props['zoomValuesByFilter'] = true;
    expect(v.props['zoomValuesByFilter'], true);
    expect(v.getOptions(true).look['zoomValuesByFilter'], true);

    v.props['zoomValuesByFilter'] = false;
    expect(v.props['zoomValuesByFilter'], false);
    expect(v.getOptions(true).look['zoomValuesByFilter'], false);
  });

  test('toggling under an active filter on the value column does not throw', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.box({value: 'height', category: 'race'});
    var threw = false;
    try {
      v.props['zoomValuesByFilter'] = true;
      // Filter half of the rows out; rows are filtered (not rows of a particular column),
      // but the value column is `height`, exercising the same code path the bug report hits.
      for (var i = 0; i < df.rowCount; i++)
        df.filter.set(i, i % 2 === 0);
      v.props['zoomValuesByFilter'] = false;
      v.props['zoomValuesByFilter'] = true;
      // Reset the filter and toggle once more.
      df.filter.setAll(true);
      v.props['zoomValuesByFilter'] = false;
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    expect(v.getOptions(true).look['zoomValuesByFilter'], false);
  });

  test('zoomValuesByFilter property exposes a non-empty description', async () => {
    const df = grok.data.demo.demog(20);
    const v = df.plot.box({value: 'height', category: 'race'});
    const props = v.getProperties() as DG.Property[];
    expect(Array.isArray(props), true);
    expect(props.length > 0, true);
    var found: DG.Property | null = null;
    for (var p of props) {
      if (p.name === 'zoomValuesByFilter') {
        found = p;
        break;
      }
    }
    expect(found != null, true);
    const desc = found!.description;
    expect(typeof desc === 'string', true);
    expect((desc ?? '').trim().length > 0, true);
  });
});
