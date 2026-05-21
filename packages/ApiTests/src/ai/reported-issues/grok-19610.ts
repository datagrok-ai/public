import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolToggle, expectNoThrow, findProp, look} from '../helpers';

// Regression coverage for GROK-19610: Box plot 'Zoom Values by Filter'.
// The user-facing fix was a tooltip clarification on the property; the test
// pins the property key (`zoomValuesByFilter`), confirms it round-trips
// through getOptions(true).look, that toggling under an active filter on
// the value column does not throw, and that the property carries a
// non-empty description (the actual remediation in the ticket).
category('AI: GROK-19610: Box plot zoom-values-by-filter toggle', () => {
  const v = (): DG.Viewer => demog().plot.box({value: 'height', category: 'race'});

  test('zoomValuesByFilter round-trips through props and look', async () => {
    expectBoolToggle(v(), 'zoomValuesByFilter', [true, false]);
  });

  test('toggling under an active filter on the value column does not throw', async () => {
    const df = demog();
    const c = df.plot.box({value: 'height', category: 'race'});
    expectNoThrow(() => {
      c.props['zoomValuesByFilter'] = true;
      for (var i = 0; i < df.rowCount; i++)
        df.filter.set(i, i % 2 === 0);
      c.props['zoomValuesByFilter'] = false;
      c.props['zoomValuesByFilter'] = true;
      df.filter.setAll(true);
      c.props['zoomValuesByFilter'] = false;
    });
    expect(look(c)['zoomValuesByFilter'], false);
  });

  test('zoomValuesByFilter property exposes a non-empty description', async () => {
    const p = findProp(demog(20).plot.box({value: 'height', category: 'race'}), 'zoomValuesByFilter');
    expect(p != null, true);
    expect((p!.description ?? '').trim().length > 0, true);
  });
});
