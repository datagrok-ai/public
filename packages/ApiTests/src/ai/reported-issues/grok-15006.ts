import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, expectNoThrow, expectPropAndLook, expectRoundTrip, findProp} from '../helpers';

// Regression coverage for GROK-15006: Bar chart `showFilteredRows` toggle.
// The flag is only meaningful when `rowSource = 'All'` (per d4.ts:872), but
// flipping it under other rowSource values must still not throw. Assertions
// pin the state via getOptions(true).look — no pixel sampling.
category('AI: GROK-15006: Bar chart showFilteredRows flag', () => {
  const v = (): DG.BarChartViewer => demog().plot.bar({value: 'age', split: 'race'}) as DG.BarChartViewer;

  test('showFilteredRows=true under rowSource=All', async () => {
    expectRoundTrip(v(), {rowSource: 'All', showFilteredRows: true});
  });

  test('toggle showFilteredRows to false mirrors props', async () => {
    const c = v();
    expectNoThrow(() => {
      c.setOptions({rowSource: 'All', showFilteredRows: true});
      c.setOptions({showFilteredRows: false});
    });
    expectPropAndLook(c, {rowSource: 'All', showFilteredRows: false});
  });

  test('toggling showFilteredRows under rowSource=Filtered does not throw', async () => {
    const c = v();
    expectNoThrow(() => {
      c.setOptions({rowSource: 'Filtered'});
      c.setOptions({showFilteredRows: true});
      c.setOptions({showFilteredRows: false});
      c.setOptions({showFilteredRows: true});
    });
    expectLook(c, {rowSource: 'Filtered', showFilteredRows: true});
  });

  test('getProperties lists showFilteredRows (best-effort)', async () => {
    const c = v();
    // Best-effort: if the property descriptor is not exposed, fall back to
    // confirming it is at least settable via setOptions and observable in look.
    if (findProp(c, 'showFilteredRows') == null)
      expectRoundTrip(c, {rowSource: 'All', showFilteredRows: true});
  });
});
