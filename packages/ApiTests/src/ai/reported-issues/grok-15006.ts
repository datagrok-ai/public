import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, expectNoThrow, expectPropAndLook, expectRoundTrip, findProp} from '../helpers';

// Regression coverage for GROK-15006: Bar chart showFilteredRows toggle.
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

  test('getProperties lists showFilteredRows', async () => {
    expect(findProp(v(), 'showFilteredRows') != null, true);
  });
});
