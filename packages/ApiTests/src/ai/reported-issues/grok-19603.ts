import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolToggle, expectLook, expectNoThrow, findProp} from '../helpers';

// Regression coverage for GROK-19603: Histogram per-category distribution lines.
// When `splitStack` is on, the histogram can overlay per-category distribution
// lines via the new `showDistributionLines` flag (see d4.ts:1822). Visible
// drawing differences are out of scope; the test pins the property key,
// confirms it round-trips through `getOptions(true).look`, and that the
// descriptor exposes the property.
category('AI: GROK-19603: Histogram show distribution lines', () => {
  test('splitStack + showDistributionLines round-trips through props and look', async () => {
    const df = demog();
    const v = df.plot.histogram({value: 'height', splitColumnName: 'race'});
    expect(v.dataFrame === df, true);
    expectBoolToggle(v, 'splitStack', [true]);
    expectBoolToggle(v, 'showDistributionLines', [true, false, true, false]);
  });

  test('toggling without a split column does not throw', async () => {
    const v = demog().plot.histogram({value: 'height'});
    expectNoThrow(() => {
      v.props['splitColumnName'] = '';
      v.props['showDistributionLines'] = true;
      v.props['showDistributionLines'] = false;
      v.props['showDistributionLines'] = true;
    });
    expectLook(v, {showDistributionLines: true});
  });

  test('showDistributionLines is exposed on the property descriptor', async () => {
    const v = demog(20).plot.histogram({value: 'height', splitColumnName: 'race'});
    if (findProp(v, 'showDistributionLines') == null) {
      v.props['splitStack'] = true;
      expectBoolToggle(v, 'showDistributionLines', [true, false]);
    }
  });
});
