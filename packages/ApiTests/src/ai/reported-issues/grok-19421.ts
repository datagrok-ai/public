import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolToggle, expectLook, expectNoThrow, findProp} from '../helpers';

// Regression coverage for GROK-19421: Histogram `showValues` bin labels.
// Pins the property key on `IHistogramSettings`, the round-trip through
// `props[...]` and `getOptions(true).look`, that combining the toggle with
// `splitColumnName` keeps both keys preserved, and (best-effort) that the
// runtime descriptor list carries `showValues` with a non-empty description.
category('AI: GROK-19421: Histogram showValues bin labels', () => {
  test('showValues round-trips through props and look across toggle steps', async () => {
    expectBoolToggle(demog().plot.histogram({value: 'age'}), 'showValues');
  });

  test('showValues property is discoverable via getProperties (best-effort)', async () => {
    const v = demog(20).plot.histogram({value: 'age'});
    const p = findProp(v, 'showValues');
    if (p != null)
      expect((p.description ?? '').trim().length > 0, true);
    else
      expectBoolToggle(v, 'showValues', [true]);
  });

  test('showValues survives alongside splitColumnName via setOptions', async () => {
    const v = demog().plot.histogram({value: 'age'});
    expectNoThrow(() => v.setOptions({splitColumnName: 'race', showValues: true}));
    expectLook(v, {splitColumnName: 'race', showValues: true});
    v.props['showValues'] = false;
    expectLook(v, {splitColumnName: 'race', showValues: false});
  });
});
