import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolToggle, expectLook, expectNoThrow, findProp} from '../helpers';

// Regression coverage for GROK-19421: Histogram showValues bin labels.
category('AI: GROK-19421: Histogram showValues bin labels', () => {
  test('showValues round-trips through props and look across toggle steps', async () => {
    expectBoolToggle(demog().plot.histogram({value: 'age'}), 'showValues');
  });

  test('showValues property is discoverable via getProperties', async () => {
    const v = demog(20).plot.histogram({value: 'age'});
    expect(findProp(v, 'showValues') != null, true);
  });

  test('showValues survives alongside splitColumnName via setOptions', async () => {
    const v = demog().plot.histogram({value: 'age'});
    expectNoThrow(() => v.setOptions({splitColumnName: 'race', showValues: true}));
    expectLook(v, {splitColumnName: 'race', showValues: true});
    v.props['showValues'] = false;
    expectLook(v, {splitColumnName: 'race', showValues: false});
  });
});
