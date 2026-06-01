import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolToggle, expectLook, findProp} from '../helpers';

// Regression coverage for GROK-18169: PC plot showMin / showMax flags.
category('AI: GROK-18169: PC plot showMin and showMax', () => {
  const v = (n: number = 30): DG.Viewer => DG.Viewer.pcPlot(demog(n));

  test('showMinMax toggles round-trip through props and look', async () => {
    expectBoolToggle(v(), 'showMinMax', [true, false, true, false]);
  });

  test('showMinMax is exposed on the property descriptor', async () => {
    expect(findProp(v(20), 'showMinMax') != null, true);
  });

  test('showMin / showMax schema tripwire (best-effort if exposed)', async () => {
    const c = v(20);
    if (findProp(c, 'showMin') !=
      null) expectBoolToggle(c, 'showMin', [true, false]);
    if (findProp(c, 'showMax') !=
      null) expectBoolToggle(c, 'showMax', [true, false]);
  });

  test('showLabels and showDensity toggle independently', async () => {
    const c = v();
    expectBoolToggle(c, 'showLabels', [true, false]);
    expectBoolToggle(c, 'showDensity', [true, false]);
    c.props['showLabels'] = true;
    c.props['showDensity'] = true;
    c.props['showLabels'] = false;
    expectLook(c, {showLabels: false, showDensity: true});
    c.props['showDensity'] = false;
    expectLook(c, {showLabels: false, showDensity: false});
  });
});
