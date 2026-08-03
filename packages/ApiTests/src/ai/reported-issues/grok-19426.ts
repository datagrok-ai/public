import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolToggle, expectLook, expectNoThrow, findProp, look} from '../helpers';

// Regression coverage for GROK-19426: Pie chart showValue / showPercentage (renamed from showInnerPercent).
category('AI: GROK-19426: Pie chart showValue and showPercentage', () => {
  const v = (n: number = 50): DG.PieChartViewer =>
    DG.Viewer.pieChart(demog(n), {category: 'race'}) as DG.PieChartViewer;

  test('showValue round-trips through props and look across toggle steps', async () => {
    expectBoolToggle(v(), 'showValue');
  });

  test('showPercentage round-trips independently of showValue', async () => {
    const c = v();
    c.props['showValue'] = false;
    c.props['showPercentage'] = true;
    expectLook(c, {showPercentage: true, showValue: false});
    c.props['showPercentage'] = false;
    expectLook(c, {showPercentage: false, showValue: false});
    c.props['showPercentage'] = true;
    expectLook(c, {showPercentage: true, showValue: false});
  });

  test('both showValue and showPercentage true: no throw, both keys in look', async () => {
    const c = v();
    expectNoThrow(() => c.setOptions({showValue: true, showPercentage: true}));
    expectLook(c, {showValue: true, showPercentage: true});
  });

  test('getProperties lists showValue and showPercentage', async () => {
    const c = v(20);
    for (const name of ['showValue', 'showPercentage'])
      expect(findProp(c, name) != null, true);
  });

  test('legacy showInnerPercent does not surface in look (renamed to showPercentage)', async () => {
    const c = v();
    expect(look(c)['showInnerPercent'] === undefined, true);
    c.props['showPercentage'] = true;
    expectLook(c, {showPercentage: true});
    expect(look(c)['showInnerPercent'] === undefined, true);
  });
});
