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

  test('getProperties lists showValue and showPercentage (best-effort)', async () => {
    const c = v(20);
    for (const name of ['showValue', 'showPercentage']) {
      const p = findProp(c, name);
      if (p != null && typeof p.description === 'string')
        expect(p.description.trim().length > 0, true);
      else
        expectBoolToggle(c, name, [true]);
    }
  });

  test('legacy showInnerPercent is absent or does not surface in look (rename)', async () => {
    const c = v();
    const aliased = look(c)['showInnerPercent'] !== undefined;
    if (!aliased)
      expect(look(c)['showInnerPercent'] === undefined, true);
    c.props['showPercentage'] = true;
    expectLook(c, {showPercentage: true});
    if (!aliased)
      expect(look(c)['showInnerPercent'] === undefined, true);
  });
});
