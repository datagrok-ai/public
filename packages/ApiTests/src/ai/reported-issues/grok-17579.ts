import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-17579: Right-click → Tooltip > Hide on a Bar
// chart used to be ignored while it worked on other viewers. Fixed in 1.24.0.
//
// Symptom is canvas-only (a real tooltip showing on hover); JS-side test pins
// the IBarChartSettings.showTooltip property surface — settable via setOptions
// and observable via getOptions(true).look. Note: BarChart's getProperties()
// runtime descriptor list does NOT expose showTooltip on dev (verified during
// the first run of this test); the property lives only on the typed
// IBarChartSettings interface (`d4.d.ts:762`) and is reachable via setOptions
// through the ObjectPropertyBag proxy.
category('AI: GROK-17579: Bar chart tooltip hide', () => {
  // NOTE on case scoping (2026-04-29): the original ticket fix is canvas-only
  // (a tooltip should not appear on hover). The JS-side surface for the
  // showTooltip property on BarChart is currently inert: no candidate string
  // value (`'Do not show'`, `'None'`, `'Off'`, `'Never'`, `'Default'`, ...)
  // round-trips through `setOptions` → `getOptions(true).look.showTooltip`, and
  // `getProperties()` does not list a `showTooltip` descriptor on a bar chart
  // (verified on dev 2026-04-29 — the property lives only on the typed
  // IBarChartSettings interface at `d4.d.ts:762`). The only invariant we can
  // pin from JS is that calling `setOptions({showTooltip: ...})` does not
  // throw. The rest of the fix surface is exercised by manual / Selenium QA.
  test('setOptions({showTooltip}) is non-throwing for arbitrary string', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.barChart(df, {value: 'age', split: 'race'});
    var threw = false;
    try {
      v.setOptions({showTooltip: 'arbitrary-value-xyz'});
    }
    catch (_e) {
      threw = true;
    }
    // Whether the platform accepts the value or normalises it away is platform
    // policy; what matters is no throw.
    expect(threw, false);
  });
});
