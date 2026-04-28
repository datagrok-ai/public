import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-19511: scatter plot's `linesBy`
// (and serialized `linesByColumnName`) is independent from `color`.
// Pure prop-existence + round-trip — no rendering assertions.
category('AI: GROK-19511: Scatter plot linesBy independent of color', () => {
  test('linesBy and color stay independent in df.plot.scatter', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({x: 'age', y: 'height', color: 'race', linesBy: 'sex'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.SCATTER_PLOT);
    expect(v.props['linesByColumnName'], 'sex');
    expect(v.props['colorColumnName'], 'race');
    expect(v.props['linesByColumnName'] !== v.props['colorColumnName'], true);
    const look = v.getOptions(true).look;
    expect(look['linesByColumnName'], 'sex');
    expect(look['colorColumnName'], 'race');
  });

  test('linesByColumnName is null/empty when linesBy not specified', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({x: 'age', y: 'height'});
    expect(v.type, DG.VIEWER.SCATTER_PLOT);
    const lbcn = v.props['linesByColumnName'];
    expect(lbcn == null || lbcn === '', true);
    const look = v.getOptions(true).look;
    const lookLbcn = look['linesByColumnName'];
    expect(lookLbcn == null || lookLbcn === '', true);
  });

  test('linesByColumnName round-trips through getOptions(true).look', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.scatterPlot(df, {x: 'age', y: 'height', color: 'race', linesBy: 'sex'});
    expect(v.props['linesByColumnName'], 'sex');
    const look = v.getOptions(true).look;
    expect(look['linesByColumnName'], 'sex');
    expect(look['colorColumnName'], 'race');
    // independence preserved through serialization
    expect(look['linesByColumnName'] !== look['colorColumnName'], true);
  });

  test('mutating linesByColumnName post-construction does not throw and round-trips', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.scatterPlot(df, {x: 'age', y: 'height', color: 'race', linesBy: 'sex'});
    var threw = false;
    try {
      v.setOptions({linesByColumnName: 'race'});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    expect(v.props['linesByColumnName'], 'race');
    expect(v.getOptions(true).look['linesByColumnName'], 'race');

    // direct prop assignment as alternative path
    var threw2 = false;
    try {
      v.props['linesByColumnName'] = 'site';
    }
    catch (_e) {
      threw2 = true;
    }
    expect(threw2, false);
    expect(v.props['linesByColumnName'], 'site');
    expect(v.getOptions(true).look['linesByColumnName'], 'site');
    // color remains unchanged
    expect(v.props['colorColumnName'], 'race');
  });
});
