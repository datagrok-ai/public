import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectCleared, expectLook, expectNoThrow, expectPropAndLook, look} from '../helpers';

// Regression coverage for GROK-19511: scatter plot's `linesBy`
// (and serialized `linesByColumnName`) is independent from `color`.
// Pure prop-existence + round-trip — no rendering assertions.
category('AI: GROK-19511: Scatter plot linesBy independent of color', () => {
  test('linesBy and color stay independent in df.plot.scatter', async () => {
    const v = demog().plot.scatter({x: 'age', y: 'height', color: 'race', linesBy: 'sex'});
    expectPropAndLook(v, {linesByColumnName: 'sex', colorColumnName: 'race'});
    expect(v.props['linesByColumnName'] !== v.props['colorColumnName'], true);
  });

  test('linesByColumnName is null/empty when linesBy not specified', async () => {
    const v = demog().plot.scatter({x: 'age', y: 'height'});
    expectCleared(v.props['linesByColumnName']);
    expectCleared(look(v)['linesByColumnName']);
  });

  test('linesByColumnName round-trips through getOptions(true).look', async () => {
    const v = DG.Viewer.scatterPlot(demog(), {x: 'age', y: 'height', color: 'race', linesBy: 'sex'});
    expectLook(v, {linesByColumnName: 'sex', colorColumnName: 'race'});
    expect(look(v)['linesByColumnName'] !== look(v)['colorColumnName'], true);
  });

  test('mutating linesByColumnName post-construction does not throw and round-trips', async () => {
    const v = DG.Viewer.scatterPlot(demog(), {x: 'age', y: 'height', color: 'race', linesBy: 'sex'});
    expectNoThrow(() => v.setOptions({linesByColumnName: 'race'}));
    expectPropAndLook(v, {linesByColumnName: 'race'});
    expectNoThrow(() => {v.props['linesByColumnName'] = 'site';});
    expectPropAndLook(v, {linesByColumnName: 'site'});
    expect(v.props['colorColumnName'], 'race');
  });
});
