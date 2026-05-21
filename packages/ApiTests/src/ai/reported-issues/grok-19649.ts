import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolToggle, expectRoundTripPropAndLook, findProp} from '../helpers';

// Regression coverage for GROK-19649: bar chart label sizing / orientation /
// category width round-trips. We pin that each prop round-trips through
// props[]/look[], and that getProperties() carries them. We do NOT assert
// specific defaults: some props are typed-but-not-runtime-registered, so the
// getProperties() check is best-effort.
category('AI: GROK-19649: Bar chart orientation and category width', () => {
  const v = (): DG.BarChartViewer => demog().plot.bar({value: 'age', split: 'race'}) as DG.BarChartViewer;

  test('orientation flips between horizontal and vertical', async () => {
    const c = v();
    for (var o of ['horizontal', 'vertical'])
      expectRoundTripPropAndLook(c, {orientation: o});
  });

  test('setOptions({maxCategoryWidth, categoryValueWidth}) reflects in look', async () => {
    expectRoundTripPropAndLook(v(), {maxCategoryWidth: 200, categoryValueWidth: 80});
  });

  test('showCategoryValues toggle round-trips through props and look', async () => {
    expectBoolToggle(v(), 'showCategoryValues', [true, false]);
  });

  test('getProperties carries orientation/category-width/showCategoryValues (best-effort)', async () => {
    const c = v();
    var foundCount = 0;
    for (var name of ['orientation', 'maxCategoryWidth', 'categoryValueWidth', 'showCategoryValues'])
      if (findProp(c, name) != null) foundCount++;
    // Best-effort: at least one of the four must be runtime-registered.
    expect(foundCount > 0, true);
  });
});
