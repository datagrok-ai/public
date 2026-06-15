import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolToggle, expectRoundTripPropAndLook, findProp} from '../helpers';

// Regression coverage for GROK-19649: bar chart orientation and category width round-trips.
category('AI: GROK-19649: Bar chart orientation and category width', () => {
  const v = (): DG.BarChartViewer => demog().plot.bar({value: 'age', split: 'race'}) as DG.BarChartViewer;

  test('orientation flips between horizontal and vertical', async () => {
    const c = v();
    for (const o of ['horizontal', 'vertical'])
      expectRoundTripPropAndLook(c, {orientation: o});
  });

  test('setOptions({maxCategoryWidth, categoryValueWidth}) reflects in look', async () => {
    expectRoundTripPropAndLook(v(), {maxCategoryWidth: 200, categoryValueWidth: 80});
  });

  test('showCategoryValues toggle round-trips through props and look', async () => {
    expectBoolToggle(v(), 'showCategoryValues', [true, false]);
  });

  test('getProperties carries orientation/category-width/showCategoryValues', async () => {
    const c = v();
    for (const name of ['orientation', 'maxCategoryWidth', 'categoryValueWidth', 'showCategoryValues'])
      expect(findProp(c, name) != null, true);
  });
});
