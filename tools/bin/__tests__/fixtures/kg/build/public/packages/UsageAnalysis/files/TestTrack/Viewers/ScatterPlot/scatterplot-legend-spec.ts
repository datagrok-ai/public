import {test} from '../../shared-page';
import {specTestOptions} from '../../spec-login';
test.use(specTestOptions);

test.describe.serial('Scatter plot legend', () => {
  test('color legend', async ({page}) => {
    test.setTimeout(600_000);
  });

  test('marker legend', async ({page}) => {});
});
