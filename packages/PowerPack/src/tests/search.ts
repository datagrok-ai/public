import {category, test} from '@datagrok-libraries/test/src/test';
import {initTemplates} from '../search/templates-search';

category('Search', () => {
  test('InitTemplates', async () => await initTemplates());
});
