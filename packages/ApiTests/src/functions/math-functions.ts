import {category, test} from '@datagrok-libraries/utils/src/test';
import {check} from './utils';


category('Math functions', () => {
  test('Add', () => check({
    'Add(2, 10)': 12,
  }));
});
