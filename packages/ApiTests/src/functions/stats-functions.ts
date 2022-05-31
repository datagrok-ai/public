import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/utils/src/test';
import {check} from './utils';


category('Statistical functions', () => {
  test('Avg', () => check({
    'Avg([1, 2, 3, 4])': 2.5,
    'Avg([null, 0.3, 0.7])': 0.5,
    'Avg([-1, -2, 3])': 0,
    'Avg([1, null, 2, 3])': 2,
    'Avg([])': undefined,
    'Avg([null, null])': DG.FLOAT_NULL,
  }));

  test('Max', () => check({
    'Max([1, 2, 4, 3])': 4,
    'Max([1, null, 2, 3])': 3,
    'Max([null, 1, 0.7, 0.3])': 1,
    'Max([1.5, -2, 1.9])': 1.9,
    'Max([null, null])': DG.FLOAT_NULL,
  }));

});
