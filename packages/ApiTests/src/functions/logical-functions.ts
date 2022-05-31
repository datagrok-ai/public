import {category, test} from '@datagrok-libraries/utils/src/test';
import {check} from './utils';


category('Logical functions', () => {
  test('And', () => check({
    'And(true, true)': true,
    'And(true, false)': false,
    'And(false, true)': false,
    'And(false, false)': false,
    'And(5 == 5, 10 < 20)': true,
    'And(2, 3)': 3,
    'And(0, 1)': 0,
    'And(1, 1)': 1,
  }));

  test('Not', () => check({
    'Not(true)': false,
    'Not(false)': true,
    'Not(1)': false,
    'Not(0)': true,
  }));

  test('Or', () => check({
    'Or(true, true)': true,
    'Or(true, false)': true,
    'Or(false, true)': true,
    'Or(false, false)': false,
    'Or(5 == 6, 20 < 10)': false,
    'Or(2, 3)': 2,
    'Or(0, 1)': 1,
    'Or(1, 1)': 1,
  }));

  test('Xor', () => check({
    'Xor(true, true)': false,
    'Xor(true, false)': true,
    'Xor(false, true)': true,
    'Xor(false, false)': false,
    'Xor(5 == 6, 20 < 10)': false,
    'Xor(5 == 5, 10 < 20)': false,
    'Xor(2, 3)': false,
    'Xor(2, 2)': false,
    'Xor(1, 0)': true,
    'Xor(1, 1)': false,
  }));
});
