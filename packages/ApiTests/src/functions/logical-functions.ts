import {category, test} from '@datagrok-libraries/test/src/test';
import {check} from './utils';


category('Functions: Logical', () => {
  test('And', () => check({
    'And(true, true)': true,
    'And(true, false)': false,
    'And(false, true)': false,
    'And(false, false)': false,
    'And(5 == 5, 10 < 20)': true,
    'And(2, 3)': undefined,
    'And(0, 1)': undefined,
    'And(1, 1)': undefined,
  }));

  test('Not', () => check({
    'Not(true)': false,
    'Not(false)': true,
    'Not(1)': undefined,
    'Not(0)': undefined,
  }));

  test('Or', () => check({
    'Or(true, true)': true,
    'Or(true, false)': true,
    'Or(false, true)': true,
    'Or(false, false)': false,
    'Or(5 == 6, 20 < 10)': false,
    'Or(2, 3)': undefined,
    'Or(0, 1)': undefined,
    'Or(1, 1)': undefined,
  }));

  test('Xor', () => check({
    'Xor(true, true)': false,
    'Xor(true, false)': true,
    'Xor(false, true)': true,
    'Xor(false, false)': false,
    'Xor(5 == 6, 20 < 10)': false,
    'Xor(5 == 5, 10 < 20)': false,
    'Xor(2, 3)': undefined,
    'Xor(2, 2)': undefined,
    'Xor(1, 0)': undefined,
    'Xor(1, 1)': undefined,
  }));
});
