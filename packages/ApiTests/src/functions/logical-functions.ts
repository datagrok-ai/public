import {category, test} from '@datagrok-libraries/utils/src/test';
import {check} from './utils';


category('Logical functions', () => {
  test('And', () => check({
    'And(true, true)': true,
    'And(true, false)': false,
    'And(false, true)': false,
    'And(false, false)': false,
    'And(5 == 5, 10 < 20)': true,
    'And(2, 3)': null,
    'And(0, 1)': null,
    'And(1, 1)': null,
  }));

  test('Not', () => check({
    'Not(true)': false,
    'Not(false)': true,
    'Not(1)': null,
    'Not(0)': null,
  }));

  test('Or', () => check({
    'Or(true, true)': true,
    'Or(true, false)': true,
    'Or(false, true)': true,
    'Or(false, false)': false,
    'Or(5 == 6, 20 < 10)': false,
    'Or(2, 3)': null,
    'Or(0, 1)': null,
    'Or(1, 1)': null,
  }));

  test('Xor', () => check({
    'Xor(true, true)': false,
    'Xor(true, false)': true,
    'Xor(false, true)': true,
    'Xor(false, false)': false,
    'Xor(5 == 6, 20 < 10)': false,
    'Xor(5 == 5, 10 < 20)': false,
    'Xor(2, 3)': null,
    'Xor(2, 2)': null,
    'Xor(1, 0)': null,
    'Xor(1, 1)': null,
  }));
});
