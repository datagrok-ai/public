import {category, test} from '@datagrok-libraries/test/src/test';
import {check, checkRandomInt} from './utils';


category('Functions: Math', () => {
  test('Abs', () => check({
    'Abs(-10)': 10,
    'Abs(10)': 10,
    'Abs(0)': 0,
  }));

  test('Acos', () => check({
    'Acos(-1)': 3.142,
    'Acos(0)': 1.571,
    'Acos(0.5)': 1.047,
    'Acos(1)': 0,
  }));

  test('Add', () => check({
    'Add(2, 10)': 12,
    'Add(24.06, 100.101)': 124.161,
    'Add(-70, 25)': -45,
    'Add(-11, -22)': -33,
    'Add(999, 0)': 999,
    'Add(0, -12)': -12,
  }));

  test('Asin', () => check({
    'Asin(-1)': -1.571,
    'Asin(0)': 0,
    'Asin(0.5)': 0.524,
    'Asin(1)': 1.571,
  }));

  test('Atan', () => check({
    'Atan(0)': 0,
    'Atan(1)': 0.785,
    'Atan(1.5)': 0.983,
    'Atan(-1)': -0.785,
  }));

  test('Atan2', () => check({
    'Atan2(1, 0)': 1.5707,
    // 'Atan2(0, -0)': 0, don't want to test undefined behavior in different browsers
    'Atan2(0, 1)': 0,
    'Atan2(-1, 2)': -0.464,
    'Atan2(1, -2)': 2.678,
  }));

  test('Ceil', () => check({
    'Ceil(3.5)': 4,
    'Ceil(-3.5)': -3,
    'Ceil(1.1)': 2,
    'Ceil(-1.1)': -1,
    'Ceil(10)': 10,
  }));

  test('Cos', () => check({
    'Cos(0)': 1,
    'Cos(PI / 6)': 0.866,
    'Cos(PI / 4)': 0.707,
    'Cos(PI / 3)': 0.5,
    'Cos(PI / 2)': 0.0001, // enables float comparison
  }));

  test('Div', () => check({
    'Div(7.542, 2)': 3.771,
    'Div(-12, 3)': -4,
    'Div(4321, 1)': 4321,
    'Div(0, 4)': 0,
    'Div(1, 0)': Infinity,
  }));

  test('Eq', () => check({
    'Eq(5, 5)': true,
    'Eq(-1, 1)': false,
    'Eq(-0, 0)': true,
    'Eq(3.142, 3.142)': true,
    'Eq("1", 1)': false,
    'Eq(0, null)': false,
    'Eq(null, null)': true,
  }));

  test('Exp', () => check({
    'Exp(1)': 2.718,
    'Exp(2)': 7.389,
    'Eq(Exp(1), E)': true,
  }));

  test('Floor', () => check({
    'Floor(3.5)': 3,
    'Floor(-3.5)': -4,
    'Floor(1.1)': 1,
    'Floor(-1.1)': -2,
    'Floor(10)': 10,
  }));

  test('Greater', () => check({
    'Greater(5, 4)': true,
    'Greater(-5, -10)': true,
    'Greater(5, 5)': false,
    'Greater(4, 5)': false,
  }));

  test('If', () => check({
    'If(true, "a", "b")': 'a',
    'If(false, "a", "b")': 'b',
    'If(true, If(true, "a", "b"), "c")': 'a',
    'If(false, "a", If(false, "b", "c"))': 'c',
    'If(Eq(10, 10), 1, 0)': 1,
    'If(Eq(10, 50), 1, 0)': 0,
    'If(Boolean(1), Boolean(1), Boolean(0))': true,
    'If(Boolean(0), Boolean(1), Boolean(0))': false,
  }));

  test('Ln', () => check({
    'Ln(1)': 0,
    'Ln(E)': 1,
    'Ln(5)': 1.609,
  }));

  test('Log', () => check({
    'Log(25, 5)': 2,
    'Log(100, 10)': 2,
    'Log(2, 10)': 0.301,
  }));

  test('Log10', () => check({
    'Log10(1)': 0,
    'Log10(10)': 1,
    'Log10(100)': 2,
    'Log10(5)': 0.699,
    'Log10(32)': 1.505,
  }));

  test('Mod', () => check({
    'Mod(8, 3)': 2,
    'Mod(9, 3)': 0,
    'Mod(1, 3)': 1,
    'Mod(7, 7)': 0,
  }));

  test('Mul', () => check({
    'Mul(10, 1.5)': 15,
    'Mul(12, 30)': 360,
    'Mul(0, 13.27)': 0,
    'Mul(123, 1)': 123,
  }));

  test('Neg', () => check({
    'Neg(-5)': 5,
    'Neg(11)': -11,
    'Neg(-0)': 0,
    'Neg(+1)': -1,
  }));

  test('NotEq', () => check({
    'NotEq(5, 5)': false,
    'NotEq(-1, 1)': true,
    'NotEq(-0, 0)': false,
    'NotEq(3.142, 3.142)': false,
  }));

  test('NotGreater', () => check({
    'NotGreater(4, 5)': true,
    'NotGreater(5, 5)': true,
    'NotGreater(6, 5)': false,
    'NotGreater(-5, -7)': false,
  }));

  test('NotSmaller', () => check({
    'NotSmaller(5, 5)': true,
    'NotSmaller(6, 5)': true,
    'NotSmaller(5, 6)': false,
    'NotSmaller(-1, 0)': false,
  }));

  test('Pow', () => check({
    'Pow(2, 0)': 1,
    'Pow(7, 1)': 7,
    'Pow(2, 3)': 8,
    'Pow(2, -2)': 0.25,
  }));

  test('Qualifier', () => check({
    'Qualifier(Qnum(1.5, "="))': '=',
    'Qualifier(Qnum(1.5, "<"))': '<',
    'Qualifier(Qnum(1.5, ">"))': '>',
    'Qualifier(1)': '=',
  }));

  test('RandBetween', () => checkRandomInt({
    'RandBetween(5, 7)': [5, 7],
    'RandBetween(-2, 2)': [-2, 2],
    'RandBetween(0, 35)': [0, 35],
    'RandBetween(-100, -50)': [-100, -50],
    'RandBetween(1, 2)': [1, 2],
  }));

  test('Rnd', () => checkRandomInt({
    'Rnd(80)': [0, 80],
    'Rnd(2)': [0, 2],
    'Rnd(-2)': [0, 2],
  }));

  test('Round', () => check({
    'Round(3.4)': 3,
    'Round(3.5)': 4,
    'Round(1)': 1,
    'Round(-3.5)': -4,
    'Round(-3.4)': -3,
  }));

  test('RoundFloat', () => check({
    'RoundFloat(12345.12345, 6)': 12345.12345,
    'RoundFloat(12345.12345, 4)': 12345.1235,
    'RoundFloat(12345.12345, 0)': 12345,
    'RoundFloat(PI, 2)': 3.14,
    'RoundFloat(0.5, 0)': 1,
    'RoundFloat(0.3, 0)': 0,
    'RoundFloat(-0.5, 0)': -1,
    'RoundFloat(175, -1)': 180,
    'RoundFloat(170, -1)': 170,
    'RoundFloat(175, -2)': 200,
    'RoundFloat(125, -2)': 100,
    'RoundFloat(12340.12345, -3.8)': 12000,
    'RoundFloat(12340.12345, -4.2)': 10000,
    'RoundFloat(12340.12345, -5)': 0,
    'RoundFloat(null, 2)': undefined,
  }));

  test('Sin', () => check({
    'Sin(0)': 0,
    'Sin(PI / 6)': 0.5,
    'Sin(PI / 4)': 0.707,
    'Sin(PI / 3)': 0.866,
    'Sin(PI / 2)': 1,
  }));

  test('Smaller', () => check({
    'Smaller(4, 5)': true,
    'Smaller(5, 5)': false,
    'Smaller(5, 4)': false,
    'Smaller(-10, -5)': true,
  }));

  test('Sqrt', () => check({
    'Sqrt(1)': 1,
    'Sqrt(4)': 2,
    'Sqrt(9)': 3,
    'Sqrt(6.25)': 2.5,
  }));

  test('Sub', () => check({
    'Sub(10, 3)': 7,
    'Sub(124.161, 24.06)': 100.101,
    'Sub(-50, 25)': -75,
    'Sub(-11, -22)': 11,
    'Sub(999, 0)': 999,
    'Sub(0, -12)': 12,
  }));

  test('Tan', () => check({
    'Tan(0)': 0,
    'Tan(PI / 6)': 0.577,
    'Tan(PI / 4)': 1.0001, // enables float comparison
    'Tan(PI / 3)': 1.732,
  }));
});
