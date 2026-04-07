import {category, test} from '@datagrok-libraries/test/src/test';
import {check} from './utils';


category('Functions: Statistical', () => {
  test('Avg', () => check({
    'Avg([1, 2, 3, 4])': 2.5,
    'Avg([null, 0.3, 0.7])': 0.5,
    'Avg([5.123])': 5.123,
    'Avg([-1, -2, 3])': 0,
    'Avg([1, null, 2, 3])': 2,
    'Avg([])': undefined,
    'Avg([null, null])': undefined,
  }));

  test('Kurt', () => check({
    'Kurt([1, 2, 3])': -1.5,
    'Kurt([1, 2, 6])': -1.5,
    'Kurt([1, null, 2, 3])': -1.5,
    'Kurt([null, null])': 0,
  }));

  test('Max', () => check({
    'Max([1, 2, 4, 3])': 4,
    'Max([1, null, 2, 3])': 3,
    'Max([null, 1, 0.7, 0.3])': 1,
    'Max([1.5, -2, 1.9])': 1.9,
    'Max([null, null])': undefined,
  }));

  test('Med', () => check({
    'Med([1, 2, 3])': 2,
    'Med([1, 2, 3, 4])': 2.5,
    'Med([1, null, 2, 3])': 2,
    'Med([-100, null, 100])': 0,
    'Med([null, 0.3, 0.7])': 0.5,
    'Med([null, 1, 0.7, 0.3])': 0.7,
    'Med([null, null])': 0,
  }));

  test('Min', () => check({
    'Min([1.5, -2, 1.9])': -2,
    'Min([1, 2, 4, 3])': 1,
    'Min([2, null, 0, 3])': 0,
    'Min([null, 1, 0.7, 0.3])': 0.3,
    'Min([null, null])': undefined,
  }));

  test('MissingValueCount', () => check({
    'MissingValueCount([1, 2, 3])': 0,
    'MissingValueCount([10, null, 7])': 1,
    'MissingValueCount([2, null, 0, null, 3])': 2,
    'MissingValueCount([null, null, null])': 3,
  }));

  test('Percentile', () => check({
    'Percentile([10, 9, 8, 7, 6, 5, 4, 3, 2, 1], 0.25)': 3,
    'Percentile([1, 2, 3, 4], 0.25)': 2,
    'Percentile([1, 2, 3, 4], 0.40)': 2,
    'Percentile([1, 2, 3, 4], 0.75)': 4,
    'Percentile([1, 2, null, 3, 4], 0.25)': 2,
    'Percentile([null], 0.4)': undefined,
    'Percentile([], 0.4)': undefined,
  }));

  test('Q1', () => check({
    'Q1([-5, -3, -1, 0, 1, 3, 5])': -3,
    'Q1([5, -5, 1, -1, 3, -3, 0])': -3,
    'Q1([1, 2, 3])': 1,
    'Q1([1, null, 2, 3])': 1,
    'Q1([null, null])': 0,
  }));

  test('Q2', () => check({
    'Q2([1, 2, 3])': 2,
    'Q2([1, 2, 3, 4])': 2.5,
    'Q2([1, null, 2, 3])': 2,
    'Q2([-100, null, 100])': 0,
    'Q2([null, 0.3, 0.7])': 0.5,
    'Q2([null, 1, 0.7, 0.3])': 0.7,
    'Q2([null, null])': 0,
  }));

  test('Q3', () => check({
    'Q3([-5, -3, -1, 0, 1, 3, 5])': 3,
    'Q3([5, -5, 1, -1, 3, -3, 0])': 3,
    'Q3([1, 2, 3])': 3,
    'Q3([1, null, 2, 3])': 3,
    'Q3([null, null])': 0,
  }));

  test('Skew', () => check({
    'Skew([1, 2, 3])': 0,
    'Skew([1, 2, 6])': 0.595,
    'Skew([1, null, 2, 3])': 0,
    'Skew([null, null])': 0,
  }));

  test('StDev', () => check({
    'StDev([1, 2, 3])': 1,
    'StDev([1, null, 2, 3])': 1,
    'StDev([null, null])': 0,
    'StDev([7, 14, 21])': 7,
    'StDev([-15, -5, 5, 15])': 12.91,
  }));

  test('Sum', () => check({
    'Sum([1, 2, 4, 3])': 10,
    'Sum([-1, 4, 12, 5])': 20,
    'Sum([2, null, 0, 3])': 5,
    'Sum([null, 1, 0.7, 0.3])': 2,
    'Sum([null, null])': 0,
    'Sum([-0])': 0,
  }));

  test('TotalCount', () => check({
    'TotalCount([2, null, 0, 3])': 4,
    'TotalCount([1, 2, 4])': 3,
    'TotalCount([null, null])': 2,
    'TotalCount([100])': 1,
    'TotalCount([])': 0,
    'TotalCount(null)': undefined,
  }));

  test('ValueCount', () => check({
    'ValueCount([1, 2, 4, 3])': 4,
    'ValueCount([2, null, 0, 3])': 3,
    'ValueCount([1, 2, 4])': 3,
    'ValueCount([null, null])': 0,
    'ValueCount([100])': 1,
    'ValueCount([])': 0,
    'ValueCount(null)': undefined,
  }));

  test('Variance', () => check({
    'Variance([1, 2, 3])': 1,
    'Variance([1, null, 2, 3])': 1,
    'Variance([null, null])': 0,
    'Variance([7, 14, 21])': 49,
    'Variance([-15, -5, 5, 15])': 166.667,
  }));
});
