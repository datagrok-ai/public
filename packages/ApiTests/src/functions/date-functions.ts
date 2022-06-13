import {category, test} from '@datagrok-libraries/utils/src/test';
import {check} from './utils';


category('Date functions', () => {
  test('DateDiff', () => check({
    'DateDiff(Date(2020, 1, 2), Date(2020, 1, 1))': 86400000,
    'DateDiff(Date(2020, 1, 1), Date(2020, 1, 2))': -86400000,
    'DateDiff(Date(2020, 1, 1), Date(2020, 1, 1))': 0,
    'DateDiff(DateTime(2020, 1, 1, 0, 0, 0, 125), Date(2020, 1, 1))': 125,
  }));

  test('DayOfMonth', () => check({
    'DayOfMonth(Date(2020, 6, 15))': 15,
    'DayOfMonth(Date(2020, 6, 1))': 1,
    'DayOfMonth(Date(2020, 6, 30))': 30,
  }));

  test('DayOfWeek', () => check({
    'DayOfWeek(Date(2020, 12, 31))': 4,
    'DayOfWeek(Date(2020, 1, 1))': 3,
    'DayOfWeek(Date(2020, 1, 5))': 7,
  }));

  test('DayOfYear', () => check({
    'DayOfYear(Date(2020, 1, 1))': 1,
    'DayOfYear(Date(2020, 2, 25))': 56,
    'DayOfYear(Date(2020, 12, 31))': 366,
  }));

  test('Hour', () => check({
    'Hour(DateTime(2020, 1, 1, 23, 59, 45, 999))': 23,
    'Hour(DateTime(2020, 1, 1, 0, 59, 45, 999))': 0,
    'Hour(DateTime(2020, 1, 1, 12, 0, 0, 0))': 12,
  }));

  test('Millisecond', () => check({
    'Millisecond(DateTime(2020, 1, 1, 0, 0, 0, 0))': 0,
    'Millisecond(DateTime(2020, 1, 1, 0, 0, 0, 100))': 100,
    'Millisecond(DateTime(2020, 1, 1, 0, 0, 0, 909))': 909,
  }));

  test('Minute', () => check({
    'Minute(DateTime(2020, 1, 1, 0, 0, 0, 0))': 0,
    'Minute(DateTime(2020, 1, 1, 0, 5, 50, 50))': 5,
    'Minute(DateTime(2020, 1, 1, 0, 59, 0, 0))': 59,
  }));

  test('Month', () => check({
    'Month(Date(2020, 1, 1))': 1,
    'Month(Date(2020, 3, 15))': 3,
    'Month(Date(2020, 12, 15))': 12,
  }));

  test('Quarter', () => check({
    'Quarter(Date(2020, 1, 1))': 1,
    'Quarter(Date(2020, 2, 10))': 1,
    'Quarter(Date(2020, 3, 15))': 1,
    'Quarter(Date(2020, 4, 1))': 2,
    'Quarter(Date(2020, 5, 1))': 2,
    'Quarter(Date(2020, 6, 1))': 2,
    'Quarter(Date(2020, 7, 1))': 3,
    'Quarter(Date(2020, 8, 1))': 3,
    'Quarter(Date(2020, 9, 1))': 3,
    'Quarter(Date(2020, 10, 1))': 4,
    'Quarter(Date(2020, 11, 1))': 4,
    'Quarter(Date(2020, 12, 1))': 4,
  }));

  test('Second', () => check({
    'Second(DateTime(2020, 1, 1, 23, 59, 59, 999))': 59,
    'Second(DateTime(2020, 1, 1, 23, 59, 45, 999))': 45,
    'Second(DateTime(2020, 1, 1, 23, 59, 0, 999))': 0,
  }));

  test('Weeknum', () => check({
    'Weeknum(Date(2020, 1, 1))': 0,
    'Weeknum(Date(2020, 2, 3))': 5,
    'Weeknum(Date(2020, 12, 31))': 53,
  }));

  test('Year', () => check({
    'Year(Date(2020, 1, 1))': 2020,
    'Year(Date(1999, 1, 1))': 1999,
    'Year(Date(2100, 1, 1))': 2100,
  }));
});
