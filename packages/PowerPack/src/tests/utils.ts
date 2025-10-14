import dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc';
dayjs.extend(utc);

export const FUNC_TESTS: {[f: string]: {[test: string]: any}} = {
  Boolean: {
    'Boolean(DateParse("20200131T132700"))': true,
  },
  ParseQnum: {
    'ParseQnum("100")': 100.00000000000003,
  },
  DateAdd: {
    'DateAdd(Date(2020, 1, 1), 3605000)': dayjs.utc('2020-01-01').add(1, 'h').add(5, 's'),
  },
  Avg: {
    'Avg($[age])': 55.599998474121094,
  },
  Abs: {
    'Abs(${height})': 162,
  },
};

export const FUNC_VALIDATION: { [f: string]: string } = {
  'Abs(num)': `Variable "num" not found`,
  'Abs()': `Abs: 1 required input parameters, 0 passed`,
  'Abs(\'i\')': `FormatException: Invalid double
i`,
  'DateParse(1)': `Function DateParse 's' param should be string type instead of number`,
  'DateParse(Abs(1))': `Function DateParse 's' param should be string type instead of num`,
  'Abs(': `Possible syntax error`,
  'Abs(${age1})': `Column age1 is missing`,
  'BinBySpecificLimits(${Age}, [18, 30, 45, 60, 75])': '',
};


export const FUNC_HINTS: { [f: string]: string } = {
  'DateParse(string)': `DateParse(s:string): datetime`,
  'Add(Abs(num), dynamic)': `Abs(x:num): num`,
};
