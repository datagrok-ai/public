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
};
