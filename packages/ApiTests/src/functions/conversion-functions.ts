import {category, test} from '@datagrok-libraries/utils/src/test';
import {check} from './utils';


category('Conversion functions', () => {
  test('Boolean', () => check({
    'Boolean(true)': true,
    'Boolean("true")': true,
    'Boolean("y")': true,
    'Boolean(1)': true,
    'Boolean(10)': true,
    'Boolean(DateParse("20200131T132700"))': true,
    'Boolean(false)': false,
    'Boolean("false")': false,
    'Boolean("n")': false,
    'Boolean("abc")': false,
    'Boolean("")': false,
    'Boolean(null)': false,
    'Boolean(0)': false,
  }));

  test('ParseFloat', () => check({
    'ParseFloat("2025")': 2025,
    'ParseFloat("12.78")': 12.78,
    'ParseFloat("-012.150")': -12.15,
  }));

  test('ParseInt', () => check({
    'ParseInt("2025")': 2025,
    'ParseInt("-012")': -12,
    'ParseInt(" 0101 ")': 101,
  }));


  test('ToString', () => check({
    'ToString(1)': '1',
    'ToString(3.14)': '3.14',
    'ToString(true)': 'true',
  }));
});
