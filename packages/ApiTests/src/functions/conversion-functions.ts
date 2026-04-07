import {category, test} from '@datagrok-libraries/test/src/test';
import {check} from './utils';


category('Functions: Conversions', () => {
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

  test('ParseQnum', () => check({
    'ParseQnum("10")': 10.000000000000004,
    'ParseQnum("<10")': 10.000000000000002,
    'ParseQnum(">10")': 10.000000000000005,
    'ParseQnum(" < 10")': 10.000000000000002,
    'ParseQnum(">  10")': 10.000000000000005,
  }));

  test('Qnum', () => check({
    'ToString(Qnum(1.5, "="))': '1.5000000000000004',
    'ToString(Qnum(1.5, "<"))': '1.5000000000000002',
    'ToString(Qnum(1.5, ">"))': '1.5000000000000007',
    'ToString(Qnum(-1, "="))': '-1.0000000000000004',
    'ToString(Qnum(-1, "<"))': '-1.0000000000000002',
    'ToString(Qnum(-1, ">"))': '-1.0000000000000007',
  }));

  test('QnumToString', () => check({
    'QnumToString(Qnum(1.5, "="))': '1.50',
    'QnumToString(Qnum(1.5, "<"))': '<1.50',
    'QnumToString(Qnum(1.5, ">"))': '>1.50',
    'QnumToString(Qnum(1.115, "="))': '1.11',
    'QnumToString(Qnum(1.115, "<"))': '<1.11',
    'QnumToString(Qnum(1.115, ">"))': '>1.11',
  }));

  test('ToString', () => check({
    'ToString(1)': '1',
    'ToString(3.14)': '3.14',
    'ToString(true)': 'true',
  }));
});
