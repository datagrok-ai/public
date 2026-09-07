import {category, test} from '@datagrok-libraries/test/src/test';
import {check} from './utils';


category('Functions: Text', () => {
  test('Add', () => check({
    'Add("bitter", "sweet")': 'bittersweet',
    'Add("", "right")': 'right',
    'Add("left", "")': 'left',
    'Add("", "")': '',
  }));

  test('Contains', () => check({
    'Contains("stormy weather", "weather")': true,
    'Contains("stormy weather", "sunny")': false,
    'Contains("empty", "")': true,
    'Contains("Case", "case")': false,
  }));

  test('EndsWith', () => check({
    'EndsWith("White Christmas", "Christmas")': true,
    'EndsWith("White Christmas", "White")': false,
    'EndsWith("White Christmas", "")': true,
  }));

  test('Eq', () => check({
    'Eq("sky", "sky")': true,
    'Eq("", "")': true,
    'Eq("SKY", "sky")': false,
    'Eq(" sky ", "sky")': false,
    'Eq("1", 1)': false,
    'Eq("", null)': false,
    'Eq(null, null)': true,
  }));

  test('IsEmpty', () => check({
    'IsEmpty("")': true,
    'IsEmpty(null)': true,
    'IsEmpty("dream")': false,
    'IsEmpty("     ")': false,
  }));

  test('IsNotEmpty', () => check({
    'IsNotEmpty("")': false,
    'IsNotEmpty(null)': false,
    'IsNotEmpty("dream")': true,
    'IsNotEmpty("     ")': true,
  }));

  test('Length', () => check({
    'Length("supercalifragilisticexpialidocious")': 34,
    'Length("bright dawn")': 11,
    'Length("")': 0,
  }));

  test('NotEq', () => check({
    'NotEq("sky", "sky")': false,
    'NotEq("", "")': false,
    'NotEq("SKY", "sky")': true,
    'NotEq(" sky ", "sky")': true,
  }));

  test('RegExpContains', () => check({
    'RegExpContains("stormy weather", "weather")': true,
    'RegExpContains("stormy weather", "sunny")': false,
    'RegExpContains("empty", "")': true,
    'RegExpContains("Case", "case")': false,
    'RegExpContains("name@gmail.com", "[\w.\-]{0,25}@(hotmail|gmail)\.com")': true,
  }));

  test('RegExpExtract', () => check({
    'RegExpExtract("Hello, world!", "l+", 0)': 'll',
    'RegExpExtract("Hello, world!", "l+", 1)': 'l',
  }));

  test('RegExpReplace', () => check({
    'RegExpReplace("Hello, world!", "l+", "LL")': 'HeLLo, worLLd!',
  }));

  test('ReplaceAll', () => check({
    'ReplaceAll("New York", "York", "Orleans")': 'New Orleans',
    'ReplaceAll("day-to-day", "-", " ")': 'day to day',
    'ReplaceAll("spaceless", " ", "-")': 'spaceless',
  }));

  test('SplitString', () => check({
    'SplitString("a,b,c,d", ",", 0)': 'a',
    'SplitString("devil-may-care", "-", 1)': 'may',
  }));

  test('StartsWith', () => check({
    'StartsWith("Sunrise", "Sun")': true,
    'StartsWith("Sunrise", "sun")': false,
    'StartsWith("Sunrise", "moon")': false,
    'StartsWith("Sunrise", "")': true,
  }));

  test('StrFind', () => check({
    'StrFind("Hello, world!", "Hello")': 0,
    'StrFind("Hello, world!", "hello")': -1,
    'StrFind("Hello, world!", "world")': 7,
    'StrFind("Hello, world!", "sun")': -1,
    'StrFind("Hello, world!", "Hello, world!")': 0,
    'StrFind("", "moon")': -1,
    'StrFind("moon", "")': -1,
    'StrFind("", "")': -1,
  }));

  test('StrLeft', () => check({
    'StrLeft("crystal", 100)': 'crystal',
    'StrLeft("crystal", 3)': 'cry',
    'StrLeft("crystal", 0)': '',
    'StrLeft("", 3)': '',
    'StrLeft("crystal", -1)': 'crysta',
    'StrLeft("crystal", -7)': '',
    'StrLeft("crystal", -10)': '',
  }));

  test('StrRight', () => check({
    'StrRight("crystal", 100)': 'crystal',
    'StrRight("crystal", 3)': 'tal',
    'StrRight("crystal", 0)': '',
    'StrRight("", 3)': '',
    'StrRight("crystal", -1)': 'rystal',
    'StrRight("crystal", -7)': '',
    'StrRight("crystal", -10)': '',
  }));

  test('StrRepeat', () => check({
    'StrRepeat("a", 5)': 'aaaaa',
    'StrRepeat("ma", 2)': 'mama',
    'StrRepeat("", 2)': '',
    'StrRepeat("day", 1)': 'day',
    'StrRepeat("light", 0)': '',
  }));

  test('Substring', () => check({
    'Substring("Snow storm", 5, 10)': 'storm',
    'Substring("stars", 0, 5)': 'stars',
    'Substring("galaxy", 0, 0)': '',
    'Substring("", 0, 0)': '',
  }));

  test('ToLowerCase', () => check({
    'ToLowerCase("ICE")': 'ice',
    'ToLowerCase("snow")': 'snow',
    'ToLowerCase("Wind")': 'wind',
  }));

  test('ToUpperCase', () => check({
    'ToUpperCase("home")': 'HOME',
    'ToUpperCase("CAT")': 'CAT',
    'ToUpperCase("Toy")': 'TOY',
  }));

  test('Trim', () => check({
    'Trim("   Outer space   ")': 'Outer space',
    'Trim("spacecraft")': 'spacecraft',
    'Trim("   ")': '',
    'Trim("")': '',
  }));
});
