import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/utils/src/test';


category('Text functions', () => {
  const gfe = grok.functions.eval;

  test('Add', async () => {
    expect(await gfe('Add("bitter", "sweet")'), 'bittersweet');
    expect(await gfe('Add("", "right")'), 'right');
    expect(await gfe('Add("left", "")'), 'left');
    expect(await gfe('Add("", "")'), '');
  });

  test('Contains', async () => {
    expect(await gfe('Contains("stormy weather", "weather")'), true);
    expect(await gfe('Contains("stormy weather", "sunny")'), false);
    expect(await gfe('Contains("empty", "")'), true);
    expect(await gfe('Contains("Case", "case")'), false);
  });

  test('EndsWith', async () => {
    expect(await gfe('EndsWith("White Christmas", "Christmas")'), true);
    expect(await gfe('EndsWith("White Christmas", "White")'), false);
    expect(await gfe('EndsWith("White Christmas", "")'), true);
  });

  test('Eq', async () => {
    expect(await gfe('Eq("sky", "sky")'), true);
    expect(await gfe('Eq("", "")'), true);
    expect(await gfe('Eq("SKY", "sky")'), false);
    expect(await gfe('Eq(" sky ", "sky")'), false);
  });

  test('IsEmpty', async () => {
    expect(await gfe('IsEmpty("")'), true);
    expect(await gfe('IsEmpty(null)'), true);
    expect(await gfe('IsEmpty("dream")'), false);
    expect(await gfe('IsEmpty("     ")'), false);
  });

  test('IsNotEmpty', async () => {
    expect(await gfe('IsNotEmpty("")'), false);
    expect(await gfe('IsNotEmpty(null)'), false);
    expect(await gfe('IsNotEmpty("dream")'), true);
    expect(await gfe('IsNotEmpty("     ")'), true);
  });

  test('Length', async () => {
    expect(await gfe('Length("supercalifragilisticexpialidocious")'), 34);
    expect(await gfe('Length("bright dawn")'), 11);
    expect(await gfe('Length("")'), 0);
  });

  test('NotEq', async () => {
    expect(await gfe('NotEq("sky", "sky")'), false);
    expect(await gfe('NotEq("", "")'), false);
    expect(await gfe('NotEq("SKY", "sky")'), true);
    expect(await gfe('NotEq(" sky ", "sky")'), true);
  });

  test('RegExpContains', async () => {
    expect(await gfe('RegExpContains("stormy weather", "weather")'), true);
    expect(await gfe('RegExpContains("stormy weather", "sunny")'), false);
    expect(await gfe('RegExpContains("empty", "")'), true);
    expect(await gfe('RegExpContains("Case", "case")'), false);
    expect(await gfe('RegExpContains("name@gmail.com", "[\w.\-]{0,25}@(hotmail|gmail)\.com")'), true);
  });

  test('RegExpExtract', async () => {
    expect(await gfe('RegExpExtract("Hello, world!", "l+", 0)'), 'll');
    expect(await gfe('RegExpExtract("Hello, world!", "l+", 1)'), 'l');
  });

  test('RegExpReplace', async () => {
    expect(await gfe('RegExpReplace("Hello, world!", "l+", "LL")'), 'HeLLo, worLLd!');
  });

  test('ReplaceAll', async () => {
    expect(await gfe('ReplaceAll("New York", "York", "Orleans")'), 'New Orleans');
    expect(await gfe('ReplaceAll("every", "", ".")'), '.e.v.e.r.y.');
    expect(await gfe('ReplaceAll("day-to-day", "-", " ")'), 'day to day');
    expect(await gfe('ReplaceAll("spaceless", " ", "-")'), 'spaceless');
  });

  test('SplitString', async () => {
    expect(await gfe('SplitString("a,b,c,d", ",", 0)'), 'a');
    expect(await gfe('SplitString("devil-may-care", "-", 1)'), 'may');
  });

  test('StartsWith', async () => {
    expect(await gfe('StartsWith("Sunrise", "Sun")'), true);
    expect(await gfe('StartsWith("Sunrise", "sun")'), false);
    expect(await gfe('StartsWith("Sunrise", "moon")'), false);
    expect(await gfe('StartsWith("Sunrise", "")'), true);
  });

  test('StrFind', async () => {
    expect(await gfe('StrFind("Hello, world!", "Hello")'), 0);
    expect(await gfe('StrFind("Hello, world!", "hello")'), -1);
    expect(await gfe('StrFind("Hello, world!", "world")'), 7);
    expect(await gfe('StrFind("Hello, world!", "sun")'), -1);
    expect(await gfe('StrFind("Hello, world!", "Hello, world!")'), 0);
    expect(await gfe('StrFind("", "moon")'), -1);
    expect(await gfe('StrFind("moon", "")'), -1);
    expect(await gfe('StrFind("", "")'), -1);
  });

  test('StrLeft', async () => {
    expect(await gfe('StrLeft("crystal", 100)'), 'crystal');
    expect(await gfe('StrLeft("crystal", 3)'), 'cry');
    expect(await gfe('StrLeft("crystal", 0)'), '');
    expect(await gfe('StrLeft("", 3)'), '');
    expect(await gfe('StrLeft("crystal", -1)'), 'crysta');
    expect(await gfe('StrLeft("crystal", -7)'), '');
    expect(await gfe('StrLeft("crystal", -10)'), '');
  });

  test('StrRight', async () => {
    expect(await gfe('StrRight("crystal", 100)'), 'crystal');
    expect(await gfe('StrRight("crystal", 3)'), 'tal');
    expect(await gfe('StrRight("crystal", 0)'), '');
    expect(await gfe('StrRight("", 3)'), '');
    expect(await gfe('StrRight("crystal", -1)'), 'rystal');
    expect(await gfe('StrRight("crystal", -7)'), '');
    expect(await gfe('StrRight("crystal", -10)'), '');
  });

  test('StrRepeat', async () => {
    expect(await gfe('StrRepeat("a", 5)'), 'aaaaa');
    expect(await gfe('StrRepeat("ma", 2)'), 'mama');
    expect(await gfe('StrRepeat("", 2)'), '');
    expect(await gfe('StrRepeat("day", 1)'), 'day');
    expect(await gfe('StrRepeat("light", 0)'), '');
  });

  test('Substring', async () => {
    expect(await gfe('Substring("Snow storm", 5, 10)'), 'storm');
    expect(await gfe('Substring("stars", 0, 5)'), 'stars');
    expect(await gfe('Substring("galaxy", 0, 0)'), '');
    expect(await gfe('Substring("", 0, 0)'), '');
  });

  test('ToLowerCase', async () => {
    expect(await gfe('ToLowerCase("ICE")'), 'ice');
    expect(await gfe('ToLowerCase("snow")'), 'snow');
    expect(await gfe('ToLowerCase("Wind")'), 'wind');
  });

  test('ToUpperCase', async () => {
    expect(await gfe('ToUpperCase("home")'), 'HOME');
    expect(await gfe('ToUpperCase("CAT")'), 'CAT');
    expect(await gfe('ToUpperCase("Toy")'), 'TOY');
  });

  test('Trim', async () => {
    expect(await gfe('Trim("   Outer space   ")'), 'Outer space');
    expect(await gfe('Trim("spacecraft")'), 'spacecraft');
    expect(await gfe('Trim("   ")'), '');
    expect(await gfe('Trim("")'), '');
  });
});
