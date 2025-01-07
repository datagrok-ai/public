import * as DG from 'datagrok-api/dg';
// import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/utils/src/test';

category('ValueMatcher', () => {
  test('numerical pattern equals', async () => {
    const matcher = DG.ValueMatcher.numerical('5');
    expect(matcher.operator, '=');
    expect(matcher.match(5).toString(), 'true');
    expect(matcher.validate(10).toString(), '"10" does not match "5"');
  });

  test('numerical pattern equals =', async () => {
    const matcher = DG.ValueMatcher.numerical('= 5');
    expect(matcher.operator, '=');
    expect(matcher.match(5).toString(), 'true');
    expect(matcher.validate(10).toString(), '"10" does not match "= 5"');
  });

  test('numerical pattern not equals', async () => {
    const matcher = DG.ValueMatcher.numerical('!= 5');
    expect(matcher.operator, '!=');
    expect(matcher.match(5).toString(), 'false');
    expect(matcher.validate(5).toString(), '"5" does not match "!= 5"');
  });

  test('numerical pattern greater than', async () => {
    const matcher = DG.ValueMatcher.numerical('> 30');
    expect(matcher.operator, '>');
    expect(matcher.match(40).toString(), 'true');
    expect(matcher.validate(20).toString(), '"20" does not match "> 30"');
  });

  test('numerical pattern greater than or equals', async () => {
    const matcher = DG.ValueMatcher.numerical('>= 5');
    expect(matcher.operator, '>=');
    expect(matcher.match(10).toString(), 'true');
    expect(matcher.validate(1).toString(), '"1" does not match ">= 5"');
  });

  test('numerical pattern less than', async () => {
    const matcher = DG.ValueMatcher.numerical('< 10');
    expect(matcher.operator, '<');
    expect(matcher.match(15).toString(), 'false');
    expect(matcher.validate(20).toString(), '"20" does not match "< 10"');
  });

  test('numerical pattern less than or equals', async () => {
    const matcher = DG.ValueMatcher.numerical('<= 10');
    expect(matcher.operator, '<=');
    expect(matcher.match(5).toString(), 'true');
    expect(matcher.validate(15).toString(), '"15" does not match "<= 10"');
  });

  test('numerical pattern range(inclusive)', async () => {
    const matcher = DG.ValueMatcher.numerical('10-20');
    expect(matcher.operator, '-');
    expect(matcher.match(25).toString(), 'false');
    expect(matcher.validate(30).toString(), '"30" does not match "10-20"');
  });

  test('numerical pattern in', async () => {
    const matcher = DG.ValueMatcher.numerical('in(5,10)');
    expect(matcher.operator, 'in');
    expect(matcher.match(5).toString(), 'true');
    expect(matcher.validate(11).toString(), '"11" does not match "in(5,10)"');
  });

  test('numerical pattern not in', async () => {
    const matcher = DG.ValueMatcher.numerical('not in(5,10)');
    expect(matcher.operator, 'not in');
    expect(matcher.match(12).toString(), 'true');
    expect(matcher.validate(5).toString(), '"5" does not match "not in(5,10)"');
  });

  test('string pattern contains', async () => {
    const matcher = DG.ValueMatcher.string('contains foo');
    expect(matcher.operator, 'contains');
    expect(matcher.match('foofd').toString(), 'true');
    expect(matcher.validate('ffff').toString(), '"ffff" does not match "contains foo"');
  });

  test('string pattern starts with', async () => {
    const matcher = DG.ValueMatcher.string('starts with Charles');
    expect(matcher.operator, 'starts with');
    expect(matcher.match('hhhCharles').toString(), 'false');
    expect(matcher.validate('fggdgdgCharles').toString(), '"fggdgdgCharles" does not match "starts with Charles"');
  });

  test('string pattern ends with', async () => {
    const matcher = DG.ValueMatcher.string('ends with District');
    expect(matcher.operator, 'ends with');
    expect(matcher.match('helloDistrict').toString(), 'true');
    expect(matcher.validate('ffff').toString(), '"ffff" does not match "ends with District"');
  });

  test('string pattern in', async () => {
    const matcher = DG.ValueMatcher.string('in (asian, other)');
    expect(matcher.operator, 'in');
    expect(matcher.match('asian').toString(), 'true');
    expect(matcher.validate('nnnn').toString(), '"nnnn" does not match "in (asian, other)"');
  });

  test('string pattern not in', async () => {
    const matcher = DG.ValueMatcher.string('not in (m, f)');
    expect(matcher.operator, 'not in');
    expect(matcher.match('m').toString(), 'false');
    expect(matcher.validate('f').toString(), '"f" does not match "not in (m, f)"');
  });

  test('datetime pattern after the specified date', async () => {
    const matcher = DG.ValueMatcher.dateTime('after 10/17/2019');
    expect(matcher.operator, 'after');
    expect(matcher.match('12/03/2022').toString(), 'true');
    expect(matcher.validate('12/03/2022') ==undefined);
  });
}, { owner: 'oserhiienko@datagrok.ai' });
