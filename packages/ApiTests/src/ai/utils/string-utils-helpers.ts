import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// JS API source: public/js-api/src/helpers.ts:185
// Two static helpers on DG.StringUtils: hashCode(s) — pure-JS 32-bit string
// hash with the (h<<5) - h + chr recurrence — and camelCaseToSentence(s, opts?)
// — Dart-backed converter that splits a camelCase identifier into a
// space-separated phrase, with toggles for capitalising the first word and
// conjunctions (and / or / than / if / but / so / as / that). Note:
// splitCamelCase uses ASCII A-Z bounds, so unicode behavior is undefined and
// not asserted.
category('AI: Utils: string-utils helpers', () => {
  test('hashCode empty string', async () => {
    expect(DG.StringUtils.hashCode(''), 0);
  });

  test('hashCode single char', async () => {
    expect(DG.StringUtils.hashCode('a'), 'a'.charCodeAt(0));
  });

  test('hashCode is deterministic and discriminates inputs', async () => {
    const samples = ['hello', 'world', 'datagrok', 'a', 'ab', 'ba'];
    for (const s of samples)
      expect(DG.StringUtils.hashCode(s), DG.StringUtils.hashCode(s));
    const set = new Set<number>();
    for (const s of samples)
      set.add(DG.StringUtils.hashCode(s));
    expect(set.size, samples.length);
  });

  test('hashCode returns 32-bit signed integer', async () => {
    const long = 'the quick brown fox jumps over the lazy dog 0123456789';
    const h = DG.StringUtils.hashCode(long);
    expect(Number.isInteger(h), true);
    expect(h >= -(2 ** 31), true);
    expect(h < 2 ** 31, true);
    expect(h | 0, h);
  });

  test('camelCaseToSentence default', async () => {
    expect(DG.StringUtils.camelCaseToSentence('helloWorld'), 'Hello World');
  });

  test('camelCaseToSentence capitalizeFirst=false', async () => {
    expect(DG.StringUtils.camelCaseToSentence('helloWorld', {capitalizeFirst: false}), 'hello World');
  });

  test('camelCaseToSentence capitalizeConjunctions', async () => {
    expect(DG.StringUtils.camelCaseToSentence('fooAndBar'), 'Foo and Bar');
    expect(DG.StringUtils.camelCaseToSentence('fooAndBar', {capitalizeConjunctions: true}), 'Foo And Bar');
  });

  test('camelCaseToSentence pass-through', async () => {
    expect(DG.StringUtils.camelCaseToSentence('ID'), 'ID');
    expect(DG.StringUtils.camelCaseToSentence('hello world'), 'hello world');
  });
});
