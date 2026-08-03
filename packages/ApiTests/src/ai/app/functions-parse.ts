import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectNoThrow} from '../helpers';

category('AI: App: Functions Parse', () => {
  test('parse builds a FuncCall from a core expression', async () => {
    const fc = grok.functions.parse('Sin(1)');
    expect(fc instanceof DG.FuncCall, true);
    expect(fc.func != null, true);
  });

  test('parse with safe=false does not throw on a valid expression', async () => {
    expectNoThrow(() => grok.functions.parse('Sin(1)', false));
  });

  test('handleOuterBracketsInColName escape/unescape round-trips', async () => {
    // Per the Dart impl, escape:true turns a nested-column reference `${age}` into
    // `$\{age\}` (escapes `${`->`$\{` and `}`->`\}`); escape:false reverses it.
    const fns = grok.functions;
    const nested = '${age}';
    const escaped = fns.handleOuterBracketsInColName(nested, true);
    expect(escaped !== nested, true);
    const unescaped = fns.handleOuterBracketsInColName(escaped, false);
    expect(unescaped, nested);
    // A plain name has no nested brackets, so escape:true leaves it unchanged.
    expect(fns.handleOuterBracketsInColName('age', true), 'age');
  });
}, {owner: 'agolovko@datagrok.ai'});
