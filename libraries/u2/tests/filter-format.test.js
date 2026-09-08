/* `Filters.format` (schema-filters WO-2, D11): the canonical string — quoting and escapes,
   bracketed names, every operator spelling, lists, spans, refs, dates, nested groups, `not`,
   single-node inlining, empties — and its fixed-point property through the grammar. Imports
   through the package barrel on purpose: the only filter suite that exercises `src/index.ts`. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import './dom-shim.js';
import {Filters} from '../src/index.js';

const g = (op, nodes, extra) => Filters.group(op, nodes, extra);
const k = (p, o, v, extra) => Filters.cond(p, o, v, extra);

function reparses(text) {
  const {tree, errors} = Filters.parseTree(text);
  assert.deepEqual(errors, [], `${text} must parse`);
  assert.equal(Filters.format(Filters.fromDomainTree(tree)), text, 'a fixed point');
}

test('format: values — strings with every escape, numbers, booleans, null, @current, lists', () => {
  assert.equal(Filters.format(k('s', '=', 'a"b\'c\\d\ne\rf\tg\bh\fij')),
    's = "a\\"b\\\'c\\\\d\\ne\\rf\\tg\\bh\\fi\\u0001j"');
  reparses('s = "a\\"b\\\'c\\\\d\\ne\\rf\\tg\\bh\\fi\\u0001j"');
  assert.equal(Filters.format(k('n', '=', 1.5)), 'n = 1.5');
  assert.equal(Filters.format(k('n', '=', -0)), 'n = 0');
  assert.equal(Filters.format(k('b', '=', true)), 'b = true');
  assert.equal(Filters.format(k('b', '!=', false)), 'b != false');
  assert.equal(Filters.format(k('x', '=', null)), 'x = null');
  assert.equal(Filters.format(k('o', '=', '@current')), 'o = @current');
  assert.equal(Filters.format(k('s', 'in', ['a', 1, true, null])), 's in ("a", 1, true, null)');
  assert.equal(Filters.format(k('s', 'not in', ['a'])), 's not in ("a")');
  assert.equal(Filters.format(k('s', '=')), 's =', 'a missing value prints nothing');
});

test('format: dates, spans, refs', () => {
  assert.equal(Filters.format(k('d', '>', new Date('2026-01-02T03:04:05.000Z'))), 'd > "2026-01-02T03:04:05.000Z"');
  reparses('d > "2026-01-02T03:04:05.000Z"');
  assert.equal(Filters.format(k('d', '>', {span: '-1w'})), 'd > -1w');
  assert.equal(Filters.format(k('d', '<=', {span: 'now'})), 'd <= now');
  assert.equal(Filters.format(k('d', '>', Filters.parseTree('d > -2d').tree[0].value)), 'd > -2d', 'a parsed span');
  assert.equal(Filters.format(k('d', 'between', [{span: '-1w'}, {span: 'now'}])), 'd between -1w and now');
  assert.equal(Filters.format(k('u', '=', {type: 'Core.users', id: 'u1', name: 'Alice'})), 'u = "u1"');
  assert.equal(Filters.format(k('u', 'in', [{type: 'Core.users', id: 'u1'}, {type: 'Core.users', id: 'u2'}])),
    'u in ("u1", "u2")');
});

test('format: property names — bare paths, bracketed otherwise, with \\] and \\\\ escapes', () => {
  assert.equal(Filters.format(k('a.b_1.c2', '=', 1)), 'a.b_1.c2 = 1');
  assert.equal(Filters.format(k('Mol Weight', '=', 1)), '[Mol Weight] = 1');
  assert.equal(Filters.format(k('a]b\\c', '=', 1)), '[a\\]b\\\\c] = 1');
  assert.equal(Filters.format(k('1st', '=', 1)), '[1st] = 1');
  assert.equal(Filters.format(k('a.', '=', 1)), '[a.] = 1');
  assert.equal(Filters.format(k('a', '=', 1)), '[a\\u0001] = 1');
  reparses('[a\\]b\\\\c] = 1');
  reparses('[a\\u0001] = 1');
  assert.equal(Filters.fromDomainTree(Filters.parseTree('[a\\]b\\\\c] = 1').tree).nodes[0].property, 'a]b\\c');
});

test('format: every operator spelling', () => {
  const cases = [
    [k('a', '=', 1), 'a = 1'], [k('a', '!=', 1), 'a != 1'], [k('a', '>', 1), 'a > 1'], [k('a', '>=', 1), 'a >= 1'],
    [k('a', '<', 1), 'a < 1'], [k('a', '<=', 1), 'a <= 1'], [k('a', 'like', 'x'), 'a like "x"'],
    [k('a', '!like', 'x'), 'a !like "x"'], [k('a', 'starts', 'x'), 'a starts "x"'], [k('a', 'ends', 'x'), 'a ends "x"'],
    [k('a', 'matches', 'x'), 'a matches "x"'], [k('a', '!matches', 'x'), 'a !matches "x"'],
    [k('a', 'fuzzy', 'x'), 'a fuzzy "x"'], [k('a', 'fuzzy', 'x', {options: {threshold: 0.6}}), 'a fuzzy(0.6) "x"'],
    [k('a', 'in', [1, 2]), 'a in (1, 2)'], [k('a', 'not in', [1, 2]), 'a not in (1, 2)'],
    [k('a', 'between', [1, 2]), 'a between 1 and 2'], [k('a', 'is null'), 'a = null'], [k('a', 'is not null'), 'a != null'],
    [k('a', 'like', '%a%b%', {options: {raw: true}}), 'a like "a%b"'],
    [k('a', 'like', 'a%b', {options: {raw: true}}), 'a like "a%b"'],
    [k('a', '!like', '%a_b%', {options: {raw: true}}), 'a !like "a_b"'],
  ];
  for (const [cond, expected] of cases) {
    assert.equal(Filters.format(cond), expected);
    reparses(expected);
  }
});

test('format: a literal % or _ in a like shape is escaped, so the text reads back as the literal', () => {
  const cases = [
    [k('o', 'like', '50%'), 'o like "50\\\\%"'], [k('o', 'starts', 'a_b'), 'o starts "a\\\\_b"'],
    [k('o', 'ends', '50%'), 'o ends "50\\\\%"'], [k('o', '!like', 'a\\b'), 'o !like "a\\\\\\\\b"'],
    [k('o', 'like', '50%_\\'), 'o like "50\\\\%\\\\_\\\\\\\\"'],
  ];
  for (const [cond, expected] of cases) {
    const text = Filters.format(cond);
    assert.equal(text, expected);
    reparses(text);
    const back = Filters.fromDomainTree(Filters.parseTree(text).tree).nodes[0];
    assert.equal(back.value, cond.value, `${text} keeps the literal`);
    assert.equal(back.options, undefined, `${text} is not a raw pattern`);
  }
  assert.equal(Filters.format(k('o', 'like', '%50%%', {options: {raw: true}})), 'o like "50%"', 'raw stays raw');
});

test('format: groups — nesting parenthesized, single-node groups inlined, empties dropped, not', () => {
  const a = k('a', '=', 1);
  const b = k('b', '=', 2);
  const c = k('c', '=', 3);
  assert.equal(Filters.format(g('and', [])), '');
  assert.equal(Filters.format(g('or', [a])), 'a = 1');
  assert.equal(Filters.format(g('and', [a, b])), 'a = 1 and b = 2');
  assert.equal(Filters.format(g('or', [a, b])), 'a = 1 or b = 2');
  assert.equal(Filters.format(g('and', [a, g('or', [b, c])])), 'a = 1 and (b = 2 or c = 3)');
  assert.equal(Filters.format(g('or', [g('and', [a, b]), c])), '(a = 1 and b = 2) or c = 3');
  assert.equal(Filters.format(g('and', [a, g('or', [b])])), 'a = 1 and b = 2');
  assert.equal(Filters.format(g('and', [a, g('or', [])])), 'a = 1');
  assert.equal(Filters.format(g('and', [g('or', [])])), '');
  assert.equal(Filters.format(g('and', [a, g('or', [g('and', [b, c])])])), 'a = 1 and (b = 2 and c = 3)');
  assert.equal(Filters.format(g('and', [a, b], {not: true})), 'not (a = 1 and b = 2)');
  assert.equal(Filters.format(g('and', [a], {not: true})), 'not (a = 1)');
  assert.equal(Filters.format(g('and', [a, g('or', [b, c], {not: true})])), 'a = 1 and not (b = 2 or c = 3)');
  assert.equal(Filters.format(g('and', [g('or', [b, c])], {not: true})), 'not (b = 2 or c = 3)');
  assert.equal(Filters.format(g('and', [g('or', [b, c], {not: true})], {not: true})), 'not (not (b = 2 or c = 3))');
  assert.equal(Filters.format(g('and', [], {not: true})), '');
});

test('format: what the grammar reads back is the pushed-down tree; formatting that is the fixed point', () => {
  const root = g('and', [k('a', '=', 1), g('or', [k('b', '=', 2), k('c', 'like', 'x')], {not: true})]);
  const text = Filters.format(root);
  assert.equal(text, 'a = 1 and not (b = 2 or c like "x")');
  const back = Filters.fromDomainTree(Filters.parseTree(text).tree);
  assert.equal(Filters.format(back), 'a = 1 and (b != 2 and c !like "x")', 'the pushed-down group stays a group');
  reparses('a = 1 and (b != 2 and c !like "x")');
  assert.deepEqual(Filters.toDomainTree(root), Filters.toDomainTree(back));
});

test('formatValue and formatProperty are the spellings format uses', () => {
  assert.equal(Filters.formatValue('a "b"'), '"a \\"b\\""');
  assert.equal(Filters.formatValue(['x', 2, null, true]), '("x", 2, null, true)');
  assert.equal(Filters.formatValue({span: '-1w'}), '-1w');
  assert.equal(Filters.formatValue({type: 'Core.users', id: 'u1'}), '"u1"');
  assert.equal(Filters.formatValue('@current'), '@current');
  assert.equal(Filters.formatValue(undefined), '');
  assert.equal(Filters.formatProperty('mol.weight'), 'mol.weight');
  assert.equal(Filters.formatProperty('mol weight]'), '[mol weight\\]]');
  assert.equal(Filters.format(k('mol weight]', '=', 'a "b"')), '[mol weight\\]] = "a \\"b\\""');
});
