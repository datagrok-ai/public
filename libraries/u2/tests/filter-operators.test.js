import {test} from 'node:test';
import assert from 'node:assert/strict';
import {BitArray, SEMTYPE, TYPE} from 'datagrok-api/u2core';
import {Filters, FilterError, OperatorRegistry, CORE_OPERATORS, KIND} from '../src/core/filter/index.js';

const NOW = new Date('2026-09-04T12:00:00.000Z');
const ctx = {now: NOW};

function ids(prop) {
  return Filters.operators.for(prop).map((o) => o.id);
}

test('for(prop): the core operators of each kind', () => {
  assert.deepEqual(ids({name: 's', type: TYPE.STRING}),
    ['=', '!=', 'in', 'not in', 'like', '!like', 'starts', 'ends', 'matches', '!matches', 'fuzzy', 'is null', 'is not null']);
  assert.deepEqual(ids({name: 'i', type: TYPE.INT}),
    ['=', '!=', '>', '>=', '<', '<=', 'between', 'in', 'not in', 'is null', 'is not null']);
  assert.deepEqual(ids({name: 'f', type: TYPE.FLOAT}), ids({name: 'i', type: TYPE.INT}));
  assert.deepEqual(ids({name: 'b', type: TYPE.BIG_INT}), ids({name: 'i', type: TYPE.INT}));
  assert.deepEqual(ids({name: 'd', type: TYPE.DATE_TIME}),
    ['=', '!=', '>', '>=', '<', '<=', 'between', 'is null', 'is not null']);
  assert.deepEqual(ids({name: 'b', type: TYPE.BOOL}), ['=', '!=', 'is null', 'is not null']);
  assert.deepEqual(ids({name: 'l', type: TYPE.LIST}), ['like', '!like', 'is null', 'is not null']);
  assert.deepEqual(ids({name: 'r', type: TYPE.STRING, ref: 'Core.users'}), ['=', '!=', 'in', 'not in', 'is null', 'is not null']);
});

test('semType set wins over core for that semType, and only there', () => {
  const mol = {name: 'smiles', type: TYPE.STRING, semType: SEMTYPE.MOLECULE};
  const contains = {id: 'Contains', label: 'Contains', arity: 1, kinds: [KIND.STRING], semType: SEMTYPE.MOLECULE, editor: 'default',
    bitset: async () => ({bits: new Uint32Array(0), length: 0})};
  const off = Filters.operators.register(contains);
  try {
    assert.equal(ids(mol)[0], 'Contains');
    assert.equal(ids(mol).includes('like'), true, 'core operators follow the semType set');
    assert.equal(ids({name: 's', type: TYPE.STRING}).includes('Contains'), false);
    assert.equal(ids({name: 'n', type: TYPE.INT, semType: SEMTYPE.MOLECULE}).includes('Contains'), false, 'kinds still apply');
    assert.equal(Filters.operators.get('Contains', mol), contains);
    assert.equal(Filters.operators.get('Contains', {name: 's', type: TYPE.STRING}), undefined);
    assert.equal(Filters.operators.get('Contains'), contains, 'without a property any registration answers');
    assert.equal(Filters.operators.get('=').semType, undefined, 'without a property core comes first');
  } finally {
    off();
  }
  assert.equal(Filters.operators.get('Contains'), undefined);
  assert.equal(Filters.operators.all().length, CORE_OPERATORS.length);
});

test('register: the same id + semType replaces in place; unregister removes only what is still there', () => {
  const reg = new OperatorRegistry();
  reg.register(CORE_OPERATORS);
  const molEq1 = {id: '=', label: 'same as', arity: 1, kinds: [KIND.STRING], semType: SEMTYPE.MOLECULE, editor: 'default'};
  const molEq2 = {...molEq1, label: 'identical to'};
  const off1 = reg.register(molEq1);
  reg.register(molEq2);
  const mol = {name: 'smiles', type: TYPE.STRING, semType: SEMTYPE.MOLECULE};
  assert.equal(reg.all().length, CORE_OPERATORS.length + 1);
  assert.equal(reg.get('=', mol), molEq2);
  assert.equal(reg.get('=', {name: 's', type: TYPE.STRING}).label, 'equals', 'the core entry is untouched');
  off1();
  assert.equal(reg.get('=', mol), molEq2, 'the replaced entry was already gone; the replacement stays');
  const fuzzy2 = {...reg.get('fuzzy'), label: 'sounds like'};
  reg.register(fuzzy2);
  assert.equal(reg.all().length, CORE_OPERATORS.length + 1);
  assert.equal(reg.get('fuzzy').label, 'sounds like');
  assert.equal(reg.all().indexOf(fuzzy2), CORE_OPERATORS.findIndex((o) => o.id === 'fuzzy'), 'in place');
});

test('kindOf: every IProperty.type spelling, propertyType first, ref and an explicit kind win', () => {
  const kind = (type) => Filters.kindOf({name: 'p', type});
  assert.equal(kind(TYPE.INT), KIND.INT);
  assert.equal(kind(TYPE.BIG_INT), KIND.BIG_INT);
  assert.equal(kind(TYPE.FLOAT), KIND.FLOAT);
  assert.equal(kind(KIND.FLOAT), KIND.FLOAT);
  assert.equal(kind(TYPE.NUM), KIND.FLOAT);
  assert.equal(kind(TYPE.QNUM), KIND.FLOAT);
  assert.equal(kind(TYPE.BOOL), KIND.BOOL);
  assert.equal(kind(TYPE.DATE_TIME), KIND.DATE_TIME);
  assert.equal(kind(TYPE.STRING), KIND.STRING);
  assert.equal(kind(TYPE.LIST), KIND.STRING_LIST);
  assert.equal(kind(TYPE.STRING_LIST), KIND.STRING_LIST);
  assert.equal(kind(TYPE.OBJECT), KIND.STRING);
  assert.equal(kind(undefined), KIND.STRING);
  assert.equal(Filters.kindOf({name: 'p', type: TYPE.STRING, propertyType: TYPE.INT}), KIND.INT);
  assert.equal(Filters.kindOf({name: 'p', type: TYPE.STRING, ref: 'Core.users'}), KIND.REF);
  assert.equal(Filters.kindOf({name: 'p', type: TYPE.INT, kind: KIND.DATE_TIME}), KIND.DATE_TIME);
});

test('resolveSpan: signed spans against now, m = 31 days, y = 365; a bad span throws', () => {
  const at = (span) => Filters.resolveSpan(span, NOW).getTime() - NOW.getTime();
  assert.equal(at('now'), 0);
  assert.equal(at('-1h'), -3600e3);
  assert.equal(at('2d'), 2 * 86400e3);
  assert.equal(at('-1w'), -7 * 86400e3);
  assert.equal(at('1m'), 31 * 86400e3);
  assert.equal(at('-1y'), -365 * 86400e3);
  assert.throws(() => Filters.resolveSpan('1x', NOW), FilterError);
});

test('domain: values leave as the wire carries them', () => {
  const str = {name: 's', type: TYPE.STRING};
  const dom = (c, prop = str) => Filters.operators.get(c.operator, prop).domain(c, prop, ctx);
  const c = (operator, value, options) => Filters.cond('p', operator, value, options ? {options} : undefined);
  assert.deepEqual(dom(c('=', 'x')), {property: 'p', operator: '=', value: 'x'});
  assert.deepEqual(dom(c('!=', 5), {name: 'i', type: TYPE.INT}), {property: 'p', operator: '!=', value: 5});
  assert.deepEqual(dom(c('>', NOW), {name: 'd', type: TYPE.DATE_TIME}),
    {property: 'p', operator: '>', value: '2026-09-04T12:00:00.000Z'});
  assert.deepEqual(dom(c('>=', {span: '-1d'}), {name: 'd', type: TYPE.DATE_TIME}),
    {property: 'p', operator: '>=', value: '2026-09-03T12:00:00.000Z'}, 'a span resolves against ctx.now');
  assert.deepEqual(dom(c('=', {type: 'User', id: '42', name: 'Alice'}), {name: 'u', type: TYPE.STRING, ref: 'Core.users'}),
    {property: 'p', operator: '=', value: '42'});
  assert.deepEqual(dom(c('=', '@current'), {name: 'u', type: TYPE.STRING, ref: 'Core.users'}),
    {property: 'p', operator: '=', value: '@current'});
  assert.deepEqual(dom(c('in', ['a', 1]), {name: 'i', type: TYPE.INT}), {property: 'p', operator: '=', value: ['a', 1]});
  assert.deepEqual(dom(c('not in', [{type: 'User', id: '1'}]), {name: 'u', type: TYPE.STRING, ref: 'Core.users'}),
    {property: 'p', operator: '!=', value: ['1']});
  assert.deepEqual(dom(c('between', [1, {span: 'now'}]), {name: 'd', type: TYPE.DATE_TIME}),
    [{property: 'p', operator: '>=', value: 1}, 'and', {property: 'p', operator: '<=', value: '2026-09-04T12:00:00.000Z'}]);
  assert.deepEqual(dom(c('like', '50%_a\\b')), {property: 'p', operator: 'like', value: '%50\\%\\_a\\\\b%'});
  assert.deepEqual(dom(c('like', 'a%', {raw: true})), {property: 'p', operator: 'like', value: 'a%'});
  assert.deepEqual(dom(c('!like', 'x')), {property: 'p', operator: 'not like', value: '%x%'});
  assert.deepEqual(dom(c('starts', 'x')), {property: 'p', operator: 'like', value: 'x%'});
  assert.deepEqual(dom(c('ends', 'x')), {property: 'p', operator: 'like', value: '%x'});
  assert.deepEqual(dom(c('matches', '^a.*')), {property: 'p', operator: '~*', value: '^a.*'});
  assert.deepEqual(dom(c('!matches', '^a.*')), {property: 'p', operator: '!~*', value: '^a.*'});
  assert.deepEqual(dom(c('fuzzy', 'asp')), [{property: 'p', operator: 'fuzzy', threshold: null, value: 'asp'}, 'or',
    {property: 'p', operator: 'like', value: '%asp%'}]);
  assert.deepEqual(dom(c('fuzzy', 'asp', {threshold: 0.6}))[0].threshold, 0.6);
  assert.deepEqual(dom(c('is null')), {property: 'p', operator: '=', value: null});
  assert.deepEqual(dom(c('is not null')), {property: 'p', operator: '!=', value: null});
});

test('the core table: every entry has a domain form and, except fuzzy, a mask', () => {
  for (const op of CORE_OPERATORS) {
    assert.equal(typeof op.domain, 'function', op.id);
    assert.equal(typeof op.mask, op.id === 'fuzzy' ? 'undefined' : 'function', op.id);
    assert.equal(op.bitset, undefined, op.id);
    assert.equal(op.semType, undefined, op.id);
  }
  assert.deepEqual(CORE_OPERATORS.filter((o) => o.arity === 0).map((o) => o.id), ['is null', 'is not null']);
  assert.deepEqual(CORE_OPERATORS.filter((o) => o.editor === 'list').map((o) => o.id), ['in', 'not in']);
  assert.deepEqual(CORE_OPERATORS.filter((o) => o.editor === 'range').map((o) => o.id), ['between']);
});

test('an exclusive semType set leaves only the null tests of the core beside it', () => {
  const mol = {name: 'smiles', type: TYPE.STRING, semType: SEMTYPE.MOLECULE};
  const has = {id: 'Contains', label: 'has substructure', arity: 1, kinds: [KIND.STRING], semType: SEMTYPE.MOLECULE,
    exclusive: true, editor: 'default', bitset: async () => ({bits: new Uint32Array(0), length: 0})};
  const off = Filters.operators.register(has);
  try {
    assert.deepEqual(ids(mol), ['Contains', 'is null', 'is not null']);
    assert.equal(Filters.operators.get('like', mol), undefined, 'string contains is gone for the molecule');
    assert.equal(Filters.operators.get('is null', mol).label, 'is empty');
    assert.equal(ids({name: 's', type: TYPE.STRING}).includes('like'), true, 'plain strings keep theirs');
    assert.equal(Filters.operators.get('like'), CORE_OPERATORS.find((o) => o.id === 'like'));
  } finally {
    off();
  }
  assert.equal(ids(mol).includes('like'), true, 'unregistered: the core set is back');
});

test('registerSet: a provider descriptor becomes an exclusive semType set; bits copy into a trimmed BitArray once', async () => {
  const mol = {name: 'smiles', type: TYPE.STRING, semType: SEMTYPE.MOLECULE};
  const seen = [];
  const words = new Uint32Array([0b11111101]);
  const set = {
    semType: SEMTYPE.MOLECULE, exclusive: true,
    operators: [
      {id: 'Contains', label: 'has substructure', arity: 1, bitset: async (col, c, signal) => {
        seen.push([col.name, c.value, signal instanceof AbortSignal]);
        return {bits: words, length: 6};
      }},
      {id: 'Similar', label: 'is similar to', arity: 1, editor: 'default', kinds: [KIND.STRING],
        bitset: async () => ({bits: new Uint32Array([0b101]).buffer, length: 3})},
    ],
  };
  const off = Filters.operators.registerSet(set);
  try {
    assert.deepEqual(ids(mol), ['Contains', 'Similar', 'is null', 'is not null']);
    assert.equal(ids({name: 's', type: TYPE.STRING}).includes('Contains'), false);
    const contains = Filters.operators.get('Contains', mol);
    assert.equal(contains.semType, SEMTYPE.MOLECULE);
    assert.equal(contains.exclusive, true);
    assert.deepEqual(contains.kinds, [KIND.STRING]);
    assert.equal(contains.editor, 'default');
    assert.equal(contains.mask, undefined);
    const col = {name: 'smiles', type: TYPE.STRING, semType: SEMTYPE.MOLECULE, length: 6, getRawData: () => new Int32Array(6),
      categories: ['']};
    const mask = await contains.bitset(col, Filters.cond('smiles', 'Contains', 'c1ccccc1'), new AbortController().signal);
    assert.ok(mask instanceof BitArray);
    assert.deepEqual([...mask.getBuffer()], [0b111101], 'the tail beyond length is cleared');
    assert.equal(mask.length, 6);
    assert.notEqual(mask.getBuffer(), words, 'the provider keeps its own words');
    assert.deepEqual(seen, [['smiles', 'c1ccccc1', true]]);
    const similar = await Filters.operators.get('Similar', mol).bitset(col, Filters.cond('smiles', 'Similar', 'x'),
      new AbortController().signal);
    assert.deepEqual([...similar.getBuffer()], [0b101], 'an ArrayBuffer is read as words');
    const frame = {rowCount: 6, column: (name) => name === 'smiles' ? col : null};
    assert.deepEqual(Array.from((await Filters.toMask(frame,
      Filters.group('and', [Filters.cond('smiles', 'Contains', 'c1ccccc1')]))).getSelectedIndexes()), [0, 2, 3, 4, 5]);
  } finally {
    off();
  }
  assert.equal(Filters.operators.get('Contains'), undefined, 'the unregister removes the whole set');
  assert.equal(Filters.operators.get('Similar'), undefined);
  assert.equal(Filters.operators.all().length, CORE_OPERATORS.length);
});

test('registerSet: a provider answering fewer words than its length rejects with a FilterError on the leaf', async () => {
  const off = Filters.operators.registerSet({semType: SEMTYPE.MOLECULE, operators: [
    {id: 'Contains', label: 'has substructure', arity: 1, bitset: async () => ({bits: new Uint32Array(1), length: 40})}]});
  try {
    const c = Filters.cond('smiles', 'Contains', 'C');
    await assert.rejects(Filters.operators.get('Contains', {name: 'smiles', type: TYPE.STRING, semType: SEMTYPE.MOLECULE})
      .bitset({name: 'smiles'}, c, new AbortController().signal), (e) => e instanceof FilterError &&
      e.problems[0].nodeId === c.id && e.problems[0].code === 'evaluation' && /Molecule.*"Contains"/.test(e.message));
  } finally {
    off();
  }
});

test('registerSet: a bad descriptor is rejected whole, checkSet names every fault', () => {
  const before = Filters.operators.all().length;
  const bitset = async () => ({bits: new Uint32Array(0), length: 0});
  assert.throws(() => Filters.operators.registerSet(null), FilterError);
  assert.throws(() => Filters.operators.registerSet({semType: SEMTYPE.MOLECULE, operators: []}), FilterError);
  assert.throws(() => Filters.operators.registerSet({semType: SEMTYPE.MOLECULE,
    operators: [{id: 'ok', label: 'ok', arity: 1, bitset}, {id: '', label: 'x', arity: 3}]}), FilterError);
  assert.equal(Filters.operators.all().length, before, 'nothing of a rejected set is registered');
  assert.deepEqual(OperatorRegistry.checkSet({semType: SEMTYPE.MOLECULE, operators: [{id: 'ok', label: 'ok', arity: 1, bitset}]}),
    []);
  assert.deepEqual(OperatorRegistry.checkSet(null), ['not an object']);
  assert.deepEqual(OperatorRegistry.checkSet({semType: '', exclusive: 'yes'}),
    ['semType must be a non-empty string', 'exclusive must be a boolean', 'operators must be a non-empty array']);
  assert.deepEqual(OperatorRegistry.checkSet({semType: SEMTYPE.MOLECULE,
    operators: [null, {id: '', label: 1, arity: 3, kinds: ['molecule'], editor: 'sketcher'}]}), [
    'operators[0]: not an object',
    'operators[1]: id must be a non-empty string',
    'operators[1]: label must be a string',
    "operators[1]: arity must be 0, 1, 2 or 'n'",
    'operators[1]: kinds must list string, int, float, bigint, datetime, bool, ref, string_list',
    'operators[1]: editor must be default, range, list, none',
    'operators[1]: bitset must be a function',
  ]);
});
