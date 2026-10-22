/* Host operators: answered by the host instead of the target (a structure search over a domain table), so they
   validate for 'domain' without a domain form, parse and format by their id, carry an optional numeric
   parameter in `options`, and a builder row edits that parameter after the value. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, Filters, FilterBuilder} from '../src/index.js';
import {TYPE} from 'datagrok-api/u2core';

const SCHEMA = Filters.schema([
  {name: 'structure', type: TYPE.STRING, semType: 'Molecule'},
  {name: 'name', type: TYPE.STRING},
]);

function hosted(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const off = Filters.operators.register([
      {id: 'has substructure', label: 'has substructure', arity: 1, kinds: [Filters.KIND.STRING], semType: 'Molecule',
        exclusive: true, editor: 'default', host: true},
      {id: 'is similar to', label: 'is similar to', arity: 1, kinds: [Filters.KIND.STRING], semType: 'Molecule',
        exclusive: true, editor: 'default', host: true,
        parameter: {name: 'threshold', label: 'Similarity at least', default: 0.7, min: 0, max: 1}},
    ]);
    try {
      Filters.resetIds('');
      await body();
    } finally {
      off();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

hosted('a host operator is offered on its semType alone and validates for a domain target', () => {
  const structure = Filters.property(SCHEMA, 'structure');
  assert.deepEqual(Filters.operators.for(structure).map((o) => o.id),
    ['has substructure', 'is similar to', 'is null', 'is not null']);
  const root = Filters.group('and', [Filters.cond('structure', 'has substructure', 'c1ccccc1')]);
  assert.deepEqual(Filters.validate(root, SCHEMA, 'domain'), []);
  assert.throws(() => Filters.toDomainTree(root), 'the host answers it; the server never sees it');
});

hosted('text names a host operator by its id, its parameter in parentheses, and formats back', () => {
  for (const text of ['structure has substructure "c1ccccc1"', 'structure is similar to(0.8) "CCO" and name = "x"']) {
    const {root, problems} = Filters.parse(text, SCHEMA, 'domain');
    assert.deepEqual(problems, []);
    assert.equal(Filters.format(root), text);
  }
  const {root} = Filters.parse('structure is similar to(0.8) "CCO"', SCHEMA);
  assert.deepEqual(root.nodes[0].options, {threshold: 0.8});
  assert.equal(Filters.parse('structure is similar to "CCO"', SCHEMA).root.nodes[0].options, undefined,
    'no parameter written: the operator\'s default applies');
});

hosted('a builder row edits the parameter after the value, into the condition\'s options', async () => {
  const fb = new FilterBuilder({schema: SCHEMA,
    value: Filters.group('and', [Filters.cond('structure', 'is similar to', 'CCO')])});
  document.body.append(fb.root);
  const parameter = fb.root.querySelector('[data-u2-part="parameter"] input');
  assert.ok(parameter !== null);
  assert.equal(parameter.value, '0.7');
  parameter.value = '0.9';
  fire(parameter, 'input');
  await flush();
  assert.deepEqual(fb.value.value.nodes[0].options, undefined, 'typed digits are not each a search');
  fire(parameter, 'blur');
  await flush();
  assert.deepEqual(fb.value.value.nodes[0].options, {threshold: 0.9}, 'leaving the box takes the number');
  fb.dispose();
});

hosted('a parameter is a number within its bounds or left out; an emptied or out-of-range one is refused', () => {
  const problems = (options) => Filters.validate(Filters.group('and',
    [Filters.cond('structure', 'is similar to', 'CCO', options === undefined ? undefined : {options})]), SCHEMA, 'domain')
    .map((p) => [p.code, p.message]);
  const refused = [['invalid-value', 'Similarity at least: a number from 0 to 1']];
  assert.deepEqual(problems(undefined), [], 'left out: the default applies');
  assert.deepEqual(problems({threshold: 0.4}), []);
  assert.deepEqual(problems({threshold: 7}), refused);
  assert.deepEqual(problems({threshold: null}), refused, 'emptied is not the default');
  assert.deepEqual(Filters.parse('structure is similar to(7) "CCO"', SCHEMA, 'domain').problems.map((p) => p.message),
    ['Similarity at least: a number from 0 to 1'], 'a link carrying one is refused too');
});

hosted('the parameter box is cleared and retyped; text that is no number in range is never applied as another', async () => {
  const fb = new FilterBuilder({schema: SCHEMA, target: 'domain',
    value: Filters.group('and', [Filters.cond('structure', 'is similar to', 'CCO', {options: {threshold: 0.7}})])});
  document.body.append(fb.root);
  const box = () => fb.root.querySelector('[data-u2-part="parameter"] input');
  const problem = () => fb.root.querySelector('[data-u2="filter-row"] [data-u2-part="problem"]');
  const threshold = () => fb.value.value.nodes[0].options.threshold;
  const applied = () => Filters.pruneInvalid(fb.value.value, SCHEMA, 'domain').nodes.length;
  const type = async (text, commit = true) => {
    box().value = text;
    fire(box(), 'input');
    if (commit)
      fire(box(), 'keydown', {key: 'Enter'});
    await flush();
  };
  assert.equal(box().closest('.u2-fb-parameter').textContent.startsWith('Similarity at least'), true,
    'the label wraps with its box');
  await type('');
  assert.equal(box().value, '', 'an emptied box stays empty: the default is not written back over it');
  assert.equal(threshold(), null);
  assert.equal(applied(), 0, 'and nothing is applied meanwhile');
  assert.equal(problem().hidden, false);
  assert.equal(problem().textContent, 'Similarity at least: a number from 0 to 1');
  await type('0', false);
  await type('0.0', false);
  assert.equal(threshold(), null, 'a number on its way is not taken before Enter');
  await type('0.0');
  assert.equal(threshold(), 0);
  await type('0.0.4', false);
  assert.equal(box().value, '0.0.4');
  assert.equal(threshold(), null, 'a mistyped number is not the last good prefix');
  assert.equal(applied(), 0);
  await type('0.4');
  assert.equal(threshold(), 0.4);
  assert.equal(applied(), 1);
  assert.equal(problem().hidden, true);
  await type('7');
  assert.equal(box().value, '7', 'the text stays for the user to fix');
  assert.equal(threshold(), null);
  assert.equal(applied(), 0, 'out of range is not applied');
  assert.equal(problem().textContent, 'Similarity at least: a number from 0 to 1');
  const node = fb.value.value.nodes[0];
  fb.value.value = {...fb.value.value, nodes: [{...node, options: {threshold: 0.5}}]};
  await flush();
  assert.equal(box().value, '0.5', 'a threshold set from outside reaches the box');
  fb.dispose();
});
