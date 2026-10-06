/* The u2-filter-builder registration: the spec renders, `value` binds two-way as the context
   signal itself, `mode` binds two-way, `query` is readable through a bind path, and the manifest
   entry carries the grammar for the LLM surface. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {Filters} from '../src/core/filter/index.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {TYPE} from 'datagrok-api/u2core';

const SCHEMA = {properties: [
  {name: 'name', type: TYPE.STRING},
  {name: 'age', type: TYPE.INT, min: 0, max: 120},
  {name: 'sex', type: TYPE.STRING, choices: ['F', 'M']},
], values: {name: ['Aspirin', 'Ibuprofen']}};

function spec(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      Filters.resetIds('');
      await body();
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function registry() {
  const reg = new Registry();
  registerAll(reg);
  return reg;
}

function pick(selectEl, value) {
  selectEl.value = value;
  fire(selectEl, 'change');
}

spec('manifest: the entry states the grammar and three examples, props carry the binding tiers', () => {
  const meta = registry().manifest().components.find((c) => c.tag === 'u2-filter-builder');
  assert.notEqual(meta, undefined);
  assert.equal(meta.category, 'Inputs');
  assert.match(meta.usage, /prop op value/);
  assert.match(meta.usage, /age > 30 and sex = "F"/);
  assert.match(meta.usage, /mw between 200 and 500/);
  assert.match(meta.usage, /-1w 2d now/);
  const props = Object.fromEntries(meta.props.map((p) => [p.name, p]));
  assert.deepEqual([props.value.type, props.value.bindable, props.value.twoWay], [TYPE.OBJECT, true, true]);
  assert.deepEqual([props.mode.bindable, props.mode.twoWay, props.mode.choices], [true, true, ['simple', 'advanced']]);
  assert.deepEqual([props.query.type, props.query.bindable, props.query.twoWay], [TYPE.STRING, true, undefined]);
  assert.equal(props.schema.bindable, undefined, 'schema is re-render tier');
  assert.deepEqual(props.orientation.choices, ['vertical', 'horizontal']);
  assert.deepEqual(meta.events, ['change']);
  assert.equal(meta.example.tag, 'u2-filter-builder');
});

spec('the example renders rows from its literal schema and value, and `query` reads back', () => {
  const reg = registry();
  const meta = reg.get('u2-filter-builder');
  const instance = renderSpec({$schema: 'dg-ui/1', root: meta.example}, new SpecContext(), reg);
  document.body.append(instance.root);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  const rows = instance.root.querySelectorAll('[data-u2-node]');
  assert.equal(rows.length, 1);
  assert.equal(rows[0].querySelector('[data-u2-part="prop"] select').value, 'age');
  assert.equal(instance.root.querySelector('[data-u2-part="query"]').textContent, 'age > 30');
  instance.dispose();
});

spec('value binds two-way: the builder adopts the context signal; DOM edits reach it and writes reach the rows', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {criteria: Filters.group('and', [Filters.cond('age', '>', 30)])}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-filter-builder', name: 'fb', props: {schema: SCHEMA}, bind: {value: '$.criteria'}},
  }, ctx, reg);
  document.body.append(instance.root);
  const fb = Control.forElement(instance.root.querySelector('[data-u2="filter-builder"]'));
  assert.equal(fb.value === ctx.data.criteria, true, 'the context signal itself');

  pick(instance.root.querySelector('[data-u2-part="op"] select'), '<');
  assert.equal(ctx.data.criteria.peek().nodes[0].operator, '<');
  ctx.data.criteria.value = Filters.group('or', [Filters.cond('sex', '=', 'M'), Filters.cond('name', 'like', 'a')]);
  const rows = instance.root.querySelectorAll('[data-u2-node]');
  assert.equal(rows.length, 2);
  assert.equal(rows[0].querySelector('[data-u2-part="value"] select').value, 'M');
  assert.equal(rows[1].querySelector('[data-u2-part="value"]').dataset.u2, 'suggest-input',
    'schema values from the literal reach the string editor');
  instance.dispose();
});

spec('mode binds two-way and query is readable through a bind path on the named node', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {m: 'simple', criteria: Filters.group('and', [Filters.cond('age', '>', 30)])}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-div-v', children: [
      {tag: 'u2-filter-builder', name: 'fb', props: {schema: SCHEMA}, bind: {value: '$.criteria', mode: '$.m'}},
      {tag: 'h3', name: 'echo', bind: {text: '$.fb.query'}},
    ]},
  }, ctx, reg);
  document.body.append(instance.root);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  const fb = Control.forElement(instance.root.querySelector('[data-u2="filter-builder"]'));
  assert.equal(fb.mode === ctx.data.m, true, 'the mode signal is adopted');
  const echo = instance.root.querySelector('h3');
  assert.equal(echo.textContent, 'age > 30', 'query flows into the bound heading');
  fb.addCondition();
  assert.equal(echo.textContent, 'age > 30 and name =');
  ctx.data.m.value = 'advanced';
  assert.equal(fb.mode.peek(), 'advanced');
  assert.equal(fb.setMode('simple'), true);
  assert.equal(ctx.data.m.peek(), 'simple', 'setMode writes the bound signal');
  const {signal, writable} = instance.resolveBinding('$.fb.query');
  assert.equal(signal.value, 'age > 30 and name =');
  assert.equal(writable, false);
  instance.dispose();
});
