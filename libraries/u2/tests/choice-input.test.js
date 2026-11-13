import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, Usage} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {ChoiceInput, MultiChoiceInput} from '../src/components/inputs/choice-input.js';

function choices(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** Runs `body` in production's mode, where a misuse is reported and u2 goes on. */
function inProduction(body) {
  const mode = Usage.mode;
  Usage.mode = 'report';
  try {
    body();
  } finally {
    Usage.mode = mode;
  }
}

const nullItem = {name: 'U2UsageError', code: 'item-null'};
const reported = () => Usage.reports().filter((r) => r.code === 'item-null').length;
const options = (input) => [...input.root.querySelectorAll('option')].map((o) => o.value);
const boxes = (input) => [...input.root.querySelectorAll('input[type=checkbox]')].map((b) => b.value);

choices('ChoiceInput: a null item throws in strict mode; in production it is reported and skipped, and the rest ' +
  'render — built, set, and from a signal', () => {
  const strict = new ChoiceInput({label: 'Table', items: ['demog']});
  assert.throws(() => strict.setItems(['cars', null]), nullItem);
  strict.dispose();

  inProduction(() => {
    const before = reported();
    const input = new ChoiceInput({label: 'Table', items: ['demog', undefined, 'cars', null]});
    assert.deepEqual(options(input), ['', 'demog', 'cars']);
    assert.deepEqual(input.getItems().list.map((i) => i.key), ['demog', 'cars']);
    assert.equal(reported(), before + 1, 'said once, where it was built');

    input.setItems([null, 'geo', {value: 'refs', label: 'Refs'}]);
    assert.deepEqual(options(input), ['', 'geo', 'refs']);

    const items = signal(['a']);
    const live = new ChoiceInput({label: 'Live', items});
    items.value = ['b', undefined];
    assert.deepEqual(options(live), ['', 'b']);
    assert.ok(reported() > before + 1);

    const select = input.root.querySelector('select');
    select.value = 'geo';
    fire(select, 'change');
    assert.equal(input.getState().value, 'geo', 'and picks as ever');
    input.dispose();
    live.dispose();
  });
});

choices('MultiChoiceInput: a null item throws in strict mode; in production it is reported and skipped, and the ' +
  'checked ones that remain stay checked', () => {
  const strict = new MultiChoiceInput({label: 'Tables', items: ['demog']});
  assert.throws(() => strict.setItems([undefined, 'cars']), nullItem);
  strict.dispose();

  inProduction(() => {
    const before = reported();
    const input = new MultiChoiceInput({label: 'Tables', items: [undefined, 'demog', 'cars'], value: ['cars']});
    assert.deepEqual(boxes(input), ['demog', 'cars']);
    input.setItems(['cars', null, 'geo']);
    assert.deepEqual(boxes(input), ['cars', 'geo']);
    assert.deepEqual(input.getState().value, ['cars']);
    assert.ok(reported() > before);
    input.dispose();
  });
});
