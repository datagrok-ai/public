import {test as nodeTest} from 'node:test';
import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {flush} from './dom-shim.js';
import {Scope, Form, TextInput, ChoiceChip, button, computed, signal} from '../src/index.js';
import {Usage} from '../src/core/scope.js';
import {Section} from '../src/components/containers/section.js';
import {dfBindings} from '../src/sources/df-bindings.js';
import {DataFrame} from './platform-doubles.mjs';
import {u2test, mount, user} from '../src/testing/index.js';

function reportOnly(body) {
  const was = Usage.mode;
  Usage.mode = 'report';
  try {
    return body();
  } finally {
    Usage.mode = was;
    Usage.clear();
  }
}

test('disabled: false is on, set or given; disabled: true is reported (reason-boolean) and says something', () => {
  const given = new TextInput({label: 'Name', disabled: false});
  assert.equal(given.disabled, null);
  assert.equal(given.root.querySelector('input').disabled, false);
  assert.equal(given.validates, true, 'its form still counts it');
  const set = new TextInput({label: 'Name'});
  set.disabled = false;
  assert.equal(set.disabled, null);
  assert.equal(set.root.querySelector('input').disabled, false);
  assert.throws(() => new TextInput({label: 'Name', disabled: true}), {code: 'reason-boolean'});
  assert.throws(() => button('Save', () => {}, {disabled: true}), {code: 'reason-boolean'});
  const off = reportOnly(() => new TextInput({label: 'Name', disabled: true}));
  assert.equal(off.disabled, 'Unavailable', 'reported only, it is off and says so — never "true"');
});

test('button.disabled = true with nothing said is reported (reason-missing); a button that says why is not', () => {
  const bare = button('Save', () => {});
  assert.throws(() => bare.disabled = true, {code: 'reason-missing'});
  const why = button('Submit', () => {}, {disabled: 'Name it first'});
  why.disabled = false;
  assert.doesNotThrow(() => why.disabled = true);
  assert.equal(why.getAttribute('aria-disabled'), 'true');
});

test('an input bound to a computed is reported (bind-readonly) and, reported only, shown read-only', () => {
  const first = signal('Ada');
  assert.throws(() => new TextInput({label: 'Full name', bind: computed(() => `${first.value} Lovelace`)}),
    {code: 'bind-readonly'});
  const shown = new TextInput({label: 'Full name', bind: computed(() => first.value), readOnly: true});
  assert.equal(shown.readOnly, true, 'readOnly: true says it is shown, not edited: no report');
  const quiet = reportOnly(() => new TextInput({label: 'Full name', bind: computed(() => first.value)}));
  assert.equal(quiet.readOnly, true, 'no edit that cannot land');
  assert.equal(quiet.root.querySelector('input').readOnly, true);
});

test('a required field out of reach holds no form: in a section off with a reason, or under a host hidden later',
  async () => {
    const off = signal('Turn on Advanced to edit these');
    const form = new Form();
    form.add(new TextInput({label: 'Name', value: 'run 1'}));
    const advanced = new Section({title: 'Advanced', disabled: off});
    form.add(new TextInput({label: 'Threshold', required: true}), advanced.body);
    form.root.append(advanced.root);
    const host = document.createElement('div');
    form.add(new TextInput({label: 'Extra', value: 'x', required: true}), host);
    form.root.append(host);
    document.body.append(form.root);
    try {
      assert.equal(form.validity.value, null, 'the section is inert: nothing in it can be filled');
      off.value = null;
      await flush();
      assert.notEqual(form.validity.value, null, 'within reach again, the empty Threshold counts');
      form.inputs[1].value = '0.5';
      assert.equal(form.validity.value, null);
      form.inputs[2].value = '';
      assert.notEqual(form.validity.value, null);
      host.hidden = true;
      await flush();
      assert.equal(form.validity.value, null, 'a host hidden by hand is seen, not only on some next read');
    } finally {
      form.root.remove();
    }
  });

test('the reach watch wakes for a hidden or inert change on an input or around one, never for one elsewhere',
  async () => {
    const input = new TextInput({label: 'Name'});
    const around = document.createElement('div');
    const elsewhere = document.createElement('div');
    around.append(input.root);
    document.body.append(around, elsewhere);
    try {
      let runs = 0;
      const judged = computed(() => {
        runs++;
        return input.validates;
      });
      assert.equal(judged.value, true);
      elsewhere.hidden = true;
      elsewhere.setAttribute('inert', '');
      await flush();
      assert.equal(judged.value, true);
      assert.equal(runs, 1, 'a status line or a view hidden away from every input judges nothing again');
      input.root.hidden = true;
      await flush();
      assert.equal(judged.value, false, 'the input\'s own root hidden by hand counts');
      input.root.hidden = false;
      around.setAttribute('inert', '');
      await flush();
      assert.equal(judged.value, false, 'and so does an element around it');
      assert.equal(runs, 3);
    } finally {
      around.remove();
      elsewhere.remove();
    }
  });

nodeTest('owner.subscribe takes a signal from plain JS: the handler follows it until the owner goes, then nothing', () => {
  const errors = [];
  const was = Scope.onUnhandledError;
  Scope.onUnhandledError = (e) => errors.push(e);
  try {
    const s = signal(1);
    const seen = [];
    const owner = Scope.root('subscribes a signal');
    owner.subscribe(s, (v) => seen.push(v));
    s.value = 2;
    owner.dispose();
    s.value = 3;
    assert.deepEqual(seen, [1, 2]);
    assert.deepEqual(errors, [], 'dispose unsubscribes without a throw');
    assert.equal(s._targets, undefined, 'nothing is left subscribed to the signal');
  } finally {
    Scope.onUnhandledError = was;
  }
});

class Day {
  constructor(ms) { this.ms = ms; }
  valueOf() { return this.ms; }
}

/** A frame whose datetime cells answer a fresh wrapper per read, as `DG.DataFrame.get` does (a new `dayjs`). */
class FreshDates extends DataFrame {
  get(column, row) {
    const v = super.get(column, row);
    return column === 'when' && v != null ? new Day(v) : v;
  }

  set(column, row, value) {
    if (this.dart.writes.length >= 20)
      throw new Error('written back in a loop');
    super.set(column, row, column === 'when' && value != null ? +value : value);
  }
}

nodeTest('currentRow over a datetime cell that reads fresh each time: binding and moving write nothing, an edit once',
  () => {
    const df = new FreshDates([{name: 'name', type: 'string'}, {name: 'when', type: 'datetime'}],
      [{name: 'a', when: 1700000000000}, {name: 'b', when: 1700000100000}]);
    const scope = Scope.root('fresh dates');
    const errors = [];
    scope.onError = (e) => errors.push(e);
    const when = dfBindings(signal(df), scope).bindStep('currentRow').bindStep('when');
    assert.equal(+when.value, 1700000000000);
    assert.deepEqual(df.dart.writes, [], 'binding reads, never writes');
    df.currentRowIdx = 1;
    assert.equal(+when.value, 1700000100000);
    assert.deepEqual(df.dart.writes, [], 'nor does moving the row');
    when.value = new Day(1800000000000);
    assert.deepEqual(df.dart.writes, [['when', 1, 1800000000000]], 'the user\'s edit lands once');
    assert.deepEqual(errors, []);
    scope.dispose();
  });

nodeTest('choiceChip: re-picking the value it holds is no input; picking another is', u2test(() => {
  const picked = [];
  const chip = mount(() => new ChoiceChip({label: 'R1', items: ['H', 'OH'], value: 'H',
    onInput: (v) => picked.push(v)}));
  user.click(chip);
  user.pick(chip, 'H');
  assert.deepEqual(picked, []);
  user.click(chip);
  user.pick(chip, 'OH');
  assert.deepEqual(picked, ['OH']);
}));
