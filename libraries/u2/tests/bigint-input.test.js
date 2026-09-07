/* BigIntInput: text in, an exact bigint out — however many digits, with no double in between. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {BigIntInput} from '../src/components/inputs/bigint-input.js';

function smoke(name, body) {
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

function mount(component) {
  document.body.append(component.root);
  return component;
}

function editor(input) {
  return input.root.querySelector('input');
}

function type(input, text) {
  const el = editor(input);
  el.value = text;
  fire(el, 'input');
}

const HUGE = '123456789012345678901234567890';

test('parse: blank is null, an integer is exact, anything else is undefined', () => {
  assert.equal(BigIntInput.parse(''), null);
  assert.equal(BigIntInput.parse('   '), null);
  assert.equal(BigIntInput.parse(HUGE), BigInt(HUGE));
  assert.equal(BigIntInput.parse(' -42 '), -42n);
  assert.equal(BigIntInput.parse('+7'), 7n);
  assert.equal(BigIntInput.parse('1.5'), undefined);
  assert.equal(BigIntInput.parse('1e3'), undefined, 'exponents are not integer text');
  assert.equal(BigIntInput.parse('0x10'), undefined, 'and neither is hex, which BigInt() would take');
  assert.equal(BigIntInput.parse('abc'), undefined);
});

smoke('a thirty-digit number round-trips through the editor unchanged', async () => {
  const input = mount(new BigIntInput({label: 'Molecules'}));
  assert.equal(input.value.value, null);

  type(input, HUGE);
  assert.equal(input.value.value, BigInt(HUGE));
  assert.equal(String(input.value.value), HUGE, 'every digit survives');
  assert.equal(input.validity.value, null);
  assert.notEqual(Number(HUGE).toString(), HUGE, 'which a double could not have done');

  input.value.value = BigInt(HUGE) + 1n;
  assert.equal(editor(input).value, '123456789012345678901234567891');
  input.dispose();
});

smoke('non-integer text stays on screen and marks the input invalid', async () => {
  const input = mount(new BigIntInput({label: 'Molecules', value: 7n}));
  assert.equal(editor(input).value, '7');

  type(input, '1.5');
  assert.equal(input.validity.value, 'Not an integer');
  assert.equal(editor(input).value, '1.5', 'the text is not taken away');
  assert.equal(input.value.value, 7n, 'the value keeps the last good number');

  fire(editor(input), 'blur');
  assert.equal(editor(input).value, '1.5', 'blur does not overwrite text that is still invalid');

  type(input, '');
  assert.equal(input.value.value, null, 'blank clears the value, as in Dart (big_int_input.dart:29)');
  assert.equal(input.validity.value, null);
  input.dispose();
});

smoke('blur normalises the committed text, and onChanged reports every commit', async () => {
  const changes = [];
  const input = mount(new BigIntInput({label: 'Count', onChanged: (v) => changes.push(String(v))}));
  type(input, ' +007 ');
  assert.equal(input.value.value, 7n);
  fire(editor(input), 'blur');
  assert.equal(editor(input).value, '7');

  input.value.value = -1n;
  assert.equal(editor(input).value, '-1');
  assert.deepEqual(changes, ['7', '-1']);
  input.dispose();
});

smoke('two editors on one signal, and dispose kills the listeners', async () => {
  const before = Scope.liveCount;
  const first = mount(new BigIntInput({label: 'From'}));
  const second = mount(new BigIntInput({label: 'Also from', bind: first.value}));

  type(first, '900719925474099100');
  assert.equal(editor(second).value, '900719925474099100');
  assert.equal(second.value.value, 900719925474099100n);

  second.dispose();
  first.dispose();
  assert.equal(Scope.liveCount, before);

  type(first, '5');
  assert.equal(first.value.value, 900719925474099100n, 'listeners died with the scope');
});
