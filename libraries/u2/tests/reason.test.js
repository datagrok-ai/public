import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {platform} from './platform-doubles.mjs';
import {Reason, Scope, TextInput, Menu, renderSpec, button, signal} from '../src/index.js';
import {PlatformInput} from '../src/dg/inputs/from-dart-input.js';

register('./dg-stub.mjs', import.meta.url);
const {dartInputFor, DartInput} = await import('../src/dg/inputs/input-bridge.js');

function off(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      platform.reset();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

off('Reason: null and \'\' are on; a signal is read now or followed; apply writes the reason with the state', () => {
  assert.equal(Reason.now(null), null);
  assert.equal(Reason.now(''), null);
  assert.equal(Reason.now(signal('Pick one')), 'Pick one');
  const el = document.createElement('div');
  Reason.apply(el, 'Pick a source file first');
  assert.equal(el.getAttribute('aria-disabled'), 'true');
  assert.equal(Reason.of(el), 'Pick a source file first');
  Reason.apply(el, '');
  assert.equal(el.hasAttribute('aria-disabled'), false);
  assert.equal(el.hasAttribute('aria-description'), false);
  assert.equal(Reason.of(el), null);
});

off('an input off with a reason: greyed, every control off, the reason on its root and in its state', () => {
  const owner = Scope.root('page');
  const reason = signal('Pick a source file first');
  const input = owner.build(() => new TextInput({label: 'Merge by', disabled: reason}));
  assert.equal(input.disabled, 'Pick a source file first');
  assert.equal(input.root.classList.contains('u2-input-disabled'), true);
  assert.equal(input.root.querySelector('input').disabled, true);
  assert.equal(input.getState().disabledReason, 'Pick a source file first');
  assert.equal(input.validates, false, 'an input that is off fails no form');
  reason.value = null;
  assert.equal(input.disabled, null);
  assert.equal(input.root.querySelector('input').disabled, false);
  assert.equal(input.getState().disabledReason, undefined);
  owner.dispose();
});

off('a menu item off says why, and a press on it runs nothing', () => {
  const ran = [];
  const menu = new Menu().item('Paste', () => ran.push('paste'), {disabled: 'The clipboard is empty'});
  menu.show({x: 0, y: 0});
  const item = document.querySelector('[data-u2-item="Paste"]');
  assert.equal(item.getAttribute('aria-disabled'), 'true');
  assert.equal(item.getAttribute('aria-description'), 'The clipboard is empty');
  item.click();
  assert.deepEqual(ran, []);
  menu.close();
});

off('a spec that still says `enabled` does not parse, naming `disabled`', () => {
  for (const node of [{tag: 'u2-text-input', name: 'name', props: {enabled: false}},
    {tag: 'u2-text-input', bind: {enabled: '$.on'}}]) {
    assert.throws(() => renderSpec({$schema: 'dg-ui/1', root: {tag: 'div', children: [node]}}),
      /prop "enabled" is now "disabled": the reason the user reads/);
  }
});

off('a platform input the platform turned off is off in u2 too, and u2\'s reason turns the platform\'s off', () => {
  const root = document.createElement('div');
  root.className = 'ui-input-root ui-input-disabled';
  root.setAttribute('data-widget', 'true');
  root.append(Object.assign(document.createElement('input'), {className: 'ui-input-editor'}));
  const dgInput = {root, caption: 'Name', value: '', enabled: false,
    onChanged: {subscribe: () => ({unsubscribe: () => {}})}};
  const input = new PlatformInput(dgInput);
  assert.equal(input.disabled, PlatformInput.NOT_HERE);
  input.disabled = null;
  assert.equal(dgInput.enabled, true, 'on in u2, on in the platform');
  input.disabled = 'Pick a table first';
  assert.equal(dgInput.enabled, false);
  input.dispose();
});

off('an editor the platform turned off says why: its tooltip, else "Not available here"', () => {
  const dart = dartInputFor(() => new TextInput({label: 'Name'}));
  dart.root.classList.add('ui-input-disabled');
  const box = dart.getInput();
  assert.equal(Reason.of(box), DartInput.NOT_HERE);
  dart.detach();
});

off('a button takes no boolean: off is a reason', () => {
  const owner = Scope.root('page');
  const save = owner.build(() => button('Save', () => {}, {disabled: 'Name it first'}));
  assert.equal(save.disabled, true);
  assert.equal(Reason.of(save), 'Name it first');
  owner.dispose();
});
