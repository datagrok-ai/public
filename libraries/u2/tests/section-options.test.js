import {test} from 'node:test';
import assert from 'node:assert/strict';
import './dom-shim.js';
import {Section, TextInput, iconButton, signal} from '../src/index.js';
import {u2test, mount, byName, expectOperable, user} from '../src/testing/index.js';

test('section: actions after the title act without folding it', u2test(() => {
  let refreshed = 0;
  const section = mount(() => new Section({title: 'Campaigns',
    actions: [iconButton('sync', () => refreshed++, {tooltip: 'Refresh'})]}));
  const actions = section.header.querySelector('[data-u2-part="actions"]');
  assert.notEqual(actions, null);
  user.click(byName(section, 'Refresh'));
  assert.equal(refreshed, 1);
  assert.equal(section.expanded.value, true, 'the press was the action\'s, not the header\'s');
  user.press(byName(section, 'Refresh'), 'Enter');
  assert.equal(section.expanded.value, true);
  user.click(section.header);
  assert.equal(section.expanded.value, false, 'the header still folds');
  expectOperable(section.root);
}));

test('section: disabled dims its body out of reach and says why after the title', u2test(() => {
  const off = signal('Pick a template first');
  const section = mount(() => new Section({title: 'Fields', disabled: off}));
  const name = new TextInput({label: 'Name', name: 'name'});
  section.add(name);
  const reason = section.header.querySelector('[data-u2-part="reason"]');
  assert.equal(reason.textContent, 'Pick a template first');
  assert.equal(section.body.hasAttribute('inert'), true);
  assert.equal(section.body.getAttribute('aria-description'), 'Pick a template first');
  assert.throws(() => user.type(byName(section, 'name'), 'KHL'), {message: /inert: Pick a template first/});
  expectOperable(section.root);
  off.value = null;
  assert.equal(reason.textContent, '');
  assert.equal(section.body.hasAttribute('inert'), false);
  user.type(byName(section, 'name'), 'KHL');
  assert.equal(name.value, 'KHL');
}));
