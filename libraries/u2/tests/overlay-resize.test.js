/* Dialog, overlay, type-ahead and filter-builder fixes from the U2Demo visual review: the platform
   dialog's footer order, OK gated on the body's form, a shown dialog kept on screen through window
   resizes, the error row cleared by an edit, anchored popups closed by a window resize, stray
   type-ahead text dropped on blur, and a new filter row quiet about its missing value until touched. */
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, Dialog, Form, TextInput, TypeAhead, Overlay, Filters, FilterBuilder, button} from '../src/index.js';
import {TYPE} from 'datagrok-api/u2core';

const mounted = [];

function smoke(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      Filters.resetIds('');
      await body();
    } finally {
      for (const control of mounted.splice(0))
        control.dispose();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function mount(control) {
  mounted.push(control);
  document.body.append(control.root);
  return control;
}

const footer = (dialog) => [...dialog.root.querySelectorAll('.u2-dialog-buttons button')]
  .map((b) => b.textContent.trim());

smoke('the footer reads as the platform\'s: CANCEL first, OK right-most and flat, other buttons between', () => {
  const dialog = Dialog.create('Order').onOK(() => {}).onCancel(() => {}).addButton('APPLY', () => {}).show();
  assert.deepEqual(footer(dialog), ['CANCEL', 'APPLY', 'OK']);
  assert.equal(dialog.getButton('OK').classList.contains('u2-btn-primary'), false);
  dialog.close();
  const late = Dialog.create('Cancel first').onCancel(() => {}).onOK(() => {}).show();
  assert.deepEqual(footer(late), ['CANCEL', 'OK']);
  late.close();
});

smoke('OK follows the validity of the form in the body, Enter included', async () => {
  const name = new TextInput({label: 'Name', value: 'u2'});
  name.addValidator((v) => v.trim() ? null : 'Required');
  let accepted = 0;
  const dialog = Dialog.create('Gate', {keepAlive: true}).add(new Form().add(name)).onOK(() => accepted++).show();
  const ok = dialog.getButton('OK');
  assert.equal(ok.disabled, false);
  name.value.value = '';
  await flush();
  assert.equal(ok.disabled, true, 'an invalid form disables OK');
  fire(dialog.root.querySelector('input'), 'keydown', {key: 'Enter'});
  assert.equal(accepted, 0, 'Enter does not accept it either');
  assert.equal(dialog.isOpen.value, true);
  name.value.value = 'fixed';
  await flush();
  assert.equal(ok.disabled, false);
  fire(ok, 'click');
  assert.equal(accepted, 1);
  dialog.dispose();
});

smoke('a shown dialog comes back on screen when the window shrinks; a centred one stays centred', () => {
  const dialog = Dialog.create('Resize').onOK(() => {}).show();
  dialog.root.rect = new DOMRect(0, 0, 200, 100);
  window.innerWidth = 1600;
  window.innerHeight = 1000;
  fire(window, 'resize', {bubbles: false});
  assert.equal(dialog.root.style.left, '700px');
  assert.equal(dialog.root.style.top, '450px');
  window.innerWidth = 500;
  window.innerHeight = 300;
  fire(window, 'resize', {bubbles: false});
  assert.equal(dialog.root.style.left, '150px');
  assert.equal(dialog.root.style.top, '100px');
  dialog.close();
  window.innerWidth = 1024;
  window.innerHeight = 768;
});

smoke('an edit in the body clears the error row of a failed OK', () => {
  const dialog = Dialog.create('Fail', {keepAlive: true}).add(new TextInput({label: 'Name'}))
    .onOK(() => { throw new Error('Server said no'); }).show();
  fire(dialog.getButton('OK'), 'click');
  assert.equal(dialog.error.value, 'Server said no');
  fire(dialog.root.querySelector('input'), 'input');
  assert.equal(dialog.error.value, null);
  dialog.dispose();
});

smoke('a busy OK clicked with the mouse keeps the focus inside the dialog', async () => {
  let release;
  const dialog = Dialog.create('Busy', {keepAlive: true}).add(new TextInput({label: 'Name'}))
    .onOK(() => new Promise((resolve) => release = resolve)).show();
  const ok = dialog.getButton('OK');
  ok.focus();
  fire(ok, 'click');
  assert.equal(dialog.busy.value, true);
  assert.equal(document.activeElement, dialog.root);
  release(false);
  await flush();
  dialog.dispose();
});

smoke('a window resize closes the anchored popups, not the dialogs', () => {
  const anchor = button('Open', () => {});
  document.body.append(anchor);
  const scope = new Scope();
  const popup = document.createElement('div');
  Overlay.show(anchor, popup, scope);
  const dialog = Dialog.create('Stays').onOK(() => {}).show();
  assert.equal(popup.isConnected, true);
  fire(window, 'resize', {bubbles: false});
  assert.equal(popup.isConnected, false, 'the popup closed');
  assert.equal(dialog.isOpen.value, true, 'the dialog did not');
  dialog.close();
  scope.dispose();
});

smoke('a type-ahead drops text that is not a pick on blur, and does not spellcheck', async () => {
  const picker = mount(new TypeAhead({source: ['All users', 'Developers'], itemText: (s) => s}));
  const input = picker.root.querySelector('input');
  assert.equal(input.spellcheck, false);
  picker.selected.value = 'All users';
  input.focus();
  input.value = 'zzzz-nobody';
  fire(input, 'input');
  assert.equal(picker.selected.value, null);
  input.blur();
  await flush();
  assert.equal(input.value, '', 'the stray text is gone, so the box reads as empty as the value');
});

smoke('a new filter row says its value is missing only once it was left, or the builder validated', async () => {
  const schema = Filters.schema([{name: 'name', type: TYPE.STRING}]);
  const fb = mount(new FilterBuilder({schema}));
  const id = fb.addCondition();
  const row = fb.root.querySelector(`[data-u2-node="${id}"]`);
  const problem = () => row.querySelector('[data-u2-part="problem"]');
  assert.equal(problem().hidden, true, 'fresh rows are quiet');
  assert.equal(row.classList.contains('u2-fb-invalid'), false);
  assert.notEqual(fb.validity.peek(), null, 'the verdict itself stands');
  assert.equal(problem().parentElement.classList.contains('u2-fb-values'), true, 'the line sits under the value');
  fire(row.querySelector('.u2-fb-values input'), 'focusout');
  await flush();
  assert.equal(problem().hidden, false, 'a row left empty says so');
  const other = fb.addCondition();
  assert.equal(fb.root.querySelector(`[data-u2-node="${other}"] [data-u2-part="problem"]`).hidden, true);
  fb.touch();
  await flush();
  assert.equal(fb.root.querySelector(`[data-u2-node="${other}"] [data-u2-part="problem"]`).hidden, false,
    'touch shows every row');
});
