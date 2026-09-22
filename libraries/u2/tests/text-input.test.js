/* TextInput variants: the password eye and the auto-resizing field. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, TextInput, signal, batch} from '../src/index.js';

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

function field(input) {
  return input.root.querySelector('input');
}

smoke('password: masked by default, the eye toggles the field both ways', async () => {
  const input = mount(new TextInput({label: 'Password', password: true, value: 'hunter2'}));
  const eye = input.root.querySelector('.u2-input-eye');
  const glyph = eye.querySelector('.u2-icon');
  assert.equal(field(input).type, 'password');
  assert.equal(eye.getAttribute('aria-label'), 'Show password');
  assert.equal(glyph.classList.contains('fa-eye'), true);

  fire(eye, 'click');
  assert.equal(field(input).type, 'text');
  assert.equal(eye.getAttribute('aria-label'), 'Hide password');
  assert.equal(glyph.classList.contains('fa-eye-slash'), true);
  assert.equal(glyph.classList.contains('fa-eye'), false);

  fire(eye, 'click');
  assert.equal(field(input).type, 'password');
  assert.equal(glyph.classList.contains('fa-eye'), true);

  field(input).value = 'hunter3';
  fire(field(input), 'input');
  assert.equal(input.value.value, 'hunter3', 'the value binding is the plain one');
  input.dispose();
});

smoke('autoResize: the width follows the measured text between minWidth and maxWidth', async () => {
  const input = mount(new TextInput({label: 'Name', autoResize: true, placeholder: 'Compound'}));
  const measure = input.root.querySelector('.u2-input-measure');
  assert.equal(measure.textContent, 'Compound', 'an empty field is measured by its placeholder');
  assert.equal(field(input).style.width, '100px', 'never below minWidth');

  measure.offsetWidth = 150;
  input.value.value = 'Acetylsalicylic acid';
  assert.equal(measure.textContent, 'Acetylsalicylic acid');
  assert.equal(field(input).style.width, '160px', 'text width plus the caret allowance');

  measure.offsetWidth = 600;
  input.value.value = 'N-(4-hydroxyphenyl)acetamide, analytical standard, 99%';
  assert.equal(field(input).style.width, '300px', 'never above maxWidth');

  const custom = mount(new TextInput({label: 'Name', autoResize: true, minWidth: 40, maxWidth: 60}));
  custom.root.querySelector('.u2-input-measure').offsetWidth = 200;
  custom.value.value = 'long';
  assert.equal(custom.root.querySelector('input').style.width, '60px');
  custom.dispose();
  input.dispose();
});

smoke('the plain and search variants are untouched by the new options', async () => {
  const plain = mount(new TextInput({label: 'Name', value: 'Aspirin'}));
  assert.equal(plain.root.querySelector('.u2-input-editor').tagName, 'INPUT');
  assert.equal(plain.root.querySelector('.u2-input-measure'), null);
  assert.equal(plain.root.querySelector('.u2-input-eye'), null);
  plain.dispose();

  const search = mount(new TextInput({label: 'Filter', search: true, value: 'ac'}));
  const clear = search.root.querySelector('.u2-input-clear');
  assert.equal(clear.hidden, false);
  fire(clear, 'click');
  assert.equal(search.value.value, '');
  assert.equal(clear.hidden, true);
  search.dispose();
});

smoke('follow: shows the model value, never echoes it to onChanged, user commits still reach it', async () => {
  const name = signal('orders');
  const changes = [];
  const input = mount(new TextInput({label: 'Name', value: '', commitOn: 'change', onChanged: (v) => changes.push(v)}))
    .follow(() => name.value);
  assert.equal(field(input).value, 'orders');

  name.value = 'order_lines';
  batch(() => name.value = 'lines');
  assert.equal(field(input).value, 'lines');
  assert.deepEqual(changes, []);

  field(input).value = 'items';
  fire(field(input), 'change');
  assert.deepEqual(changes, ['items']);
  input.dispose();
});

smoke('follow: a re-read with an unchanged value leaves the field being typed in alone', async () => {
  const revision = signal(0);
  const input = mount(new TextInput({label: 'Name', value: '', commitOn: 'change'}))
    .follow(() => (revision.value, 'orders'));
  const editor = field(input);
  editor.value = 'ord';
  revision.value++;
  assert.equal(field(input), editor);
  assert.equal(editor.value, 'ord');
  assert.equal(input.value.value, 'orders');
  input.dispose();
});
