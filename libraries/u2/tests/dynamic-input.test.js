/* DynamicInput: the renderer is the editor — every write swaps the element, nothing edits back. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {DynamicInput} from '../src/components/dynamic-input.js';

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
  return input.root.querySelector('.u2-dynamic');
}

function chip(text) {
  const el = document.createElement('span');
  el.className = 'demo-chip';
  el.textContent = text;
  return el;
}

smoke('a write re-renders, and only the newest element is in the editor', async () => {
  const rendered = [];
  const input = mount(new DynamicInput({label: 'Object', value: 'Aspirin',
    render: (value) => {
      rendered.push(value);
      return value === null ? null : chip(String(value));
    }}));

  assert.deepEqual(rendered, ['Aspirin']);
  assert.equal(editor(input).textContent, 'Aspirin');
  assert.equal(editor(input).querySelectorAll('.demo-chip').length, 1);

  input.value.value = 'Ibuprofen';
  assert.deepEqual(rendered, ['Aspirin', 'Ibuprofen']);
  assert.equal(editor(input).querySelectorAll('.demo-chip').length, 1, 'the old element is gone');
  assert.equal(editor(input).textContent, 'Ibuprofen');
  input.dispose();
});

smoke('null renders empty, and the editor fills back in on the next value', async () => {
  const input = mount(new DynamicInput({label: 'Object',
    render: (value) => value === null ? null : chip(String(value))}));
  assert.equal(editor(input).childNodes.length, 0);
  assert.equal(input.value.value, null, 'null is the default value');

  input.value.value = 42;
  assert.equal(editor(input).textContent, '42');
  input.value.value = null;
  assert.equal(editor(input).childNodes.length, 0);
  input.dispose();
});

smoke('dispose removes the rendered element and releases the effect', async () => {
  const before = Scope.liveCount;
  let renders = 0;
  const input = mount(new DynamicInput({label: 'Object', value: 'x',
    render: (value) => {
      renders++;
      return chip(String(value));
    }}));
  const element = editor(input).firstChild;
  assert.equal(element.parentNode, editor(input));

  input.dispose();
  assert.equal(element.parentNode, null, 'the rendered element is detached');
  assert.equal(Scope.liveCount, before);

  input.value.value = 'y';
  assert.equal(renders, 1, 'the render effect died with the scope');
});

smoke('size fixes the editor box, and the input carries its own data-u2 id', async () => {
  const input = mount(new DynamicInput({label: 'Preview', size: {width: 200, height: 80},
    render: () => chip('preview')}));
  assert.equal(input.root.dataset.u2, 'dynamic-input');
  assert.equal(editor(input).style.width, '200px');
  assert.equal(editor(input).style.height, '80px');
  input.dispose();
});

smoke('a validator still applies — the value is read-only, not unchecked', async () => {
  const input = mount(new DynamicInput({label: 'Object', render: () => null}));
  input.addValidator((value) => value === null ? 'Nothing selected' : null);
  assert.equal(input.validity.value, 'Nothing selected');

  input.value.value = 'Aspirin';
  assert.equal(input.validity.value, null);
  input.dispose();
});
