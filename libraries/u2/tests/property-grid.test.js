/* The list row (WO-10): a string_list descriptor edits through ListInput, and the element-wise
   comparison that keeps a rebuilt array from reading as an edit — identity would loop. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, PropertyGrid} from '../src/index.js';

function grid(name, body) {
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

/** Counts what the grid reports, the initial (null) reading included. */
function watch(scope, target) {
  const changes = [];
  scope.effect(() => changes.push(target.onChanged.value));
  return changes;
}

grid('string_list: the row edits a list, and the same list rebuilt is not an edit', () => {
  const scope = new Scope();
  const editor = new PropertyGrid();
  document.body.append(editor.root);
  editor.setProperties([{name: 'items', type: 'string_list', description: 'Crumbs'}],
    {items: ['Home', 'Projects']});
  const changes = watch(scope, editor);
  const field = editor.root.querySelector('.u2-list-text');
  assert.equal(field.value, 'Home,Projects', 'the list is joined into the field');
  assert.equal(field.placeholder, 'Crumbs', 'the description is the placeholder');
  assert.deepEqual(changes, [null]);

  field.value = 'Home,Projects,Demo';
  fire(field, 'input');
  assert.deepEqual(editor.values.peek().items, ['Home', 'Projects', 'Demo']);
  assert.deepEqual(changes[1], {name: 'items', value: ['Home', 'Projects', 'Demo']});

  editor.values.value = {items: ['Home', 'Projects', 'Demo']};
  assert.equal(changes.length, 2, 'a fresh array carrying the same list reports nothing');

  editor.values.value = {items: ['Home']};
  assert.equal(field.value, 'Home', 'a real change still reaches the field');
  assert.equal(changes.length, 2, 'and is not echoed back as an edit');

  scope.dispose();
  editor.dispose();
});

grid('string_list: anything that is not a list of strings edits as the empty list', () => {
  const editor = new PropertyGrid();
  document.body.append(editor.root);
  editor.setProperties([{name: 'items', type: 'string_list'}], {items: null});
  const field = editor.root.querySelector('.u2-list-text');
  assert.equal(field.value, '', 'null shows as an empty list, not as "null"');

  field.value = 'a';
  fire(field, 'input');
  assert.deepEqual(editor.values.peek().items, ['a']);

  editor.setProperties([{name: 'items', type: 'string_list', readonly: true}], {items: ['a', 'b']});
  assert.equal(editor.root.querySelector('.u2-propgrid-value').textContent, 'a,b');
  editor.dispose();
});
