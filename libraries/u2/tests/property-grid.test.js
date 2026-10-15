/* The list row (WO-10): a string_list descriptor edits through ListInput, and the element-wise
   comparison that keeps a rebuilt array from reading as an edit — identity would loop. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
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

/** What the grid reports through its `change` event. */
function watch(target) {
  const changes = [];
  target.on('change', (change) => changes.push(change));
  return changes;
}

grid('string_list: the row edits a list, and the same list rebuilt is not an edit', () => {
  const editor = new PropertyGrid();
  document.body.append(editor.root);
  editor.setProperties([{name: 'items', type: 'string_list', description: 'Crumbs'}],
    {items: ['Home', 'Projects']});
  const changes = watch(editor);
  const field = editor.root.querySelector('.u2-list-text');
  assert.equal(field.value, 'Home,Projects', 'the list is joined into the field');
  assert.equal(field.placeholder, 'Crumbs', 'the description is the placeholder');
  assert.deepEqual(changes, []);

  field.value = 'Home,Projects,Demo';
  fire(field, 'input');
  assert.deepEqual(editor.values.peek().items, ['Home', 'Projects', 'Demo']);
  assert.deepEqual(changes[0], {name: 'items', value: ['Home', 'Projects', 'Demo']});

  editor.values.value = {items: ['Home', 'Projects', 'Demo']};
  assert.equal(changes.length, 1, 'a fresh array carrying the same list reports nothing');

  editor.values.value = {items: ['Home']};
  assert.equal(field.value, 'Home', 'a real change still reaches the field');
  assert.equal(changes.length, 1, 'and is not echoed back as an edit');

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

grid('a row is captioned the platform way, and its tooltip is the description, never the name again', () => {
  const editor = new PropertyGrid();
  document.body.append(editor.root);
  editor.setProperties([{name: 'showLegend', type: 'bool'}, {name: 'opacity', type: 'double', description: 'From 0 to 1'}],
    {showLegend: true, opacity: 0.8});
  const names = [...editor.root.querySelectorAll('.u2-propgrid-name')];
  assert.deepEqual(names.map((n) => n.textContent), ['Show Legend', 'Opacity']);
  assert.deepEqual(names.map((n) => n.title), ['', 'From 0 to 1']);
  assert.deepEqual(names.map((n) => n.parentElement.dataset.u2Name), ['showLegend', 'opacity'],
    'the row keeps its property name for automation');
  assert.deepEqual(editor.getItems().list.map((x) => [x.key, x.label]), [['showLegend', 'Show Legend'], ['opacity', 'Opacity']]);
  editor.dispose();
});

grid('an invalid row keeps its editor: the box has a floor, so the full-width error line wraps under it', () => {
  const sheet = readFileSync(new URL('../css/property-grid.css', import.meta.url), 'utf8');
  const at = sheet.indexOf('.u2-propgrid-editor > .u2-input-root > .u2-input-box {');
  assert.notEqual(at, -1);
  assert.match(sheet.slice(at, sheet.indexOf('}', at)), /min-width: 60px;/);
});
