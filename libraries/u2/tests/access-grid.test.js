/* AccessGrid: inherited rows read-only above the editable ones, a checkbox change rewrites the
   row, remove drops it, the picker adds a principal with the default grant and leaves out the
   ones already there, a locked capability cannot be granted; registered as u2-access-grid. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, signal} from '../src/index.js';
import {AccessGrid} from '../src/components/forms/access-grid.js';
import {Registry} from '../src/spec/registry.js';
import {registerAll} from '../src/spec/registrations.js';

function ui(name, body) {
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

const CAPS = [{name: 'view', label: 'View'}, {name: 'edit', label: 'Edit'}, {name: 'delete', label: 'Delete'}];

function grid(options = {}) {
  const g = new AccessGrid({label: 'Access', capabilities: CAPS, principals: ['Sales', 'Developers', 'Chemists'],
    inherited: [{principal: 'You (creator)', can: {view: true, edit: true, delete: true}, from: ''}],
    value: [{principal: 'Sales', can: {view: true, edit: false, delete: false}}], ...options});
  document.body.append(g.root);
  return g;
}

const rows = (g) => g.root.querySelectorAll('tbody tr');
const boxes = (tr) => tr.querySelectorAll('.u2-access-grid-check');

ui('renders inherited rows read-only above the editable ones, headers from the capabilities', () => {
  const g = grid();
  assert.equal(g.root.dataset.u2, 'access-grid');
  assert.deepEqual(g.root.querySelectorAll('th').map((th) => th.textContent), ['Group', 'View', 'Edit', 'Delete', '']);
  const [creator, sales] = rows(g);
  assert.equal(creator.classList.contains('u2-access-grid-inherited'), true);
  assert.deepEqual(boxes(creator).map((b) => [b.checked, b.disabled]), [[true, true], [true, true], [true, true]]);
  assert.equal(creator.querySelector('.u2-access-grid-remove'), null, 'nothing to remove on an inherited row');
  assert.deepEqual(boxes(sales).map((b) => [b.checked, b.disabled]), [[true, false], [false, false], [false, false]]);
  assert.notEqual(sales.querySelector('.u2-access-grid-remove'), null);
  g.dispose();
});

ui('a checkbox change rewrites the row; remove drops it; the value never carries inherited rows', () => {
  const g = grid();
  const changes = [];
  g.effect(() => changes.push(g.value.value));
  const sales = rows(g)[1];
  boxes(sales)[1].click();
  assert.deepEqual(g.value.value, [{principal: 'Sales', can: {view: true, edit: true, delete: false}}]);
  fire(rows(g)[1].querySelector('.u2-access-grid-remove'), 'click');
  assert.deepEqual(g.value.value, []);
  assert.equal(rows(g).length, 1, 'the inherited row stays');
  assert.equal(changes.length, 3);
  g.dispose();
});

ui('the picker adds a principal with the default grant and offers only those not yet in a row', () => {
  const g = grid();
  const picker = g.root.querySelector('.u2-access-grid-add');
  assert.deepEqual(picker.childNodes.map((o) => o.value), ['', 'Developers', 'Chemists'],
    'Sales and the creator are taken');
  picker.value = 'Developers';
  fire(picker, 'change');
  assert.deepEqual(g.value.value, [
    {principal: 'Sales', can: {view: true, edit: false, delete: false}},
    {principal: 'Developers', can: {view: true, edit: false, delete: false}},
  ]);
  assert.deepEqual(g.root.querySelector('.u2-access-grid-add').childNodes.map((o) => o.value), ['', 'Chemists']);
  g.dispose();
});

ui('a locked capability keeps its box, disabled, and follows its signal', () => {
  const locked = signal(['edit', 'delete']);
  const g = grid({locked});
  assert.deepEqual(boxes(rows(g)[1]).map((b) => b.disabled), [false, true, true]);
  locked.value = [];
  assert.deepEqual(boxes(rows(g)[1]).map((b) => b.disabled), [false, false, false]);
  g.dispose();
});

ui('registered as u2-access-grid with the value prop two-way', () => {
  const reg = new Registry();
  registerAll(reg);
  const meta = reg.get('u2-access-grid');
  assert.equal(meta.category, 'Inputs');
  assert.equal(meta.props.find((p) => p.name === 'value').twoWay, true);
  const built = meta.create({label: 'Access', capabilities: ['view', 'edit'], principals: ['Sales'],
    value: [{principal: 'Sales', can: {view: true}}]});
  assert.deepEqual(built.root.querySelectorAll('th').map((th) => th.textContent), ['Group', 'View', 'Edit', '']);
  built.dispose();
});
