/* AccessGrid: inherited rows read-only above the editable ones, a checkbox change rewrites the
   row, remove drops it, the picker adds a principal with the first unlocked grant and leaves out the
   ones already there, a locked capability cannot be granted — not after an enable either (P98); the
   add select adds on a pick or Enter, not on an arrow key (P99); principals follow (P103); rows are
   kept by principal and addressable as items; registered as u2-access-grid. */

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

ui('a locked capability keeps its box, disabled and unchecked whatever the row holds, and follows its signal', () => {
  const locked = signal(['edit', 'delete']);
  const g = grid({locked, value: [{principal: 'Sales', can: {view: true, edit: true, delete: false}}]});
  assert.deepEqual(boxes(rows(g)[1]).map((b) => [b.checked, b.disabled]), [[true, false], [false, true], [false, true]]);
  assert.deepEqual(boxes(rows(g)[0]).map((b) => b.checked), [true, false, false], 'the inherited row too');
  assert.equal(g.value.value[0].can.edit, true, 'the value keeps what was granted');
  locked.value = [];
  assert.deepEqual(boxes(rows(g)[1]).map((b) => [b.checked, b.disabled]), [[true, false], [true, false], [false, false]]);
  g.dispose();
});

ui('the focused cell is kept and stays focused across a change; a removed row hands focus on', () => {
  const g = grid({value: [{principal: 'Sales', can: {view: true, edit: false, delete: false}},
    {principal: 'Developers', can: {view: true, edit: false, delete: false}}]});
  const edit = boxes(rows(g)[2])[1];
  const row = rows(g)[2];
  edit.focus();
  edit.click();
  assert.equal(g.value.value[1].can.edit, true);
  assert.equal(rows(g)[2] === row && boxes(row)[1] === edit, true, 'the row was patched, not rebuilt');
  assert.equal(edit.checked, true);
  assert.equal(document.activeElement === edit, true, 'Developers › Edit keeps the focus');
  const remove = rows(g)[1].querySelector('.u2-access-grid-remove');
  remove.focus();
  fire(remove, 'click');
  assert.deepEqual(g.value.value.map((r) => r.principal), ['Developers']);
  assert.equal(document.activeElement === rows(g)[1].querySelector('.u2-access-grid-remove'), true,
    'the removed row hands focus to the next row\'s remove');
  const last = rows(g)[1].querySelector('.u2-access-grid-remove');
  last.focus();
  fire(last, 'click');
  assert.deepEqual(g.value.value, []);
  assert.equal(document.activeElement === g.root.querySelector('.u2-access-grid-add'), true,
    'no editable row left (the creator row cannot take focus): the add picker');
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

const P98 = {inherited: [{principal: 'You (creator)', can: {view: true, edit: true, delete: true}, from: ''}],
  locked: ['edit', 'delete']};

ui('P98: built disabled then enabled, inherited rows and locked capabilities stay disabled', () => {
  const enabled = signal(false);
  const g = grid({...P98, enabled});
  assert.deepEqual(boxes(rows(g)[1]).map((b) => b.disabled), [true, true, true]);
  enabled.value = true;
  assert.deepEqual(boxes(rows(g)[0]).map((b) => b.disabled), [true, true, true], 'the inherited row');
  assert.deepEqual(boxes(rows(g)[1]).map((b) => b.disabled), [false, true, true], 'the locked ones');
  assert.equal(rows(g)[1].querySelector('.u2-access-grid-remove').disabled, false);
  g.enabled = true;
  assert.deepEqual(boxes(rows(g)[1]).map((b) => b.disabled), [false, true, true], 'a repeated enable sweep too');
  g.enabled = false;
  g.value.value = [...g.value.peek(), {principal: 'Developers', can: {view: true}}];
  assert.deepEqual(boxes(rows(g)[2]).map((b) => b.disabled), [true, true, true], 'a row added while disabled');
  g.dispose();
});

ui('P98: a click on a locked box, an inherited box or a disabled grid changes nothing', () => {
  const g = grid(P98);
  const edit = boxes(rows(g)[1])[1];
  edit.disabled = false;
  edit.click();
  assert.equal(edit.checked, false, 'springs back');
  assert.deepEqual(g.value.value, [{principal: 'Sales', can: {view: true, edit: false, delete: false}}]);
  const inherited = boxes(rows(g)[0])[0];
  inherited.disabled = false;
  inherited.click();
  assert.equal(inherited.checked, true);
  g.enabled = false;
  const view = boxes(rows(g)[1])[0];
  view.disabled = false;
  view.click();
  assert.equal(view.checked, true);
  fire(rows(g)[1].querySelector('.u2-access-grid-remove'), 'click');
  assert.deepEqual(g.value.value.map((r) => r.principal), ['Sales'], 'nor removes a row');
  g.dispose();
});

ui('P98: an added principal is granted the first capability not locked', () => {
  const g = grid({locked: ['view']});
  const picker = g.root.querySelector('.u2-access-grid-add');
  picker.value = 'Developers';
  fire(picker, 'change');
  assert.deepEqual(g.value.value[1], {principal: 'Developers', can: {view: false, edit: true, delete: false}});
  assert.deepEqual(boxes(rows(g)[2]).map((b) => b.checked), [false, true, false], 'what the boxes show');
  g.dispose();
});

ui('P99: a key that moves the add select\'s choice only moves it; Enter adds; leaving it resets; a later pick adds', async () => {
  const g = grid();
  const picker = g.root.querySelector('.u2-access-grid-add');
  picker.focus();
  // the browser moves a closed select's choice on ArrowDown and fires change at once
  fire(picker, 'keydown', {key: 'ArrowDown'});
  picker.value = 'Developers';
  fire(picker, 'change');
  assert.equal(g.value.value.length, 1, 'nothing added by the arrow');
  assert.equal(fire(picker, 'keydown', {key: 'Enter'}), false, 'Enter is the select\'s');
  assert.deepEqual(g.value.value.map((r) => r.principal), ['Sales', 'Developers']);
  fire(picker, 'keydown', {key: 'c'});
  picker.value = 'Chemists';
  fire(picker, 'change');
  assert.equal(g.value.value.length, 2, 'nor by a letter');
  picker.blur();
  assert.equal(picker.value, '', 'back to the placeholder');
  picker.focus();
  fire(picker, 'keydown', {key: 'ArrowDown'});
  await flush();
  picker.value = 'Chemists';
  fire(picker, 'change');
  assert.deepEqual(g.value.value.map((r) => r.principal), ['Sales', 'Developers', 'Chemists'],
    'a pick from the list the key opened (macOS), in a later task, adds');
  g.dispose();
});

ui('P103: setPrincipals and a bound principals list feed the picker and the row labels', () => {
  const g = grid();
  const offered = () => g.root.querySelector('.u2-access-grid-add').childNodes.map((o) => o.value);
  g.setPrincipals([{value: 'Sales', label: 'Sales team'}, 'Biologists']);
  assert.deepEqual(offered(), ['', 'Biologists']);
  assert.equal(rows(g)[1].querySelector('.u2-access-grid-principal').textContent, 'Sales team');
  assert.deepEqual(g.principals.map((p) => typeof p === 'string' ? p : p.value), ['Sales', 'Biologists']);
  g.dispose();
  const reg = new Registry();
  registerAll(reg);
  const principals = signal(['Sales']);
  const built = reg.get('u2-access-grid').create({label: 'Access', capabilities: ['view'], principals, value: []});
  principals.value = ['Sales', 'Chemists'];
  assert.deepEqual(built.root.querySelector('.u2-access-grid-add').childNodes.map((o) => o.value), ['', 'Sales', 'Chemists']);
  built.dispose();
});

ui('rows are items keyed by principal: data-u2-item, role row, item-label; inherited ones disabled', async () => {
  const g = grid();
  assert.deepEqual(rows(g).map((r) => [r.dataset.u2Item, r.getAttribute('role'),
    r.querySelector('[data-u2-part="item-label"]').textContent]),
  [['inherited:You (creator)', 'row', 'You (creator)'], ['Sales', 'row', 'Sales']]);
  const state = g.getState();
  assert.deepEqual(state.items.list.map((i) => [i.key, i.label, i.disabled]),
    [['inherited:You (creator)', 'You (creator)', true], ['Sales', 'Sales', undefined]]);
  assert.equal(state.readOnly, false);
  assert.equal(state.text, 'Sales', 'the editable rows, as shown');
  assert.equal(state.label, 'Access');
  assert.equal(g.root.querySelector('table').getAttribute('role'), 'grid');
  assert.equal(await g.revealItem('sales') === rows(g)[1], true);
  g.value.value = null;
  assert.deepEqual(rows(g).map((r) => r.dataset.u2Item), ['inherited:You (creator)'], 'a null value shows no rows');
  g.dispose();
});
