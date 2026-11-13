import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire} from './dom-shim.js';
import {MasterDetail, TextInput, Tooltip, computed, divV, signal} from '../src/index.js';
import {u2test, mount, byKind, byName, state, expectOperable, user} from '../src/testing/index.js';

const FUNCTIONS = ['Descriptors', 'Solubility', 'Toxicity'].map((name) => ({key: name, name, checked: signal(false)}));

function computeEditor(functions = FUNCTIONS, built = []) {
  return new MasterDetail({items: functions, keyOf: (f) => f.key, name: 'compute',
    row: (f) => ({label: f.name, check: f.checked, hint: computed(() => f.checked.value ? 'on' : '')}),
    detail: (f) => {
      built.push(f.key);
      return divV([new TextInput({label: `${f.name} threshold`, name: `${f.key}-threshold`})]);
    },
    disabled: (f) => computed(() => f.checked.value ? null : `Check ${f.name} to change its settings`)});
}

const rows = (md) => [...md.root.querySelectorAll('[role="option"]')];
const panes = (md) => [...md.root.querySelectorAll('[data-u2-part="detail"] [role="region"]')];
const shown = (md) => panes(md).filter((p) => !p.hidden);
const paneName = (p) => document.getElementById(p.getAttribute('aria-labelledby'))?.textContent ?? null;

test('masterDetail: rows with a check and a hint beside the current one\'s detail, built when first shown and kept',
  u2test(() => {
    const built = [];
    const md = mount(() => computeEditor(FUNCTIONS, built));
    assert.deepEqual(byKind(document, 'master-detail'), [md.root]);
    assert.equal(md.current.value, 'Descriptors', 'the first item is current');
    assert.deepEqual(rows(md).map((r) => [r.dataset.u2Item, r.getAttribute('aria-selected'),
      r.getAttribute('aria-checked')]), [['Descriptors', 'true', 'false'], ['Solubility', 'false', 'false'],
      ['Toxicity', 'false', 'false']]);
    assert.deepEqual(built, ['Descriptors'], 'only the current detail is built');
    user.click(rows(md)[1]);
    assert.equal(md.current.value, 'Solubility');
    assert.deepEqual(built, ['Descriptors', 'Solubility']);
    assert.deepEqual(shown(md).map(paneName), ['Solubility'], 'a detail is named by its row');
    user.click(rows(md)[0]);
    assert.deepEqual(built, ['Descriptors', 'Solubility'], 'a detail shown before is kept, not built again');
    const s = state(md.root);
    assert.equal(s.current, 'Descriptors');
    assert.deepEqual(s.selection, ['Descriptors']);
    assert.deepEqual(s.items.list.map((i) => [i.key, i.label, i.checked]),
      [['Descriptors', 'Descriptors', false], ['Solubility', 'Solubility', false], ['Toxicity', 'Toxicity', false]]);
  }));

test('masterDetail: a detail off says why, out of reach, until its item is checked', u2test(() => {
  const md = mount(() => computeEditor());
  const pane = shown(md)[0];
  const reason = pane.querySelector('[data-u2-part="reason"]');
  assert.equal(reason.textContent, 'Check Descriptors to change its settings');
  assert.equal(reason.hidden, false, 'the reason is on screen, over the dimmed settings');
  assert.throws(() => user.type(byName(md, 'Descriptors-threshold'), '5'),
    {message: /inert: Check Descriptors to change its settings/});
  expectOperable(md.root);
  user.click(rows(md)[0].querySelector('input[type="checkbox"]'));
  assert.equal(FUNCTIONS[0].checked.value, true);
  assert.equal(rows(md)[0].getAttribute('aria-checked'), 'true');
  assert.equal(reason.hidden, true);
  user.type(byName(md, 'Descriptors-threshold'), '5');
  assert.equal(rows(md)[0].querySelector('.u2-master-detail-hint').textContent, 'on', 'the hint follows');
  FUNCTIONS[0].checked.value = false;
}));

test('masterDetail: the arrow keys move the current item, Space checks it; items that leave take their detail',
  u2test(() => {
    const items = signal(FUNCTIONS.map((f) => ({...f, checked: signal(false)})));
    const md = mount(() => new MasterDetail({items, keyOf: (f) => f.key, row: (f) => ({label: f.name, check: f.checked}),
      detail: (f) => divV([new TextInput({label: f.name})])}));
    const list = md.root.querySelector('[data-u2-part="list"]');
    user.press(list, 'ArrowDown');
    assert.equal(md.current.value, 'Solubility');
    assert.equal(list.getAttribute('aria-activedescendant'), rows(md)[1].id);
    user.press(list, 'End');
    assert.equal(md.current.value, 'Toxicity');
    user.press(list, ' ');
    assert.equal(items.value[2].checked.value, true);
    user.press(list, 'Home');
    assert.equal(md.current.value, 'Descriptors');
    assert.equal(panes(md).length, 3, 'each detail the keys reached was built once');
    items.value = items.value.slice(1);
    assert.equal(md.current.value, 'Solubility', 'the current item gone, the first one left is current');
    assert.equal(rows(md).length, 2);
    assert.equal(panes(md).some((p) => paneName(p) === 'Descriptors'), false);
    items.value = [];
    assert.equal(md.current.value, null);
    assert.equal(md.root.querySelector('.u2-master-detail-empty').hidden, false);
  }));

test('masterDetail: the list takes its accessible name, a row\'s check its tooltip', u2test(() => {
  const md = mount(() => new MasterDetail({items: FUNCTIONS, keyOf: (f) => f.key, label: 'Functions',
    row: (f) => ({label: f.name, check: f.checked, checkTooltip: `Toggle calculation of ${f.name}`}),
    detail: (f) => divV([new TextInput({label: f.name})])}));
  assert.equal(md.root.querySelector('[data-u2-part="list"]').getAttribute('aria-label'), 'Functions');
  const box = rows(md)[1].querySelector('input[type="checkbox"]');
  const delay = Tooltip.delay;
  Tooltip.delay = 0;
  try {
    fire(box, 'pointerenter');
    assert.equal(Tooltip.root.textContent, 'Toggle calculation of Solubility');
    fire(box, 'pointerleave');
  }
  finally {
    Tooltip.delay = delay;
  }
  assert.equal(mount(() => computeEditor()).root.querySelector('[data-u2-part="list"]').hasAttribute('aria-label'),
    false, 'no name unless given');
}));
