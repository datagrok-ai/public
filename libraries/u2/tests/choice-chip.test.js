import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush} from './dom-shim.js';
import {ChoiceChip, computed, signal} from '../src/index.js';
import {u2test, mount, byKind, state, expectOperable, user} from '../src/testing/index.js';

const CAPS = [{value: 'H', label: 'H', description: 'hydrogen · [*:1][H]'},
  {value: 'OH', label: 'OH', description: 'hydroxyl · O[*:1]'}, {value: 'Br', label: 'Br', description: 'bromo · Br[*:1]'}];

const menu = () => document.querySelector('[data-u2="menu"]');
const items = () => [...(menu()?.querySelectorAll('[role="menuitemcheckbox"], [role="menuitem"]') ?? [])].map((el) => ({el,
  label: el.querySelector('[data-u2-part="item-label"]').textContent,
  detail: el.querySelector('[data-u2-part="item-detail"]')?.textContent ?? '', checked: el.getAttribute('aria-checked')}));
const text = (el) => el.textContent.replace(/\s+/g, ' ').trim();

test('choiceChip: reads "label value", a button that opens a menu, and says its state', u2test(() => {
  const chip = mount(() => new ChoiceChip({label: 'R1', name: 'r1', items: CAPS, value: 'OH'}));
  assert.equal(chip.root.tagName, 'BUTTON');
  assert.equal(text(chip.root), 'R1 OH');
  assert.equal(chip.root.getAttribute('aria-haspopup'), 'menu');
  assert.equal(chip.root.getAttribute('aria-expanded'), 'false');
  assert.equal(chip.root.getAttribute('aria-description'), 'hydroxyl · O[*:1]', 'the current choice says what it is');
  assert.equal(chip.root.getAttribute('aria-invalid'), 'false');
  assert.deepEqual(byKind(document, 'choice-chip'), [chip.root]);
  const s = state(chip.root);
  assert.equal(s.kind, 'choice-chip');
  assert.equal(s.name, 'r1');
  assert.equal(s.label, 'R1');
  assert.equal(s.value, 'OH');
  assert.equal(s.text, 'OH');
  assert.equal(s.valid, true);
  assert.equal(s.open, false);
  assert.deepEqual(s.selection, ['OH']);
  expectOperable(chip.root);
}));

test('choiceChip: a click opens the choices, each with what it stands for, the current one checked; a pick is the ' +
  'user\'s input and closes the menu', u2test(() => {
  const picked = [];
  const chip = mount(() => new ChoiceChip({label: 'R1', items: CAPS, value: 'OH', onInput: (v) => picked.push(v)}));
  user.click(chip);
  assert.equal(chip.root.getAttribute('aria-expanded'), 'true');
  assert.equal(chip.isOpen, true);
  assert.deepEqual(items().map((i) => [i.label, i.detail, i.checked]),
    [['H', 'hydrogen · [*:1][H]', 'false'], ['OH', 'hydroxyl · O[*:1]', 'true'], ['Br', 'bromo · Br[*:1]', 'false']]);
  assert.equal(menu().dataset.u2OwnerId, chip.root.id, 'the menu is the chip\'s popup');
  expectOperable(menu());
  user.pick(chip, 'Br');
  assert.equal(chip.value, 'Br');
  assert.deepEqual(picked, ['Br']);
  assert.equal(menu(), null, 'a pick closes the menu');
  assert.equal(chip.root.getAttribute('aria-expanded'), 'false');
  assert.equal(text(chip.root), 'R1 Br');
  chip.value = 'H';
  assert.deepEqual(picked, ['Br'], 'a write from code is no input');
  assert.equal(text(chip.root), 'R1 H');
  user.click(chip);
  user.click(chip);
  assert.equal(menu(), null, 'a second press closes the menu');
}));

test('choiceChip: the keyboard opens it with ArrowDown and Escape gives the focus back', u2test(() => {
  const chip = mount(() => new ChoiceChip({label: 'Sort', items: ['Name', 'Created'], value: 'Name'}));
  user.press(chip, 'ArrowDown');
  assert.equal(chip.isOpen, true);
  assert.equal(document.activeElement, menu(), 'the menu takes the keys');
  user.press(menu(), 'Escape');
  assert.equal(chip.isOpen, false);
  assert.equal(document.activeElement, chip.root, 'the focus is back on the chip');
  user.press(chip, 'Enter');
  assert.equal(chip.isOpen, true, 'Enter presses the button');
  chip.dispose();
  assert.equal(menu(), null, 'the menu goes with the chip');
}));

test('choiceChip: invalid shows at once, with its reason, and follows a signal', u2test(() => {
  const problem = signal('X is a placeholder, not a cap group');
  const chip = mount(() => new ChoiceChip({label: 'R4', items: CAPS, valueText: 'X', invalid: problem}));
  assert.equal(chip.root.getAttribute('aria-invalid'), 'true', 'no touch needed: the chip shows its problem');
  assert.equal(chip.root.getAttribute('aria-description'), 'X is a placeholder, not a cap group');
  assert.equal(state(chip.root).valid, false);
  assert.equal(state(chip.root).error, 'X is a placeholder, not a cap group');
  assert.equal(text(chip.root), 'R4 X');
  problem.value = null;
  assert.equal(chip.root.getAttribute('aria-invalid'), 'false');
  assert.equal(chip.root.hasAttribute('aria-description'), false, 'a value none of the items describe says nothing');
}));

test('choiceChip: off with a reason refuses the press and says why', u2test(() => {
  const off = signal('Draw R groups first');
  const chip = mount(() => new ChoiceChip({label: 'R1', items: CAPS, value: 'H', disabled: off}));
  assert.throws(() => user.click(chip), {message: 'user.click: button "R1 H" is disabled: Draw R groups first'});
  assert.equal(chip.isOpen, false);
  assert.equal(state(chip.root).disabledReason, 'Draw R groups first');
  assert.equal(chip.disabled, 'Draw R groups first');
  expectOperable(chip.root);
  off.value = null;
  user.click(chip);
  assert.equal(chip.isOpen, true);
  assert.equal(chip.root.getAttribute('aria-description'), 'hydrogen · [*:1][H]', 'on again, it describes its value');
  off.value = 'Locked';
  assert.equal(chip.isOpen, false, 'turning off closes its menu');
}));

test('choiceChip: what it shows and says can follow a model; actions follow the choices; active looks open',
  u2test(() => {
    const def = signal({cap: 'OH', note: 'read from the alternate ID'});
    const editing = signal(false);
    const chip = mount(() => new ChoiceChip({label: 'R2', items: CAPS,
      valueText: computed(() => def.value.cap), description: computed(() => `hydroxyl (${def.value.note})`),
      actions: () => [{text: editing.peek() ? 'Hide definition' : 'Edit definition…',
        run: () => { editing.value = !editing.peek(); }}],
      active: editing}));
    chip.follow(() => def.value.cap);
    assert.equal(chip.value, 'OH');
    assert.equal(chip.root.getAttribute('aria-description'), 'hydroxyl (read from the alternate ID)');
    def.value = {cap: 'Azide', note: 'as defined'};
    assert.equal(chip.value, 'Azide', 'follow keeps the value on the model');
    assert.equal(text(chip.root), 'R2 Azide');
    user.click(chip);
    assert.deepEqual(items().map((i) => i.label), ['H', 'OH', 'Br', 'Edit definition…']);
    assert.equal(items().filter((i) => i.checked === 'true').length, 0, 'a value none of the items is checks nothing');
    assert.equal(menu().querySelectorAll('[role="separator"]').length, 1);
    user.pick(chip, 'Edit definition…');
    assert.equal(editing.value, true);
    assert.equal(chip.root.classList.contains('u2-choice-chip-active'), true);
    user.click(chip);
    assert.deepEqual(items().map((i) => i.label).slice(-1), ['Hide definition'], 'asked again on every open');
    user.pick(chip, 'Hide definition');
    assert.equal(chip.root.classList.contains('u2-choice-chip-active'), false);
  }));

test('choiceChip: items given as a signal are read on every open', u2test(async () => {
  const caps = signal(CAPS.slice(0, 2));
  const chip = mount(() => new ChoiceChip({label: 'R1', items: caps, value: 'Br'}));
  assert.equal(text(chip.root), 'R1 –', 'a value none of the items shows the dash');
  caps.value = CAPS;
  assert.equal(text(chip.root), 'R1 Br', 'and its label once its item arrives');
  user.click(chip);
  assert.deepEqual(items().map((i) => i.label), ['H', 'OH', 'Br']);
  assert.deepEqual(chip.getItems().list.map((i) => [i.key, i.selected]), [['H', false], ['OH', false], ['Br', true]]);
  user.click(chip);
  await flush();
}));
