import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire} from './dom-shim.js';
import {layout} from './layout.mjs';
import {operability} from 'datagrok-api/u2core';
import {TextInput, DynamicInput, ColorInput, Input, Reason, button, iconButton, div} from '../src/index.js';
import {SliderInput} from '../src/components/inputs/slider-input.js';
import {u2test, mount, expectOperable} from '../src/testing/index.js';

const NO_CASCADE = typeof CSSStyleSheet !== 'function' && 'needs a DOM with the CSS cascade';

const findings = (root, options) => operability(root, options).map((f) => [f.rule, f.detail]);

function campaignName() {
  const input = mount(() => new TextInput({label: 'Campaign Name', value: 'KHL-3', readOnly: true}));
  const edit = iconButton('pencil', () => {}, {tooltip: 'Edit Campaign Name'});
  input.addOptions(edit);
  return {input, edit};
}

test('a read-only input\'s own option stays operable under today\'s rail rule', u2test(() => {
  const {input} = campaignName();
  assert.deepEqual(findings(input.root), []);
  expectOperable(input);
}));

test('no-pointer: the old read-only rail rule leaves the pencil dead, and the probe names who took the pointer',
  {skip: NO_CASCADE}, u2test(() => {
    const {input, edit} = campaignName();
    const rule = document.createElement('style');
    rule.textContent = '.u2-input-readonly .u2-input-options { pointer-events: none; }';
    document.head.append(rule);
    try {
      const found = operability(input.root);
      assert.deepEqual(found.map((f) => f.rule), ['no-pointer']);
      assert.equal(found[0].element, edit);
      assert.equal(found[0].detail, 'button "Edit Campaign Name" takes no pointer: pointer-events: none on ' +
        'div.u2-input-options');
      assert.throws(() => expectOperable(input), /no-pointer: button "Edit Campaign Name" takes no pointer/);
    }
    finally {
      rule.remove();
    }
  }));

test('an option that edits is off on a read-only input, saying why, and on again once it is editable',
  {skip: NO_CASCADE}, u2test(() => {
    const {input} = campaignName();
    const step = iconButton('plus', () => {}, {tooltip: 'Add one'});
    input.addOptions(step, {edits: true});
    assert.deepEqual(findings(input.root), []);
    assert.equal(Reason.of(step), Input.READ_ONLY);
    expectOperable(input);
    input.readOnly = false;
    assert.equal(Reason.of(step), null);
    const own = iconButton('download', () => {}, {tooltip: 'Download', disabled: 'Nothing to download yet'});
    input.addOptions(own, {edits: true});
    input.readOnly = true;
    input.readOnly = false;
    assert.equal(Reason.of(own), 'Nothing to download yet', 'an option off for its own reason keeps it');
  }));

test('a read-only slider takes the pointer, so its value shows on hover; a drag or a key does not move it',
  {skip: NO_CASCADE}, u2test(() => {
    const slider = mount(() => new SliderInput({label: 'Dose', min: 0, max: 100, value: 25, readOnly: true}));
    const track = slider.root.querySelector('input[type="range"]');
    assert.deepEqual(findings(slider.root), []);
    assert.notEqual(getComputedStyle(track).pointerEvents, 'none');
    assert.equal(fire(track, 'mousedown'), false, 'the press that starts a drag is cancelled');
    assert.equal(fire(track, 'keydown', {key: 'ArrowRight'}), false, 'a stepping key is cancelled');
    track.value = '60';
    fire(track, 'input');
    assert.equal(slider.peek(), 25, 'the value stays');
  }));

test('a read-only color input: its swatch is off, saying why, not a dead button', {skip: NO_CASCADE}, u2test(() => {
  const color = mount(() => new ColorInput({label: 'Color', value: '#ff0000', readOnly: true}));
  assert.deepEqual(findings(color.root), []);
  assert.equal(Reason.of(color.root.querySelector('button')), Input.READ_ONLY);
  color.readOnly = false;
  assert.deepEqual(findings(color.root), []);
  assert.equal(Reason.of(color.root.querySelector('button')), null);
}));

test('covered: with hit testing, an element whose centre the page gives to another is covered', u2test(() => {
  const {input, edit} = campaignName();
  layout(edit, {rect: new DOMRect(100, 10, 16, 28)});
  const row = input.root.querySelector('.u2-input-box');
  const elementFromPoint = document.elementFromPoint;
  document.elementFromPoint = (x, y) => x === 108 && y === 24 ? row : null;
  try {
    assert.deepEqual(findings(input.root), [], 'no hit test without layout');
    assert.deepEqual(findings(input.root, {hitTest: true}),
      [['covered', 'button "Edit Campaign Name" is covered by div.u2-input-box']]);
    document.elementFromPoint = () => edit.firstChild;
    assert.deepEqual(findings(input.root, {hitTest: true}), [], 'a hit inside the element is the element');
  } finally {
    document.elementFromPoint = elementFromPoint;
  }
}));

test('oversized-option: a 160×100 drawing on the 28px rail is flagged, naming DynamicInput', u2test(() => {
  const input = new TextInput({label: 'Reference compound', readOnly: true});
  const drawing = div([], 'molecule');
  layout(drawing, {rect: new DOMRect(790, 593, 160, 100)});
  input.addOptions(drawing);
  layout(input.root.querySelector('[data-u2-part="options"]'), {rect: new DOMRect(790, 593, 160, 28)});
  assert.deepEqual(findings(input.root), [['oversized-option', 'div.molecule is 100px tall on a 28px options rail: ' +
    'an element shown as the value belongs in a DynamicInput row']]);
  const shown = new DynamicInput({label: 'Reference compound', value: 'c1ccccc1',
    render: () => div([], 'molecule')});
  assert.deepEqual(findings(shown.root), [], 'the value in a DynamicInput row is no rail option');
}));

test('inert: a region made inert with no reason is reported once, with a reason not at all', u2test(() => {
  const pane = mount(() => div([button('Run', () => {}), button('Edit', () => {})], 'pane'));
  pane.setAttribute('inert', '');
  assert.deepEqual(findings(pane), [['inert', 'div.pane is inert and says no reason: button "Run" does nothing']]);
  pane.title = 'Check Compute to change its settings';
  assert.deepEqual(findings(pane), []);
}));

test('disabled-without-reason: a reasonless off button and input are reported; a reasoned one is not', u2test(() => {
  const off = Object.assign(button('Save', () => {}), {disabled: true});
  const why = button('Submit', () => {}, {disabled: 'Name the campaign first'});
  const input = new TextInput({label: 'Name'});
  // the API takes no off without a reason: one taken off its root by hand is the misuse left
  input.disabled = 'Name the campaign first';
  input.root.removeAttribute('aria-description');
  const box = mount(() => div([off, why, input.root]));
  const found = operability(box);
  assert.deepEqual(found.map((f) => [f.rule, f.element]),
    [['disabled-without-reason', off], ['disabled-without-reason', input.root]],
    'one finding for the input, not one per control inside it');
  assert.equal(found[0].detail, 'button "Save" is disabled and says no reason');
}));

test('name-mismatch: an icon button whose accessible name is not its tooltip, or that has no name', u2test(() => {
  const edit = iconButton('pencil', () => {}, {tooltip: 'Edit'});
  const box = mount(() => div([edit]));
  assert.deepEqual(findings(box), []);
  edit.setAttribute('aria-label', 'Rename');
  assert.deepEqual(findings(box), [['name-mismatch',
    'button.u2-icon-btn is "Rename" to a screen reader and "Edit" on hover']]);
  edit.removeAttribute('aria-label');
  assert.deepEqual(findings(box), [['name-mismatch',
    'button.u2-icon-btn has no accessible name; its tooltip says "Edit"']]);
}));

test('what a user cannot see is not probed', u2test(() => {
  const off = Object.assign(button('Save', () => {}), {disabled: true});
  const box = div([div([off])]);
  box.firstChild.hidden = true;
  assert.deepEqual(findings(box), []);
}));
