/* BoolInput with a postfix: the box comes first in the input box and the hint after it, on the
   options rail; and the skin keeps the editor at its own width — a wide form gives every editor
   `flex: 1`, and a bool editor shrunk to nothing puts its box under the hint, where a click
   lands on the hint. Layout is not measurable here, so the rules are asserted as written. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {resetDom} from './dom-shim.js';
import {BoolInput} from '../src/components/inputs/bool-input.js';
import {signal} from '../src/core/signals.js';
import {Form} from '../src/components/forms/form.js';

const css = (name) => readFileSync(new URL(`../css/${name}`, import.meta.url), 'utf8');

/** The declarations of the rule whose selector list contains `selector`, verbatim. */
function rule(sheet, selector) {
  const at = sheet.indexOf(selector);
  assert.notEqual(at, -1, `${selector} is a rule`);
  return sheet.slice(sheet.indexOf('{', at) + 1, sheet.indexOf('}', at));
}

test('the box precedes its postfix inside the input box; the hint is on the options rail', () => {
  try {
    const form = new Form({layout: 'wide'});
    const input = new BoolInput({label: 'Writable', name: 'writable', postfix: 'users with Edit may write'});
    form.add(input);
    document.body.append(form.root);
    const box = input.root.querySelector('.u2-input-box');
    assert.deepEqual([...box.children].map((c) => c.className), ['u2-input-bool u2-input-editor', 'u2-div u2-input-options']);
    assert.equal(box.querySelector('.u2-input-editor').querySelector('.u2-input-checkbox').type, 'checkbox');
    assert.equal(box.querySelector('.u2-input-options').querySelector('.u2-input-postfix').textContent,
      'users with Edit may write');
    form.dispose();
  } finally {
    resetDom();
  }
});

test('the bool editor never shrinks: the skin says so, at a specificity the wide form\'s stretch does not reach', () => {
  assert.match(rule(css('inputs.css'), '.u2-input-box > .u2-input-editor.u2-input-bool'), /flex: none;/);
  // the wide form's exception is `flex: none` too: stopping only the growth would leave its
  // `flex: 1` basis of 0 on the radio and check lists, which have no skin rule of their own and
  // collapsed to 0px
  const wide = rule(css('form.css'), '.u2-form-wide .u2-input-root[data-u2=\'bool-input\'] .u2-input-editor');
  assert.match(wide, /flex: none;/);
  assert.doesNotMatch(wide, /flex: 1|flex-grow|flex-shrink|flex-basis/);
});

test('captionAfter: the checklist row — box first, its caption after it, tied by for; no label column', async () => {
  try {
    const form = new Form({layout: 'normal'});
    const row = new BoolInput({label: 'Include salts', name: 'salts', captionAfter: true});
    form.add(row);
    document.body.append(form.root);
    const box = row.root.querySelector('input[type="checkbox"]');
    const caption = row.root.querySelector('[data-u2-part="label"]');
    assert.equal(row.root.classList.contains('u2-input-caption-after'), true);
    assert.equal(caption.getAttribute('for'), box.id, 'a click on the caption toggles the box natively');
    const inputs = css('inputs.css');
    assert.match(rule(inputs, '.u2-input-caption-after > .u2-input-label'), /order: 1;/);
    assert.match(rule(css('form.css'), '.u2-form .u2-input-caption-after > .u2-input-label'), /min-width: 0;/);
    form.dispose();
    row.dispose();
  } finally {
    resetDom();
  }
});

test('a switch toggles from its caption too; read-only or disabled it does not', async () => {
  try {
    const sw = new BoolInput({label: 'Notify', switch: true});
    document.body.append(sw.root);
    const caption = sw.root.querySelector('[data-u2-part="label"]');
    caption.click();
    assert.equal(sw.value.peek(), true);
    sw.readOnly = true;
    caption.click();
    assert.equal(sw.value.peek(), true, 'read-only');
    sw.readOnly = false;
    sw.enabled = false;
    caption.click();
    assert.equal(sw.value.peek(), true, 'disabled');
    sw.dispose();
  } finally {
    resetDom();
  }
});

test('required: false is an answer, so a fresh box is valid; only a null value is empty', () => {
  try {
    const box = new BoolInput({label: 'Agree', required: true});
    assert.equal(box.value.peek(), false);
    assert.equal(box.getState().valid, true);
    const bound = signal(null);
    const unanswered = new BoolInput({label: 'Agree', required: true, bind: bound});
    unanswered.touch();
    assert.equal(unanswered.getState().valid, false);
    bound.value = false;
    assert.equal(unanswered.getState().valid, true);
    box.dispose();
    unanswered.dispose();
  } finally {
    resetDom();
  }
});
