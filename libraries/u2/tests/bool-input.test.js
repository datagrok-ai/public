/* BoolInput with a postfix: the box comes first in the input box and the hint after it, on the
   options rail; and the skin keeps the editor at its own width — a wide form gives every editor
   `flex: 1`, and a bool editor shrunk to nothing puts its box under the hint, where a click
   lands on the hint. Layout is not measurable here, so the rules are asserted as written. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {resetDom} from './dom-shim.js';
import {BoolInput} from '../src/components/inputs/bool-input.js';
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
  // the wide form's exception only stops the growth: the skin's `flex: none` (0,3,0) outranks
  // its `flex: 1` (0,2,0) on shrink and basis, and the (0,4,0) exception leaves both alone
  const wide = rule(css('form.css'), '.u2-form-wide .u2-input-root[data-u2=\'bool-input\'] .u2-input-editor');
  assert.match(wide, /flex-grow: 0;/);
  assert.doesNotMatch(wide, /flex: 1|flex-shrink|flex-basis/);
});
