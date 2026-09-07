/* MultiChoiceInput summary row (`showSummaryCheckbox`): a tri-state "<N> of <M>" checkbox that
   toggles all/none, with the item list collapsing on the text/chevron — like the platform's
   func-param tree group (func_param_editor.dart:754-830), but with the summary checkbox
   pixel-flush over the item checkbox column. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {MultiChoiceInput} from '../src/components/inputs/choice-input.js';

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

const ITEMS = ['Liver', 'Kidney', 'Heart', 'Lung', 'Brain'];

function summary(input) {
  return input.root.querySelector('.u2-multi-choice-summary');
}

function summaryBox(input) {
  return summary(input).querySelector('.u2-input-checkbox');
}

function summaryText(input) {
  return input.root.querySelector('.u2-multi-choice-summary-text').textContent;
}

function toggle(input) {
  return input.root.querySelector('.u2-multi-choice-summary-toggle');
}

function itemsHost(input) {
  return input.root.querySelector('.u2-multi-choice-items');
}

function itemBoxes(input) {
  return input.root.querySelectorAll('.u2-multi-choice-item .u2-input-checkbox');
}

function checkItem(input, index, checked) {
  const box = itemBoxes(input)[index];
  box.checked = checked;
  fire(box, 'change');
}

smoke('the summary reads "N of M" on init and follows checks and unchecks', () => {
  const input = mount(new MultiChoiceInput({label: 'Tissues', items: ITEMS,
    value: ['Liver', 'Kidney'], showSummaryCheckbox: true}));
  assert.equal(summaryText(input), '2 of 5');

  checkItem(input, 2, true);
  assert.equal(summaryText(input), '3 of 5');
  assert.deepEqual(input.value.value, ['Liver', 'Kidney', 'Heart']);

  checkItem(input, 0, false);
  checkItem(input, 1, false);
  checkItem(input, 2, false);
  assert.equal(summaryText(input), '0 of 5');

  input.value.value = ITEMS.slice();
  assert.equal(summaryText(input), '5 of 5', 'programmatic writes reach the summary too');
  input.dispose();
});

smoke('the summary checkbox is tri-state: unchecked, indeterminate, checked', () => {
  const input = mount(new MultiChoiceInput({items: ITEMS, showSummaryCheckbox: true}));
  assert.equal(summaryBox(input).checked, false);
  assert.equal(summaryBox(input).indeterminate, false);

  checkItem(input, 0, true);
  assert.equal(summaryBox(input).checked, false);
  assert.equal(summaryBox(input).indeterminate, true, 'some but not all');

  input.value.value = ITEMS.slice();
  assert.equal(summaryBox(input).checked, true);
  assert.equal(summaryBox(input).indeterminate, false);

  input.value.value = [];
  assert.equal(summaryBox(input).checked, false);
  assert.equal(summaryBox(input).indeterminate, false);
  input.dispose();
});

smoke('clicking the summary checkbox selects all from indeterminate, then none, as user edits', () => {
  const changes = [];
  const input = mount(new MultiChoiceInput({items: ITEMS, value: ['Heart'],
    showSummaryCheckbox: true, onChanged: (v) => changes.push(v)}));
  assert.equal(summaryBox(input).indeterminate, true);

  // a click on an indeterminate checkbox lands as checked
  summaryBox(input).checked = true;
  fire(summaryBox(input), 'change');
  assert.deepEqual(input.value.value, ITEMS, 'all selected, in item order');
  assert.equal(summaryText(input), '5 of 5');
  assert.equal(summaryBox(input).checked, true);
  for (const box of itemBoxes(input))
    assert.equal(box.checked, true);

  summaryBox(input).checked = false;
  fire(summaryBox(input), 'change');
  assert.deepEqual(input.value.value, []);
  assert.equal(summaryText(input), '0 of 5');
  for (const box of itemBoxes(input))
    assert.equal(box.checked, false);

  assert.deepEqual(changes, [ITEMS, []], 'both writes fired onChanged, as ordinary user edits');
  input.dispose();
});

smoke('the text and the chevron collapse the list; the content is hidden, never detached', () => {
  const input = mount(new MultiChoiceInput({items: ITEMS, value: ['Liver'],
    showSummaryCheckbox: true}));
  assert.equal(input.expanded.value, true);
  assert.equal(toggle(input).getAttribute('aria-expanded'), 'true');
  assert.notEqual(itemsHost(input).style.display, 'none');

  fire(toggle(input), 'click');
  assert.equal(input.expanded.value, false);
  assert.equal(toggle(input).getAttribute('aria-expanded'), 'false');
  assert.equal(itemsHost(input).style.display, 'none');
  assert.equal(itemBoxes(input).length, 5, 'rows stay in the DOM');
  assert.deepEqual(input.value.value, ['Liver'], 'the value survives collapse');

  fire(toggle(input), 'click');
  assert.equal(input.expanded.value, true);
  assert.notEqual(itemsHost(input).style.display, 'none');
  assert.equal(itemBoxes(input)[0].checked, true, 'checked state restored intact');
  input.dispose();
});

smoke('Enter and Space on the toggle collapse and expand', () => {
  const input = mount(new MultiChoiceInput({items: ITEMS, showSummaryCheckbox: true}));
  fire(toggle(input), 'keydown', {key: 'Enter'});
  assert.equal(input.expanded.value, false);
  fire(toggle(input), 'keydown', {key: ' '});
  assert.equal(input.expanded.value, true);
  fire(toggle(input), 'keydown', {key: 'Escape'});
  assert.equal(input.expanded.value, true, 'other keys are ignored');
  input.dispose();
});

smoke('setItems recomputes the text and preserves the collapse state', () => {
  const input = mount(new MultiChoiceInput({items: ITEMS, value: ['Liver', 'Heart'],
    showSummaryCheckbox: true}));
  fire(toggle(input), 'click');
  assert.equal(itemsHost(input).style.display, 'none');

  input.setItems(['Liver', 'Kidney', 'Spleen']);
  assert.equal(summaryText(input), '1 of 3', 'Heart vanished with its check');
  assert.deepEqual(input.value.value, ['Liver']);
  assert.equal(itemsHost(input).style.display, 'none', 'still collapsed');
  assert.equal(input.expanded.value, false);

  fire(toggle(input), 'click');
  assert.equal(itemBoxes(input).length, 3);
  input.dispose();
});

smoke('empty items read "0 of 0" and the summary checkbox is a no-op', () => {
  const changes = [];
  const input = mount(new MultiChoiceInput({items: [], showSummaryCheckbox: true,
    onChanged: (v) => changes.push(v)}));
  assert.equal(summaryText(input), '0 of 0');
  assert.notEqual(input.root.querySelector('.u2-multi-choice-empty'), null);

  summaryBox(input).checked = true;
  fire(summaryBox(input), 'change');
  assert.deepEqual(input.value.value, []);
  assert.equal(summaryBox(input).checked, false, 'snapped back to unchecked');
  assert.deepEqual(changes, [], 'no write, no onChanged');
  input.dispose();
});

smoke('without the option the editor is exactly the plain checkbox list', () => {
  const input = mount(new MultiChoiceInput({label: 'Tissues', items: ITEMS, value: ['Liver']}));
  assert.equal(summary(input), null);
  assert.equal(itemsHost(input), null, 'no wrapper either');
  assert.equal(input.expanded, undefined);
  const editor = input.root.querySelector('.u2-multi-choice');
  for (const child of editor.children)
    assert.equal(child.className, 'u2-multi-choice-item', 'item rows sit directly in the editor');
  assert.equal(editor.children.length, 5);

  checkItem(input, 1, true);
  assert.deepEqual(input.value.value, ['Liver', 'Kidney'], 'value semantics unchanged');
  input.dispose();
});

smoke('the summary checkbox is named by the "N of M" text, and disabling parks the toggle', () => {
  const input = mount(new MultiChoiceInput({items: ITEMS, value: ['Liver'],
    showSummaryCheckbox: true}));
  const text = input.root.querySelector('.u2-multi-choice-summary-text');
  assert.notEqual(text.id, '', 'the text span carries an id');
  assert.equal(summaryBox(input).getAttribute('aria-labelledby'), text.id);
  assert.equal(toggle(input).getAttribute('aria-label'), 'Collapse list',
    'the toggle never shares the checkbox\'s "N of M" name');
  assert.equal(toggle(input).tabIndex, 0);

  input.enabled = false;
  assert.equal(toggle(input).tabIndex, -1, 'the span gives up its tab stop with the input');
  fire(toggle(input), 'keydown', {key: 'Enter'});
  fire(toggle(input), 'click');
  assert.equal(input.expanded.value, true, 'gestures are inert while disabled');

  input.enabled = true;
  assert.equal(toggle(input).tabIndex, 0);
  fire(toggle(input), 'click');
  assert.equal(input.expanded.value, false, 'and live again');
  input.dispose();
});

smoke('the summary checkbox shares the item checkboxes\' class and column structure', () => {
  const input = mount(new MultiChoiceInput({items: ITEMS, value: ['Liver'],
    showSummaryCheckbox: true}));
  const editor = input.root.querySelector('.u2-multi-choice');
  const sBox = summaryBox(input);
  const iBox = itemBoxes(input)[0];

  assert.ok(sBox.classList.contains('u2-input-checkbox'), 'same checkbox class, same margin rule');
  assert.equal(summary(input).parentNode, editor, 'summary row hangs off the editor column');
  assert.equal(itemsHost(input).parentNode, editor, 'and so does the items wrapper — same left edge');
  assert.equal(summary(input).children[0], sBox,
    'the checkbox is the summary row\'s first element, nothing indents it');
  assert.equal(iBox.parentNode.children[0], iBox,
    'exactly as in an item row — the Dart -20px offset is deliberately not copied');
  input.dispose();
});
