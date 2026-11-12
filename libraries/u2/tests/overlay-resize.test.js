/* Dialog, overlay, type-ahead and filter-builder fixes from the U2Demo visual review: the platform
   dialog's footer order, OK gated on the body's form, a shown dialog kept on screen through window
   resizes, the error row cleared by an edit, anchored popups closed by a window resize, stray
   type-ahead text dropped on blur, and a new filter row quiet about its missing value until touched. */
import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {layout} from './layout.mjs';
import {Scope, Dialog, Form, TextInput, ChoiceInput, TypeAhead, Overlay, Filters, FilterBuilder, VirtualTree, button}
  from '../src/index.js';
import {readFileSync} from 'node:fs';
import {TYPE} from 'datagrok-api/u2core';

const mounted = [];

function smoke(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      Filters.resetIds('');
      await body();
    } finally {
      for (const control of mounted.splice(0))
        control.dispose();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function mount(control) {
  mounted.push(control);
  document.body.append(control.root);
  return control;
}

const footer = (dialog) => [...dialog.root.querySelectorAll('.u2-dialog-buttons button')]
  .map((b) => b.textContent.trim());

smoke('the footer reads as the platform\'s: CANCEL first, OK right-most and flat, other buttons between', () => {
  const dialog = Dialog.create('Order').onOK(() => {}).onCancel(() => {}).addButton('APPLY', () => {}).show();
  assert.deepEqual(footer(dialog), ['CANCEL', 'APPLY', 'OK']);
  assert.equal(dialog.getButton('OK').classList.contains('u2-btn-primary'), false);
  dialog.close();
  const late = Dialog.create('Cancel first').onCancel(() => {}).onOK(() => {}).show();
  assert.deepEqual(footer(late), ['CANCEL', 'OK']);
  late.close();
});

smoke('OK follows the validity of the form in the body, Enter included', async () => {
  const name = new TextInput({label: 'Name', value: 'u2'});
  name.addValidator((v) => v.trim() ? null : 'Required');
  let accepted = 0;
  const dialog = Dialog.create('Gate', {keepAlive: true}).add(new Form().add(name)).onOK(() => accepted++).show();
  const ok = dialog.getButton('OK');
  assert.equal(ok.disabled, false);
  name.value = '';
  await flush();
  assert.equal(ok.disabled, true, 'an invalid form disables OK');
  fire(dialog.root.querySelector('input'), 'keydown', {key: 'Enter'});
  assert.equal(accepted, 0, 'Enter does not accept it either');
  assert.equal(dialog.isOpen.value, true);
  name.value = 'fixed';
  await flush();
  assert.equal(ok.disabled, false);
  fire(ok, 'click');
  assert.equal(accepted, 1);
  dialog.dispose();
});

smoke('a shown dialog comes back on screen when the window shrinks; a centred one stays centred', () => {
  const dialog = Dialog.create('Resize').onOK(() => {}).show();
  layout(dialog.root, {rect: new DOMRect(0, 0, 200, 100)});
  window.innerWidth = 1600;
  window.innerHeight = 1000;
  fire(window, 'resize', {bubbles: false});
  assert.equal(dialog.root.style.left, '700px');
  assert.equal(dialog.root.style.top, '450px');
  window.innerWidth = 500;
  window.innerHeight = 300;
  fire(window, 'resize', {bubbles: false});
  assert.equal(dialog.root.style.left, '150px');
  assert.equal(dialog.root.style.top, '100px');
  dialog.close();
  window.innerWidth = 1024;
  window.innerHeight = 768;
});

smoke('a shown dialog whose content grows comes back on screen, as the platform\'s re-clips', () => {
  const dialog = Dialog.create('Grows').onOK(() => {}).show();
  layout(dialog.root, {rect: new DOMRect(0, 0, 200, 100)});
  const observer = ResizeObserver.instances.find((o) => o.targets.includes(dialog.root) && !o.disconnected);
  observer.callback([]);
  assert.equal(dialog.root.style.top, '334px');
  layout(dialog.root, {rect: new DOMRect(0, 0, 200, 700)});
  observer.callback([]);
  assert.equal(dialog.root.style.top, '34px', 'the taller dialog is centred again, its footer inside the window');
  dialog.close();
});

smoke('a dialog dragged right is measured at the left edge, so a narrower window does not crush it', () => {
  const dialog = Dialog.create('Dragged').onOK(() => {}).show({x: 1300, y: 100});
  // a fixed box shrinks to the room right of its left edge
  dialog.root.getBoundingClientRect = () => {
    const left = parseFloat(dialog.root.style.left) || 0;
    return new DOMRect(left, 0, Math.max(100, Math.min(300, window.innerWidth - left)), 150);
  };
  window.innerWidth = 700;
  fire(window, 'resize', {bubbles: false});
  assert.equal(dialog.root.style.left, '400px', 'its full 300px fits against the right edge');
  window.innerWidth = 1024;
  dialog.close();
});

smoke('show selects the text of the first input, as the platform dialog does', () => {
  const name = new TextInput({label: 'Name', value: 'u2'});
  const input = name.root.querySelector('input');
  let selected = 0;
  input.select = () => selected++;
  const dialog = Dialog.create('Select').add(name).onOK(() => {}).show();
  assert.equal(document.activeElement, input);
  assert.equal(selected, 1);
  dialog.close();
});

smoke('a choice straight in a dialog body takes the width a form would measure for it, as ui.dialog\'s', () => {
  const dataset = new ChoiceInput({label: 'Dataset', name: 'dataset', items: ['spgi', 'peptides'], value: 'spgi'});
  const dialog = Dialog.create('Open').add(dataset).onOK(() => {});
  const body = dialog.root.querySelector('.u2-dialog-content');
  const append = body.append.bind(body);
  body.append = (el) => {
    el.getBoundingClientRect = () => ({width: el.textContent.length * 7});
    append(el);
  };
  dialog.show();
  assert.equal(dataset.root.style.getPropertyValue('--u2-choice-width'), `${'peptides'.length * 7 + 30}px`);
  dialog.close();
});

smoke('onHelp adds the ? beside the close button and runs on its click', () => {
  let helped = 0;
  const plain = Dialog.create('Plain').onOK(() => {}).show();
  assert.equal(plain.getButton('Help'), undefined, 'no ? without a help');
  plain.close();
  const dialog = Dialog.create('Help').onHelp(() => helped++).onOK(() => {}).show();
  const help = dialog.getButton('Help');
  assert.equal(help.nextSibling, dialog.getButton('Close'));
  fire(help, 'click');
  assert.equal(helped, 1);
  assert.equal(dialog.isOpen.value, true);
  dialog.close();
});

smoke('an edit in the body clears the error row of a failed OK', () => {
  const dialog = Dialog.create('Fail', {keepAlive: true}).add(new TextInput({label: 'Name'}))
    .onOK(() => { throw new Error('Server said no'); }).show();
  fire(dialog.getButton('OK'), 'click');
  assert.equal(dialog.error.value, 'Server said no');
  fire(dialog.root.querySelector('input'), 'input');
  assert.equal(dialog.error.value, null);
  dialog.dispose();
});

smoke('a busy OK clicked with the mouse keeps the focus inside the dialog', async () => {
  let release;
  const dialog = Dialog.create('Busy', {keepAlive: true}).add(new TextInput({label: 'Name'}))
    .onOK(() => new Promise((resolve) => release = resolve)).show();
  const ok = dialog.getButton('OK');
  ok.focus();
  fire(ok, 'click');
  assert.equal(dialog.busy.value, true);
  assert.equal(document.activeElement, dialog.root);
  release(false);
  await flush();
  dialog.dispose();
});

smoke('a window resize closes the anchored popups, not the dialogs', () => {
  const anchor = button('Open', () => {});
  document.body.append(anchor);
  const scope = new Scope();
  const popup = document.createElement('div');
  Overlay.show(anchor, popup, {owner: scope});
  const dialog = Dialog.create('Stays').onOK(() => {}).show();
  assert.equal(popup.isConnected, true);
  fire(window, 'resize', {bubbles: false});
  assert.equal(popup.isConnected, false, 'the popup closed');
  assert.equal(dialog.isOpen.value, true, 'the dialog did not');
  dialog.close();
  scope.dispose();
});

smoke('a type-ahead drops text that is not a pick on blur, and does not spellcheck', async () => {
  const picker = mount(new TypeAhead({source: ['All users', 'Developers'], itemText: (s) => s}));
  const input = picker.root.querySelector('input');
  assert.equal(input.spellcheck, false);
  picker.selected.value = 'All users';
  input.focus();
  input.value = 'zzzz-nobody';
  fire(input, 'input');
  assert.equal(picker.selected.value, null);
  input.blur();
  await flush();
  assert.equal(input.value, '', 'the stray text is gone, so the box reads as empty as the value');
});

smoke('a new filter row says its value is missing only once it was left, or the builder validated', async () => {
  const schema = Filters.schema([{name: 'name', type: TYPE.STRING}]);
  const fb = mount(new FilterBuilder({schema}));
  const id = fb.addCondition();
  const row = fb.root.querySelector(`[data-u2-node="${id}"]`);
  const problem = () => row.querySelector('[data-u2-part="problem"]');
  assert.equal(problem().hidden, true, 'fresh rows are quiet');
  assert.equal(row.classList.contains('u2-fb-invalid'), false);
  assert.notEqual(fb.validity.peek(), null, 'the verdict itself stands');
  assert.equal(problem().parentElement.classList.contains('u2-fb-values'), true, 'the line sits under the value');
  fire(row.querySelector('.u2-fb-values input'), 'focusout');
  await flush();
  assert.equal(problem().hidden, false, 'a row left empty says so');
  const other = fb.addCondition();
  assert.equal(fb.root.querySelector(`[data-u2-node="${other}"] [data-u2-part="problem"]`).hidden, true);
  fb.touch();
  await flush();
  assert.equal(fb.root.querySelector(`[data-u2-node="${other}"] [data-u2-part="problem"]`).hidden, false,
    'touch shows every row');
});

smoke('a tree row reads as the platform tree\'s: 20px tall, a chevron for a twistie, 20px per level', () => {
  const tree = mount(new VirtualTree());
  layout(tree.root.querySelector('.u2-list'), {clientHeight: 200});
  tree.setRoots([{id: 'a', label: 'Alpha', children: [{id: 'a1', label: 'Alpha one'}]}]);
  const row = (i) => tree.root.querySelector(`.u2-list-row[data-index="${i}"]`);
  const twistie = (i) => row(i).querySelector('.u2-tree-twistie');
  assert.equal(twistie(0).classList.contains('fa-chevron-right'), true, 'a closed branch shows a chevron');
  assert.equal(twistie(0).textContent, '', 'not a text glyph');
  fire(twistie(0), 'click');
  assert.equal(twistie(0).classList.contains('fa-chevron-down'), true, 'an open one points down');
  assert.equal(row(1).style.height, '20px');
  assert.equal(row(1).querySelector('.u2-tree-row').style.getPropertyValue('--u2-tree-indent'), '20px');
  assert.equal(twistie(1).classList.contains('u2-tree-twistie-hidden'), true, 'a leaf keeps the box, hidden');
});

// the shim cannot lay out: the platform metrics of the display skins are pinned at the source
test('accordion / table / tree / tabs / filter css: the platform\'s metrics and states hold', () => {
  const css = (name) => readFileSync(new URL(`../css/${name}.css`, import.meta.url), 'utf8').replace(/\r\n/g, '\n');
  const block = (sheet, selector) => {
    const at = sheet.indexOf(selector);
    assert.notEqual(at, -1, `${selector} is a rule`);
    return sheet.slice(sheet.indexOf('{', at) + 1, sheet.indexOf('}', at));
  };
  const accordion = block(css('accordion'), '\n.u2-accordion-content {');
  assert.match(accordion, /padding: 2px 0 10px 20px;/, 'the open pane sits as .d4-accordion-pane-content does');
  assert.match(accordion, /color: var\(--dg-grey-5\);/);
  const table = css('table');
  assert.match(block(table, '.u2-table th,\n.u2-table td {'), /padding: 0 0 0 var\(--dg-space-l\);\n  border-bottom: none;/,
    'ui.table: 12px in, no rules');
  assert.match(block(table, '\n.u2-table th {'), /font-size: var\(--dg-font-size\);\n  font-weight: var\(--dg-font-weight-bold\);/);
  assert.match(block(css('tree'), '.u2-tree-row > .u2-tree-twistie {'), /width: 20px;/);
  const tabs = css('tabs');
  const strip = '.u2-tabs:not(.u2-tabs-document):not(.u2-tabs-vertical) > .u2-tabs-header';
  assert.match(block(tabs, `${strip} {`), /height: 34px;\n  box-shadow: none;/, 'a 34px stripe with no rule of its own');
  assert.match(block(tabs, `${strip} .u2-tabs-tab {`), /padding: var\(--dg-space-m\) var\(--dg-space-xl\);/);
  const filter = css('filter');
  assert.match(block(filter, '.u2-input-root[data-u2=\'filter-builder\'].u2-fb-blank > .u2-input-error:not(:empty) {'),
    /display: -webkit-box;/, 'a blank builder shows its verdict, whatever an older copy of the sheet says');
  assert.match(block(filter, '.u2-input-disabled .u2-fb .u2-fb-none-add {'), /color: var\(--dg-text-color-light\);/);
  assert.match(block(filter, '.u2-input-root[data-u2=\'filter-builder\'].u2-input-readonly .u2-fb-none-add {'),
    /display: none;/, 'read-only hides what edits the tree');
});
