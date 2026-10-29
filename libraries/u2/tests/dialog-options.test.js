/* Dialog show options and confirmDialog(): `visible` hides without closing, `owner` closes with it, a fixed or
   resizable height lets the body fill the dialog, the first focus skips read-only fields; confirmDialog() is CANCEL
   plus a verb and answers whether the user went ahead. Also the required asterisk (G-R1) and the AccessGrid
   headers that must not wrap (G-AG1). */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, Control, Dialog, Form, TextInput, signal, confirmDialog} from '../src/index.js';

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

const css = (sheet) => readFileSync(new URL(`../css/${sheet}`, import.meta.url), 'utf8');

smoke('visible: false hides the dialog and a modal\'s backdrop, without closing it; true shows it again', () => {
  const visible = signal(true);
  const dialog = new Dialog('Progress tracker').onOK(() => {}).show({modal: true, visible});
  const backdrop = document.querySelector('.u2-dialog-backdrop');
  assert.equal(dialog.root.hidden, false);
  assert.equal(dialog.modal, true);
  visible.value = false;
  assert.equal(dialog.root.hidden, true);
  assert.equal(backdrop.hidden, true);
  assert.equal(dialog.isOpen.value, true, 'hidden, not closed');
  assert.equal(dialog.modal, false, 'a hidden dialog blocks nothing and takes no Escape');
  visible.value = true;
  assert.equal(dialog.root.hidden, false);
  dialog.close();
  visible.value = false;
  assert.equal(dialog.root.hidden, false, 'the follow ended with the showing');
});

smoke('status: a quiet footer line says what OK does, the reason while OK is off; start buttons sit before it', () => {
  const reason = signal(null);
  const moved = signal('Drag a card to another lane to change its stage');
  const dialog = new Dialog('Progress Tracker').onOK(() => {}).okDisabled(reason).status(moved)
    .addButton('Modify Stages', () => {}, {start: true}).show();
  const footer = dialog.root.querySelector('[data-u2-part="footer"]');
  const order = () => [...footer.children].map((el) => el.classList.contains('u2-dialog-reason') ? 'line' :
    el.textContent.trim());
  assert.deepEqual(order(), ['Modify Stages', 'line', 'CANCEL', 'OK']);
  const line = dialog.root.querySelector('[data-u2-part="reason"]');
  assert.equal(line.textContent, 'Drag a card to another lane to change its stage');
  assert.equal(line.classList.contains('u2-status-line-problem'), false, 'quiet while OK can run');
  moved.value = '2 molecules moved · OK saves the campaign, CANCEL puts them back';
  assert.equal(line.textContent, '2 molecules moved · OK saves the campaign, CANCEL puts them back');
  reason.value = 'Nothing to save';
  assert.equal(line.textContent, 'Nothing to save');
  assert.equal(line.classList.contains('u2-status-line-problem'), true);
  dialog.close();
});

smoke('owner: the dialog closes when it is disposed; closed first, it leaves the owner nothing', () => {
  const session = new Control();
  const dialog = new Dialog('Merge').onOK(() => {}).show({owner: session});
  session.dispose();
  assert.equal(dialog.isOpen.value, false);
  assert.equal(dialog.scope.isDisposed, true);

  const scope = new Scope();
  const kept = new Dialog('Kept', {keepAlive: true}).onOK(() => {});
  kept.show({owner: scope});
  kept.close();
  scope.dispose();
  assert.equal(kept.scope.isDisposed, false, 'a closed dialog is no longer the owner\'s to end');
  kept.dispose();
});

smoke('the close event carries the reason', () => {
  const reasons = [];
  const dialog = new Dialog('Edit').onOK(() => {});
  dialog.on('close', (reason) => reasons.push(reason));
  dialog.show();
  fire(dialog.getButton('OK'), 'click');
  assert.deepEqual(reasons, ['ok']);
});

smoke('a fixed or resizable height lets the body fill it; the body\'s own cap applies otherwise', () => {
  const fixed = new Dialog('Tracker').onOK(() => {}).show({height: 600});
  assert.equal(fixed.root.classList.contains('u2-dialog-sized'), true);
  assert.equal(fixed.root.style.height, '600px');
  const resizable = new Dialog('Compute').onOK(() => {}).show({width: 960, height: 640, resizable: true});
  assert.equal(resizable.root.classList.contains('u2-dialog-resizable'), true);
  const plain = new Dialog('Plain').onOK(() => {}).show();
  assert.equal(plain.root.classList.contains('u2-dialog-sized'), false);
  const sheet = css('dialog.css');
  assert.match(sheet, /\.u2-dialog-sized > \.u2-dialog-content \{\s*max-height: none;/);
  assert.match(sheet, /\.u2-dialog-resizable \{[^}]*resize: both;/);
  fixed.close();
  resizable.close();
  plain.close();
});

smoke('the first focus goes to the first field that takes edits, never a read-only one', () => {
  const name = new TextInput({label: 'Name', value: 'KHL-3', readOnly: true});
  const status = new TextInput({label: 'Status', value: 'Submitted'});
  const dialog = new Dialog('Submit').add(name).add(status).onOK(() => {}).show();
  assert.ok(document.activeElement === status.root.querySelector('input'), 'the editable Status, not the read-only Name');
  dialog.close();
});

smoke('confirm: the verb submits and answers true; CANCEL answers false', async () => {
  const yes = confirmDialog({title: 'Delete campaign', message: 'Are you sure you want to delete campaign KHL-3?',
    detail: 'This cannot be undone.', confirmText: 'Delete', destructive: true, name: 'delete-campaign'});
  const dialog = Control.forElement(document.querySelector('[data-u2-name="delete-campaign"]'));
  assert.equal(dialog.getState().label, 'Delete campaign');
  assert.ok(dialog.getButton('CANCEL'), 'CANCEL comes with the verb');
  assert.ok(dialog.root.querySelector('.u2-confirm-icon.fa-trash-alt'), 'a destructive action shows the trash');
  assert.equal(dialog.root.querySelector('.u2-confirm-detail').textContent, 'This cannot be undone.');
  assert.ok(document.activeElement === dialog.getButton('CANCEL'), 'CANCEL has the first focus, never the action');
  fire(dialog.getButton('Delete'), 'click');
  assert.equal(await yes, true);

  const no = confirmDialog({title: 'Delete template', message: 'Are you sure?', confirmText: 'Delete'});
  assert.ok(document.activeElement === Control.forElement(document.querySelector('.u2-dialog')).getButton('CANCEL'),
    'whether the action is destructive or not');
  fire(Control.forElement(document.querySelector('.u2-dialog')).getButton('CANCEL'), 'click');
  assert.equal(await no, false);
});

smoke('confirm with work: busy until it settles, a failure stays in the dialog, the answer comes after', async () => {
  let fail = true;
  let finish;
  const answer = confirmDialog({title: 'Delete campaign', message: 'Sure?', confirmText: 'Delete',
    onConfirm: () => new Promise((resolve, reject) => finish = () => fail ? reject(new Error('No access')) : resolve())});
  const dialog = Control.forElement(document.querySelector('.u2-dialog'));
  fire(dialog.getButton('Delete'), 'click');
  assert.equal(dialog.busy.value, true);
  finish();
  await flush();
  assert.equal(dialog.error.value, 'No access');
  assert.equal(dialog.isOpen.value, true);
  fail = false;
  fire(dialog.getButton('Delete'), 'click');
  finish();
  assert.equal(await answer, true);
  assert.equal(dialog.isOpen.value, false);
});

smoke('closeOnly: one CLOSE and never a CANCEL; CLOSE, ✕ and Escape all close as OK, keeping what was done', () => {
  const closed = [];
  const dialog = new Dialog('Progress tracker', {closeOnly: true})
    .addButton('Modify Stages', () => {}, {start: true})
    .onClose((reason) => closed.push(reason));
  assert.throws(() => dialog.onCancel(() => closed.push('cancel ran')), {code: 'close-only-cancel'},
    'there is no CANCEL to hang it on');
  dialog.show();
  const footer = [...dialog.root.querySelectorAll('[data-u2-part="footer"] button')].map((b) => b.textContent);
  assert.deepEqual(footer, ['Modify Stages', 'CLOSE'], 'a footer button brings no CANCEL');
  assert.equal(dialog.getButton('CANCEL'), undefined);
  fire(dialog.root, 'keydown', {key: 'Escape'});
  assert.deepEqual(closed, ['ok'], 'Escape closes as OK; nothing is cancelled');

  const again = new Dialog('Progress tracker', {closeOnly: true, keepAlive: true}).onClose((r) => closed.push(r));
  again.show();
  fire(again.root.querySelector('[data-u2-part="close"]'), 'click');
  again.show();
  fire(again.getButton('CLOSE'), 'click');
  assert.deepEqual(closed, ['ok', 'ok', 'ok'], '✕ and CLOSE alike');
  again.dispose();
});

smoke('closeOnly: while CLOSE is off, Escape and ✕ close as a cancel and never run onOK; getButton(\'Close\') is CLOSE',
  () => {
    const ran = [];
    const held = new Dialog('Tracker', {closeOnly: true}).onOK(() => ran.push('ok'))
      .okDisabled('Wait for the upload to finish').onClose((r) => ran.push(r));
    held.show();
    const close = held.root.querySelector('[data-u2-part="footer"] button');
    assert.equal(held.getButton('Close'), close, 'the footer\'s CLOSE, not the title bar\'s ✕');
    fire(close, 'click');
    assert.deepEqual(ran, [], 'CLOSE is off');
    fire(held.root, 'keydown', {key: 'Escape'});
    assert.deepEqual(ran, ['cancel'], 'closed, onOK never ran');

    const saved = [];
    const form = new Dialog('Settings', {closeOnly: true}).onOK(() => saved.push('saved'));
    const f = new Form();
    f.add(new TextInput({label: 'Name', required: true}));
    form.add(f).show();
    fire(form.root.querySelector('[data-u2-part="close"]'), 'click');
    assert.deepEqual(saved, [], 'an invalid form is never saved by ✕');
    assert.equal(form.isOpen.value, false);
  });

smoke('a footer button whose own reason turned on while OK was busy keeps it when OK settles', async () => {
  const off = signal(null);
  let finish;
  const dialog = new Dialog('Run').addButton('Export', () => {}, {disabled: off})
    .onOK(() => new Promise((resolve) => finish = resolve)).show();
  const exp = dialog.getButton('Export');
  fire(dialog.getButton('OK'), 'click');
  assert.equal(exp.getAttribute('aria-description'), Dialog.WAITING);
  off.value = 'Nothing to export yet';
  finish(false);
  await flush();
  assert.equal(dialog.busy.value, false);
  assert.equal(exp.getAttribute('aria-disabled'), 'true', 'still off');
  assert.equal(exp.getAttribute('aria-description'), 'Nothing to export yet');
  off.value = null;
  assert.equal(exp.getAttribute('aria-disabled'), null, 'and on with its reason gone');
  dialog.close();
});

smoke('okDisabled takes a reason: OK is off saying why, and on again with null', () => {
  const reason = signal('Pick a viewer');
  const dialog = new Dialog('Add viewer').onOK(() => {}).okDisabled(reason).show();
  const ok = dialog.getButton('OK');
  assert.equal(ok.getAttribute('aria-disabled'), 'true');
  assert.equal(ok.getAttribute('aria-description'), 'Pick a viewer');
  reason.value = null;
  assert.equal(ok.hasAttribute('aria-disabled'), false);
  assert.equal(ok.hasAttribute('aria-description'), false);
  dialog.close();
});

smoke('required: true puts the asterisk on the caption; a live flag follows', () => {
  const required = signal(true);
  const input = new TextInput({label: 'Name', required});
  assert.equal(input.root.classList.contains('u2-input-required'), true);
  required.value = false;
  assert.equal(input.root.classList.contains('u2-input-required'), false);
  const key = new TextInput({label: 'Key', nullable: false});
  assert.equal(key.root.classList.contains('u2-input-required'), false, 'nullable: false refuses empty, unmarked');
  input.dispose();
  key.dispose();
});

test('the AccessGrid\'s headers never wrap ("Can View" on one line)', () => {
  assert.match(css('access-grid.css'), /\.u2-access-grid th \{[^}]*white-space: nowrap;/);
});
