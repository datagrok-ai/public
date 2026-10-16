/* Dialog lifecycle (R8): closing disposes the dialog and what it holds unless `keepAlive`; an async
   OK keeps it open and busy, shows a failure in the error row and aborts on Escape; button handles,
   the OK caption and `onClose`; the one overlay stack routing Escape (P9) with one listener set
   (P24); the dialog band shared with the platform's dialogs and raise on press (P76). */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {debug} from 'datagrok-api/u2core';
import {Scope, Control, Overlay, Dialog, Form, TextInput, NumberInput, Combobox, PropertyGrid, button, divH, span, signal,
  computed, Tooltip}
  from '../src/index.js';

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

const press = (dialog, text) => fire(dialog.getButton(text), 'click');

smoke('every close path disposes the dialog and the form built with it in a click handler (P10)', () => {
  const paths = [
    (d) => press(d, 'OK'),
    (d) => press(d, 'CANCEL'),
    (d) => press(d, 'Close'),
    (d) => fire(d.root, 'keydown', {key: 'Escape'}),
    (d) => fire(d.root, 'keydown', {key: 'Enter'}),
    (d) => d.close(),
  ];
  let dialog;
  const page = Control.build(() => [button('Edit', () => {
    const form = new Form().addAll([new TextInput({label: 'Name'}), new NumberInput({label: 'Dose'})]);
    dialog = Dialog.create('Edit compound').add(form).onOK(() => {}).onCancel(() => {}).show();
  })]);
  document.body.append(page.root);
  const open = page.root.querySelector('button');
  const base = Scope.liveCount;
  for (let i = 0; i < 12; i++) {
    fire(open, 'click');
    assert.equal(dialog.isOpen.value, true);
    paths[i % paths.length](dialog);
    assert.equal(dialog.isOpen.value, false, `path ${i % paths.length} closes`);
    assert.equal(dialog.scope.isDisposed, true, `path ${i % paths.length} disposes`);
  }
  assert.equal(Scope.liveCount, base, 'nothing of the twelve dialogs is left');
  assert.equal(document.querySelector('.u2-dialog'), null);
  page.dispose();
});

smoke('content built outside any scope is freed with the dialog; a control that never joins one is left alone', () => {
  const name = new TextInput({label: 'Name'});
  class PlatformWidget extends Control {
    get ambientAdoption() {
      return false;
    }
  }
  const widget = new PlatformWidget();
  Dialog.create('Top level').add(name).add(widget).show().close();
  assert.equal(name.scope.isDisposed, true);
  assert.equal(widget.scope.isDisposed, false, 'a platform widget is for its own owner to kill');
  widget.dispose();
});

smoke('keepAlive: a closed dialog shows again with its content alive; dispose frees both', () => {
  const name = new TextInput({label: 'Name', value: 'u2'});
  const dialog = Dialog.create('Again', {keepAlive: true}).add(name).onOK(() => {});
  dialog.show();
  press(dialog, 'OK');
  assert.equal(dialog.isOpen.value, false);
  assert.equal(dialog.scope.isDisposed, false, 'kept');
  name.value.value = 'v2';
  dialog.show();
  assert.equal(dialog.isOpen.value, true);
  assert.equal(name.root.querySelector('input').value, 'v2', 'its bindings still run');
  dialog.dispose();
  assert.equal(dialog.root.isConnected, false);
  assert.equal(name.scope.isDisposed, true, 'the content goes with it');
});

smoke('a disposed dialog does not show again, and says why', () => {
  const dialog = Dialog.create('Once').show();
  dialog.close();
  const warn = console.warn;
  const said = [];
  console.warn = (m) => said.push(m);
  try {
    dialog.show();
  } finally {
    console.warn = warn;
  }
  assert.equal(dialog.isOpen.value, false);
  assert.match(said[0], /keepAlive/);
});

smoke('a dialog shown from a row button outlives the row and is freed by close; never shown, it dies with it', () => {
  const row = new Scope();
  let dialog;
  let name;
  const edit = Scope.runWith(row, () => button('Edit…', () => {
    name = new TextInput({label: 'Name'});
    dialog = Dialog.create('Edit row').add(name).onOK(() => {}).show({modal: true});
  }));
  document.body.append(edit);
  fire(edit, 'click');
  row.dispose();
  assert.deepEqual([dialog.isOpen.value, dialog.scope.isDisposed, name.scope.isDisposed], [true, false, false],
    'the row scrolled away: the dialog and its content stay');
  press(dialog, 'OK');
  assert.deepEqual([dialog.scope.isDisposed, name.scope.isDisposed], [true, true], 'closing frees both');

  const page = new Scope();
  const unshown = Scope.runWith(page, () => Dialog.create('Never shown'));
  page.dispose();
  assert.equal(unshown.scope.isDisposed, true);
});

smoke('plain elements passed as add(() => …) from a row keep their handlers and live text once the row is gone', () => {
  const row = new Scope();
  const count = signal(3);
  let dialog;
  const edit = Scope.runWith(row, () => button('Edit…', () => {
    dialog = Dialog.create('Edit row')
      .add(() => divH([button('Reset', () => count.value = 0), span(count, 'count')]))
      .onOK(() => {}).show();
  }));
  document.body.append(edit);
  fire(edit, 'click');
  row.dispose();
  const text = () => dialog.root.querySelector('.count').textContent;
  count.value = 5;
  assert.equal(text(), '5', 'the text still follows its signal');
  fire(dialog.root.querySelector('.u2-dialog-content button'), 'click');
  assert.deepEqual([count.value, text()], [0, '0'], 'the button still answers');
  dialog.close();
  count.value = 7;
  assert.equal(text(), '0', 'closing freed the binding');
});

smoke('what the dialog adopted keeps its callbacks once the scope that built it is gone', () => {
  const builder = new Scope();
  const seen = [];
  const {name, grid, dialog} = Scope.runWith(builder, () => {
    const name = new TextInput({label: 'Name', onChanged: (v) => seen.push(`changed ${v}`),
      onEnter: () => seen.push('enter')});
    const grid = new PropertyGrid({onChanged: (c) => seen.push(`${c.name}=${c.value}`)});
    grid.setProperties([{name: 'dose', type: 'string'}], {dose: ''});
    return {name, grid, dialog: Dialog.create('Edit').add(name).add(grid).onOK(() => {}).show()};
  });
  builder.dispose();
  name.value.value = 'Ada';
  fire(name.root.querySelector('input'), 'keydown', {key: 'Enter'});
  const field = grid.root.querySelector('.u2-propgrid-row input');
  field.value = '5';
  fire(field, 'input');
  assert.deepEqual(seen, ['changed Ada', 'enter', 'dose=5']);
  dialog.close();
  assert.equal(name.scope.isDisposed, true);
});

smoke('a keepAlive dialog stays with the scope that built it, and closes with it', () => {
  const page = new Scope();
  const dialog = Scope.runWith(page, () => Dialog.create('Owned', {keepAlive: true}).show({modal: true}));
  assert.equal(dialog.root.isConnected, true);
  page.dispose();
  assert.deepEqual([dialog.isOpen.value, dialog.scope.isDisposed], [false, true]);
  assert.equal(document.querySelector('.u2-dialog, .u2-dialog-backdrop'), null);
});

smoke('what an OK handler builds belongs to whoever registered it, not to the dialog it closes', () => {
  const page = new Scope();
  let made;
  const dialog = Scope.runWith(page, () => Dialog.create('Make').onOK(() => made = new TextInput({label: 'Made'})));
  dialog.show();
  press(dialog, 'OK');
  assert.equal(dialog.isOpen.value, false);
  assert.equal(made.scope.isDisposed, false, 'outlives the dialog');
  page.dispose();
  assert.equal(made.scope.isDisposed, true, 'and goes with its owner');

  const registrant = new Scope();
  let ran = 0;
  const late = Dialog.create('Late');
  Scope.runWith(registrant, () => late.onOK(() => ran++));
  registrant.dispose();
  late.show();
  press(late, 'OK');
  assert.deepEqual([ran, late.isOpen.value], [1, false], 'a handler whose registrant is gone still runs');
});

smoke('an async OK keeps the dialog open and busy until it resolves, then closes (P75)', async () => {
  let resolve;
  let calls = 0;
  const dialog = Dialog.create('Save').onOK(() => {
    calls++;
    return new Promise((r) => resolve = r);
  }).onCancel(() => {}).addButton('MORE', () => {});
  dialog.show();
  const ok = dialog.getButton('OK');
  fire(ok, 'click');
  assert.equal(dialog.isOpen.value, true, 'open while saving');
  assert.equal(dialog.busy.value, true);
  assert.equal(dialog.root.getAttribute('aria-busy'), 'true');
  assert.equal(ok.getAttribute('aria-busy'), 'true');
  assert.notEqual(ok.querySelector('[data-u2-part="busy"]'), null, 'a spinner on OK');
  assert.deepEqual(['OK', 'CANCEL', 'MORE'].map((t) => dialog.getButton(t).disabled), [true, true, true]);
  assert.equal(dialog.getState().busy, true);
  fire(ok, 'click');
  fire(dialog.root, 'keydown', {key: 'Enter'});
  assert.equal(calls, 1, 'a second press while busy is ignored');
  resolve();
  await flush();
  assert.equal(dialog.isOpen.value, false);
  assert.equal(dialog.scope.isDisposed, true);
});

smoke('a rejected or throwing OK stays open with the message in the error row; the next OK clears it', async () => {
  let attempt = 0;
  const dialog = Dialog.create('Save').onOK(async () => {
    if (++attempt === 1)
      throw new Error('Server said no');
  });
  dialog.show();
  const row = dialog.root.querySelector('[data-u2-part="error-row"]');
  assert.equal(row.hidden, true);
  press(dialog, 'OK');
  await flush();
  assert.equal(dialog.isOpen.value, true);
  assert.equal(dialog.busy.value, false);
  assert.equal(dialog.error.value, 'Server said no');
  assert.equal(row.hidden, false);
  assert.equal(row.textContent, 'Server said no');
  assert.equal(row.getAttribute('role'), 'alert');
  assert.equal(dialog.getButton('OK').disabled, false, 'the buttons are back');
  assert.equal(dialog.root.hasAttribute('aria-busy'), false);
  assert.deepEqual([dialog.getState().error, dialog.getState().open], ['Server said no', true]);
  press(dialog, 'OK');
  assert.equal(dialog.error.value, null, 'cleared on the next OK');
  await flush();
  assert.equal(dialog.isOpen.value, false);

  const sync = Dialog.create('Sync').onOK(() => {
    throw new Error('Not like that');
  });
  sync.show();
  press(sync, 'OK');
  assert.deepEqual([sync.isOpen.value, sync.error.value], [true, 'Not like that']);
  sync.dispose();

  let answer = false;
  const vetoing = Dialog.create('Veto').onOK(async () => answer);
  vetoing.show();
  press(vetoing, 'OK');
  await flush();
  assert.deepEqual([vetoing.isOpen.value, vetoing.error.value], [true, null], 'resolving to false stays, quietly');
  answer = true;
  press(vetoing, 'OK');
  await flush();
  assert.equal(vetoing.isOpen.value, false);
});

smoke('Escape or ✕ while busy aborts the OK and keeps the dialog; the next Escape cancels', async () => {
  let signal;
  let cancelled = 0;
  const dialog = Dialog.create('Upload').onOK(({signal: s}) => {
    signal = s;
    return new Promise((_, reject) => s.addEventListener('abort', () => reject(s.reason)));
  }).onCancel(() => cancelled++);
  dialog.show();
  press(dialog, 'OK');
  fire(dialog.root, 'keydown', {key: 'Escape'});
  assert.equal(signal.aborted, true);
  assert.equal(cancelled, 0, 'aborting is not cancelling');
  await flush();
  assert.deepEqual([dialog.isOpen.value, dialog.busy.value, dialog.error.value], [true, false, null],
    'still open, idle, and an abort shows no error');
  press(dialog, 'OK');
  press(dialog, 'Close');
  assert.equal(signal.aborted, true, '✕ aborts too');
  await flush();
  assert.equal(dialog.isOpen.value, true);
  fire(dialog.root, 'keydown', {key: 'Escape'});
  assert.deepEqual([dialog.isOpen.value, cancelled], [false, 1]);
});

smoke('Escape in an open combobox inside a dialog closes only the popup (P9)', async () => {
  const combobox = new Combobox({items: ['alpha', 'beta']});
  const dialog = Dialog.create('Pick').add(combobox).onOK(() => {}).onCancel(() => {});
  dialog.show();
  const input = combobox.root.querySelector('input');
  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();
  assert.equal(combobox.isOpen.value, true);
  const popup = document.querySelector('.u2-combobox-popup');
  assert.equal(popup.dataset.u2Part, 'popup');
  assert.equal(input.getAttribute('aria-controls'), popup.id);
  assert.equal(input.getAttribute('aria-expanded'), 'true');
  assert.equal(fire(input, 'keydown', {key: 'Escape'}), false, 'the popup took it');
  assert.equal(combobox.isOpen.value, false);
  assert.equal(dialog.isOpen.value, true, 'the dialog never saw it');
  assert.equal(input.getAttribute('aria-expanded'), 'false');
  fire(input, 'keydown', {key: 'Escape'});
  assert.equal(dialog.isOpen.value, false, 'with the popup closed, Escape cancels the dialog');
});

smoke('getButton finds OK by its caption, CANCEL, custom buttons and ✕; addButton passes the event and the dialog', () => {
  const pressed = [];
  const dialog = Dialog.create('Rename').onOK(() => false, {text: 'SAVE'}).onCancel(() => {})
    .addButton('Apply', (e, d) => pressed.push([e.type, d]));
  assert.equal(dialog.getButton('save').textContent, 'SAVE');
  assert.equal(dialog.getButton('OK'), undefined);
  assert.equal(dialog.getButton('Cancel').textContent, 'CANCEL');
  assert.equal(dialog.getButton('close'), dialog.root.querySelector('[data-u2-part="close"]'));
  dialog.show();
  press(dialog, 'APPLY');
  assert.deepEqual(pressed, [['click', dialog]]);
  assert.equal(dialog.isOpen.value, true, 'a footer button does not close by itself');
  dialog.getButton('save').disabled = true;
  fire(dialog.root, 'keydown', {key: 'Enter'});
  assert.equal(dialog.isOpen.value, true, 'a disabled OK takes no Enter');
  dialog.dispose();
});

smoke('leaving an empty required field for CANCEL shows no verdict before the press lands', () => {
  const name = new TextInput({label: 'Name', name: 'name', required: true});
  const dialog = Dialog.create('New collection').add(new Form().add(name)).onOK(() => {});
  dialog.show();
  const cancel = dialog.getButton('CANCEL');
  const error = () => name.root.querySelector('[data-u2-part="error"]').textContent;
  fire(name.root.querySelector('input'), 'focusout', {relatedTarget: cancel});
  assert.equal(error(), '', 'the line would move CANCEL from under the pointer');
  fire(cancel, 'click');
  assert.equal(dialog.isOpen.value, false, 'one click closes');

  const other = new TextInput({label: 'Other', name: 'other', required: true});
  document.body.append(other.root);
  fire(other.root.querySelector('input'), 'focusout', {relatedTarget: document.body});
  assert.notEqual(other.root.querySelector('[data-u2-part="error"]').textContent, '', 'left for anything else, it shows');
  other.dispose();
});

smoke('closing the dialog on top gives focus to the one below when what held it is gone', () => {
  const below = Dialog.create('Below').onOK(() => {});
  below.show();
  const top = Dialog.create('Top').onOK(() => {});
  const opener = document.createElement('button');
  document.body.append(opener);
  opener.focus();
  top.show();
  opener.remove();
  top.close();
  assert.equal(document.activeElement, below.root, 'Escape still reaches the dialog below');

  // opened from the page (a ribbon icon takes no focus): the opener reads as the body
  document.activeElement?.blur?.();
  const fromPage = Dialog.create('From the ribbon').onOK(() => {});
  fromPage.show();
  fromPage.close();
  assert.equal(document.activeElement, below.root, 'not the page under both');
  below.dispose();
});

smoke('CANCEL comes with the first action, as in a platform dialog, and closes as a cancel', () => {
  const reasons = [];
  const withOk = Dialog.create('Remove').onOK(() => {}, {text: 'Remove'}).onClose((r) => reasons.push(r));
  assert.deepEqual(Array.from(withOk.root.querySelectorAll('button')).filter((b) => b.dataset.u2Part !== 'close')
    .map((b) => b.textContent).filter((t) => t), ['CANCEL', 'Remove']);
  withOk.show();
  press(withOk, 'CANCEL');
  assert.deepEqual(reasons, ['cancel']);
  const withButton = Dialog.create('Exists').addButton('Merge', () => {});
  assert.equal(withButton.getButton('Cancel').textContent, 'CANCEL');
  const bare = Dialog.create('Note');
  assert.equal(bare.getButton('Cancel'), undefined, 'a dialog without actions has no footer to cancel from');
  withButton.dispose();
  bare.dispose();
});

smoke('onClose hears why: ok, cancel (CANCEL, ✕, Escape) or code; a handler that disposes keeps its reason', () => {
  const reasons = [];
  const dialog = Dialog.create('Why', {keepAlive: true}).onOK(() => {}).onCancel(() => {})
    .onClose((r) => reasons.push(r));
  for (const close of [
    () => press(dialog, 'OK'),
    () => fire(dialog.root, 'keydown', {key: 'Enter'}),
    () => press(dialog, 'CANCEL'),
    () => press(dialog, 'Close'),
    () => fire(dialog.root, 'keydown', {key: 'Escape'}),
    () => dialog.close(),
  ]) {
    dialog.show();
    close();
  }
  assert.deepEqual(reasons, ['ok', 'ok', 'cancel', 'cancel', 'cancel', 'code']);
  dialog.dispose();
  assert.deepEqual(reasons, ['ok', 'ok', 'cancel', 'cancel', 'cancel', 'code'], 'disposing a closed dialog says nothing');

  const heard = [];
  const self = Dialog.create('Self').onOK(() => self.dispose()).onClose((r) => heard.push(r));
  self.show();
  press(self, 'OK');
  assert.deepEqual(heard, ['ok']);
});

smoke('the state block: label, open, busy, error and the validity of a form inside', () => {
  const name = new TextInput({label: 'Name'});
  name.addValidator((v) => v.trim() ? null : 'Required');
  const dialog = Dialog.create('Rename', {name: 'rename'}).add(new Form().add(name));
  const closed = dialog.getState();
  assert.deepEqual([closed.kind, closed.name, closed.label, closed.open, closed.busy, closed.error],
    ['dialog', 'rename', 'Rename', false, false, null]);
  dialog.show();
  assert.equal(dialog.getState().open, true);
  assert.equal(dialog.getState().valid, false, 'the empty required field');
  name.value.value = 'Ada';
  assert.equal(dialog.getState().valid, true);
  dialog.dispose();
});

smoke('one keydown and one pointerdown listener serve every popup and dialog, gone once all close (P24)', () => {
  debug.enable({orphanSeconds: 0});
  try {
    const scope = new Scope();
    const a = Scope.runWith(scope, () => Dialog.create('A').show());
    const b = Scope.runWith(scope, () => Dialog.create('B').show({modal: true}));
    const closes = [];
    for (let i = 0; i < 3; i++) {
      const anchor = document.createElement('button');
      b.add(anchor);
      closes.push(Overlay.show(anchor, document.createElement('div'), scope));
    }
    const global = (type) => debug.listeners().filter((l) => l.type === type);
    assert.equal(global('keydown').length, 1);
    assert.equal(global('pointerdown').length, 1);
    assert.deepEqual([global('keydown')[0].target, global('keydown')[0].capture, global('keydown')[0].label],
      ['document', true, 'singleton']);
    for (const close of closes)
      close();
    a.close();
    assert.equal(global('keydown').length, 1, 'B is still open');
    b.close();
    assert.equal(global('keydown').length + global('pointerdown').length, 0);
    scope.dispose();
  } finally {
    debug.disable();
  }
});

smoke('dialogs share the platform dialog band: above a DG.Dialog at 3500, out of the overlay host (P76)', () => {
  assert.equal(Overlay.nextDialogZ(), 3000, 'the band starts at --dg-z-dialog');
  const platform = document.createElement('div');
  platform.className = 'd4-dialog';
  platform.style.zIndex = '3500';
  document.body.append(platform);
  assert.equal(Overlay.nextDialogZ(), 3501);
  const dialog = Dialog.create('Above').show({modal: true});
  const backdrop = document.querySelector('.u2-dialog-backdrop');
  assert.equal(dialog.root.parentElement, document.body);
  assert.equal(backdrop.parentElement, document.body);
  assert.deepEqual([backdrop.style.zIndex, dialog.root.style.zIndex], ['3501', '3502']);
  assert.equal(Overlay.nextDialogZ(), 3503, 'what opens next — the platform\'s confirm included — goes above');
  dialog.dispose();
});

smoke('a press raises the lower of two dialogs, and Escape then goes to it', () => {
  const lower = Dialog.create('Lower').onCancel(() => {}).show({modal: true});
  const upper = Dialog.create('Upper').onCancel(() => {}).show({modal: true});
  assert.ok(Number(lower.root.style.zIndex) < Number(upper.root.style.zIndex));
  fire(lower.root, 'pointerdown');
  assert.ok(Number(lower.root.style.zIndex) > Number(upper.root.style.zIndex), 'raised');
  const z = lower.root.style.zIndex;
  fire(lower.root, 'pointerdown');
  assert.equal(lower.root.style.zIndex, z, 'the top one stays where it is');
  fire(document.body, 'keydown', {key: 'Escape'});
  assert.deepEqual([lower.isOpen.value, upper.isOpen.value], [false, true], 'the raised one is on top of the stack');
  upper.dispose();
});

smoke('dialog parts are stamped for automation', () => {
  const dialog = Dialog.create('Parts').add('Body').onOK(() => {});
  const part = (name) => dialog.root.querySelector(`[data-u2-part="${name}"]`);
  assert.equal(part('title').textContent, 'Parts');
  assert.equal(part('body').textContent, 'Body');
  assert.equal(part('footer').contains(dialog.getButton('OK')), true);
  assert.notEqual(part('close'), null);
  assert.notEqual(part('error-row'), null);
  dialog.dispose();
});

smoke('okDisabled: OK is off while the reason speaks, and says why in its tooltip; a form in the body adds its own', async () => {
  const picked = signal(0);
  const name = new TextInput({label: 'Name', required: true});
  const dialog = Dialog.create('New collection')
    .add(new Form().add(name))
    .okDisabled(computed(() => picked.value === 0 ? 'Select at least one monomer' : null))
    .onOK(() => {})
    .show();
  const ok = dialog.getButton('OK');
  const tip = async () => {
    fire(ok, 'pointerenter');
    await new Promise((r) => setTimeout(r, 700));
    const text = Tooltip.root.textContent;
    fire(ok, 'pointerleave');
    return text;
  };
  assert.equal(ok.disabled, true);
  assert.equal(await tip(), 'Name: Value can\'t be empty', 'the form\'s fields first, top to bottom');
  name.value.value = 'Peptides';
  assert.equal(ok.disabled, true);
  assert.equal(await tip(), 'Select at least one monomer');
  picked.value = 2;
  assert.equal(ok.disabled, false);
  dialog.close();
  name.dispose();
});
