/* The automation identity contract: data-u2 (kind) + data-u2-name (instance) + data-u2-part
   (structure) + data-u2-owner (portaled popups), and the naming nudge. Doc:
   core/docs/features/ui2/AUTOMATION.md. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, Control, Combobox, TextInput, ChoiceInput, Form, Dialog, Overlay,
  div, renderSpec, Registry, registerAll} from '../src/index.js';

function check(name, body) {
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

check('Control: setting name stamps data-u2-name; clearing removes it', () => {
  const control = new Control();
  assert.equal(control.root.dataset.u2Name, undefined, 'unnamed roots carry no automation id');

  control.name = 'org';
  assert.equal(control.root.dataset.u2Name, 'org');
  control.name = 'organization';
  assert.equal(control.root.dataset.u2Name, 'organization', 'a rename re-stamps');
  control.name = undefined;
  assert.equal(control.root.dataset.u2Name, undefined, 'clearing the name clears the id');
  control.dispose();
});

check('Control: a reassigned root carries the name over', () => {
  const control = new Control();
  control.name = 'org';
  const replacement = document.createElement('div');
  control.root = replacement;
  assert.equal(replacement.dataset.u2Name, 'org', 'the new root is stamped');
  control.dispose();
});

check('Input: an explicit name is the automation id, a label alone is only the form key', () => {
  const named = new TextInput({label: 'Organization', name: 'org'});
  assert.equal(named.name, 'org');
  assert.equal(named.root.dataset.u2Name, 'org');

  const labeled = new TextInput({label: 'Organization'});
  assert.equal(labeled.name, 'Organization', 'the label fallback still keys forms');
  assert.equal(labeled.root.dataset.u2Name, undefined,
    'a caption is not a stable identity and never stamps');

  named.dispose();
  labeled.dispose();
});

check('Input: parts are stamped and getWidgetStatus().parts answers the same elements', () => {
  const input = new TextInput({label: 'Organization', name: 'org'});
  const parts = input.getWidgetStatus().parts;
  assert.equal(parts.editor.dataset.u2Part, 'editor');
  assert.equal(parts.label.dataset.u2Part, 'label');
  assert.equal(parts.error.dataset.u2Part, 'error');
  assert.equal(input.root.querySelector('[data-u2-part="editor"]'), parts.editor,
    'the attribute and the status surface point at one element');
  assert.equal(parts.options, undefined, 'an empty rail is not a part');

  const withPostfix = new TextInput({label: 'Dose', name: 'dose', postfix: 'mg'});
  assert.equal(withPostfix.getWidgetStatus().parts.options.dataset.u2Part, 'options');

  const inline = new TextInput({inline: true, name: 'filter'});
  assert.equal(inline.getWidgetStatus().parts.label, undefined, 'inline inputs have no label part');

  input.dispose();
  withPostfix.dispose();
  inline.dispose();
});

check('the hierarchy composes: dialog › form › input resolves as one descendant selector', () => {
  const org = new TextInput({label: 'Organization', name: 'org'});
  const form = new Form();
  form.name = 'projectForm';
  form.add(org);
  const dialog = Dialog.create('Load Project', {name: 'loadProjectDialog'}).add(form);
  dialog.show();

  const hit = document.querySelector(
    '[data-u2-name="loadProjectDialog"] [data-u2-name="projectForm"] [data-u2-name="org"]');
  assert.equal(hit, org.root, 'the load-project selector lands on the input');
  assert.equal(document.querySelector('[data-u2-name="loadProjectDialog"]'), dialog.root);

  dialog.dispose();
  form.dispose();
  org.dispose();
});

check('Overlay: portaled content carries data-u2-owner from the nearest named ancestor', async () => {
  const scope = new Scope();
  const owner = div([], 'host');
  owner.dataset.u2Name = 'org';
  const anchor = document.createElement('button');
  owner.append(anchor);
  document.body.append(owner);

  const content = document.createElement('div');
  const close = Overlay.show(anchor, content, scope);
  await flush();
  assert.equal(content.dataset.u2Owner, 'org', 'the popup names its owner across the portal');
  close();

  // reused content from an unnamed anchor must not keep a stale owner
  const bare = document.createElement('button');
  document.body.append(bare);
  const reclose = Overlay.show(bare, content, scope);
  await flush();
  assert.equal(content.dataset.u2Owner, undefined, 'no named ancestor, no owner link');
  reclose();
  scope.dispose();
});

check('Combobox: the dropdown of a named control is addressable through its owner', async () => {
  const combobox = new Combobox({items: ['alpha', 'beta'], placeholder: 'Pick'});
  combobox.name = 'org';
  document.body.append(combobox.root);
  const input = combobox.root.querySelector('input');
  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();

  const popup = document.body.querySelector('.u2-combobox-popup');
  assert.equal(popup.dataset.u2Owner, 'org',
    'the popup is not a DOM descendant, but the owner edge recovers the hierarchy');
  combobox.dispose();
});

check('spec-built and hand-built UI share the selector vocabulary', () => {
  const reg = new Registry();
  registerAll(reg);
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-form', name: 'projectForm',
    children: [{tag: 'u2-text-input', name: 'org', props: {label: 'Organization'}}]}},
  undefined, reg);
  document.body.append(instance.root);
  const hit = document.querySelector('[data-u2-name="projectForm"] [data-u2-name="org"]');
  assert.equal(hit, Control.forElement(hit).root, 'the attribute and forElement agree');
  instance.dispose();
});

check('Form: an input with neither name nor label warns once per kind', () => {
  const warnings = [];
  const original = console.warn;
  console.warn = (message) => warnings.push(String(message));
  const form = new Form();
  const first = new ChoiceInput({items: ['a']});
  const second = new ChoiceInput({items: ['b']});
  const labeled = new TextInput({label: 'Named'});
  try {
    form.add(first);
    form.add(second);
    const nudges = warnings.filter((w) => w.includes('no form key and no automation id'));
    assert.equal(nudges.length, 1, 'one nudge per kind, not per instance');
    assert.ok(nudges[0].includes('choice-input'), `the nudge names the kind: ${nudges[0]}`);

    form.add(labeled);
    assert.equal(warnings.filter((w) => w.includes('no form key')).length, 1,
      'a labeled input is keyed and stays quiet');
  } finally {
    console.warn = original;
    form.dispose();
    first.dispose();
    second.dispose();
    labeled.dispose();
  }
});
