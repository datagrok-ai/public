/* WO-3 (Appearance) — the designer panel over the shared appearance group: section ranking, the
   folded pane with its count badge, clear-to-delete through the panel path (Ruling 8), the
   component-own-prop merge (R-a), the bound-appearance Bindings rows and the generalized
   "Add binding…". Pattern of tests/handler-edit.test.js; dg-stub.mjs serves the platform. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {TextInput} from '../src/components/inputs/text-input.js';
import {Registry} from '../src/spec/registry.js';
import {registerAll} from '../src/spec/registrations.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {SpecEditor} from '../src/spec/editor.js';

register('./dg-stub.mjs', import.meta.url);
const {SpecNodeRef} = await import('../src/dg/designer/node-ref.js');
const {SpecNodeHandler, propsFor, disposePanel} = await import('../src/dg/designer/handler.js');
const {bindsOf} = await import('../src/dg/designer/prop-model.js');
const {shell} = await import('datagrok-api/grok');

/** R-a: a component-own prop may declare the shared category and merges into the section. */
const OWN = {
  tag: 'u2-a-own',
  create: (props) => new TextInput({label: props.label}),
  description: 'Own prop declared into the shared Appearance category',
  props: [
    {name: 'label', type: 'string'},
    {name: 'accent', type: 'string', category: 'Appearance'},
  ],
  example: {tag: 'u2-a-own'},
};

const MISC = {
  tag: 'u2-a-misc',
  create: () => new TextInput({}),
  description: 'A Misc-category prop for the ranking test',
  props: [
    {name: 'label', type: 'string'},
    {name: 'debug', type: 'string', category: 'Misc'},
  ],
  example: {tag: 'u2-a-misc'},
};

function editable() {
  return {$schema: 'dg-ui/1', root: {tag: 'div', name: 'layout', children: [
    {tag: 'u2-text-input', name: 'plain', props: {label: 'Name'}},
    {tag: 'u2-text-input', name: 'styled', props: {label: 'Styled', color: '#ff0000', padding: '4px'}},
    {tag: 'u2-text-input', name: 'boundColor', props: {label: 'Bound'}, bind: {color: '$.c'}},
    {tag: 'u2-a-own', name: 'own', props: {accent: 'x'}},
    {tag: 'u2-a-misc', name: 'misc'},
  ]}};
}

function find(spec, name) {
  const walk = (n) => {
    if (n.name === name)
      return n;
    for (const child of n.children ?? []) {
      const found = walk(child);
      if (found)
        return found;
    }
    return null;
  };
  return walk(spec.root);
}

/** The section under an `h3`, or a folding pane's content — matched by title prefix, because the
 * Appearance pane's title carries the count badge. */
function section(panel, title) {
  const kids = [...panel.children];
  const i = kids.findIndex((el) => el.tagName === 'H3' && el.textContent === title);
  if (i >= 0)
    return kids[i + 1];
  const pane = [...panel.querySelectorAll('.u2-accordion-pane')]
    .find((p) => p.querySelector('.u2-accordion-title')?.textContent.startsWith(title));
  return pane?.querySelector('.u2-accordion-content') ?? null;
}

/** The Appearance pane as `{title, expanded}` — what the badge and the fold say. */
function appearancePane(panel) {
  const pane = [...panel.querySelectorAll('.u2-accordion-pane')]
    .find((p) => p.querySelector('.u2-accordion-title').textContent.startsWith('Appearance'));
  return pane === undefined ? null : {
    title: pane.querySelector('.u2-accordion-title').textContent,
    expanded: pane.querySelector('.u2-accordion-header').getAttribute('aria-expanded'),
  };
}

function field(panel, title, name) {
  return section(panel, title).querySelector(`[data-u2-prop="${name}"]`)
    .querySelector('input, textarea, select');
}

function type(panel, title, name, text) {
  const input = field(panel, title, name);
  input.value = text;
  fire(input, 'input');
  return input;
}

function picker() {
  const dialog = document.querySelector('.u2-bind-picker');
  const list = dialog?.querySelector('.u2-list');
  if (list) {
    list.clientHeight = 400;
    fire(list, 'scroll');
  }
  return dialog;
}

const pickerRow = (dialog, label) => [...dialog.querySelectorAll('.u2-list-row')]
  .find((row) => (row.querySelector('.u2-tree-label')?.textContent ?? '').startsWith(label));

const okButton = (dialog) => [...dialog.querySelectorAll('.u2-dialog-buttons button')]
  .find((b) => b.textContent === 'OK');

function edit(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const warn = console.warn;
    const warnings = [];
    const warning = shell.warning;
    console.warn = () => {};
    shell.warning = (message) => warnings.push(message);
    const reg = new Registry();
    registerAll(reg);
    for (const meta of [OWN, MISC])
      reg.register(meta);
    const spec = editable();
    const instance = renderSpec(spec, new SpecContext({data: {c: 'red'}}), reg);
    const editor = new SpecEditor(instance);
    const patches = [];
    const apply = editor.apply.bind(editor);
    editor.apply = (patch) => {
      patches.push(patch);
      apply(patch);
    };
    const handler = new SpecNodeHandler();
    try {
      await body({
        instance, editor, patches, warnings,
        node: (n) => find(spec, n),
        ref: (n) => new SpecNodeRef(instance, find(spec, n), editor),
        panel: (n) => handler.renderProperties(new SpecNodeRef(instance, find(spec, n), editor)),
      });
    } finally {
      disposePanel();
      instance.dispose();
      shell.warning = warning;
      console.warn = warn;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

edit('propsFor ranks Appearance after the first-seen sections, before Misc', ({ref}) => {
  assert.deepEqual(propsFor(ref('plain')).map((s) => s.title), ['Properties', 'Appearance']);
  assert.deepEqual(propsFor(ref('misc')).map((s) => s.title), ['Properties', 'Appearance', 'Misc']);
});

edit('the pane folds on an unstyled node and opens with the count on a styled one', ({panel}) => {
  assert.deepEqual(appearancePane(panel('plain')), {title: 'Appearance', expanded: 'false'});
  assert.deepEqual(appearancePane(panel('styled')), {title: 'Appearance (2)', expanded: 'true'});
});

edit('the badge tracks set and clear in place; the fold state is never forced', ({node, panel}) => {
  const shown = panel('plain');
  type(shown, 'Appearance', 'padding', '8px');
  assert.equal(node('plain').props.padding, '8px');
  assert.deepEqual(appearancePane(shown), {title: 'Appearance (1)', expanded: 'false'});
  type(shown, 'Appearance', 'padding', '');
  assert.equal('padding' in node('plain').props, false);
  assert.deepEqual(appearancePane(shown), {title: 'Appearance', expanded: 'false'});
});

edit('Ruling 8: set then clear through the panel leaves the dump byte-identical (text, choice, number)',
  ({instance, node, panel}) => {
    const before = JSON.stringify(instance.dump());
    const shown = panel('plain');

    type(shown, 'Appearance', 'padding', '8px 12px');
    assert.equal(node('plain').props.padding, '8px 12px');
    type(shown, 'Appearance', 'padding', '');
    assert.equal(JSON.stringify(instance.dump()), before, 'text: no trace');

    const cursor = field(shown, 'Appearance', 'cursor');
    assert.equal(cursor.tagName, 'SELECT');
    cursor.value = 'pointer';
    fire(cursor, 'change');
    assert.equal(node('plain').props.cursor, 'pointer');
    cursor.value = '';
    fire(cursor, 'change');
    assert.equal(JSON.stringify(instance.dump()), before, 'choice: cleared to absent, not null');

    type(shown, 'Appearance', 'opacity', '0.5');
    assert.equal(node('plain').props.opacity, 0.5);
    type(shown, 'Appearance', 'opacity', '');
    assert.equal(JSON.stringify(instance.dump()), before, 'number: no trace');
  });

edit('undo of a clear restores the value, the style and the dump', ({editor, instance, node, panel}) => {
  const shown = panel('plain');
  type(shown, 'Appearance', 'padding', '8px');
  // a different prop between set and clear, so the clear is its own undo entry
  type(shown, 'Appearance', 'opacity', '0.5');
  type(shown, 'Appearance', 'padding', '');
  assert.equal('padding' in node('plain').props, false);

  editor.undo();
  assert.equal(node('plain').props.padding, '8px');
  assert.equal(instance.nodes().get(node('plain')).root.style.padding, '8px');
  assert.equal(field(shown, 'Appearance', 'padding').value, '8px', 'the field refreshed in place');
});

edit('the color field styles the node; typed-empty is refused, so \'\' can never serialize',
  ({instance, node, panel}) => {
    const shown = panel('plain');
    const color = field(shown, 'Appearance', 'color');
    color.value = '#ff0000';
    fire(color, 'input');
    assert.equal(node('plain').props.color, '#ff0000');
    const styled = JSON.stringify(instance.dump());

    // the hex-only ColorInput keeps '' out of the value signal (color-input.ts HEX gate): the
    // clear gesture is missing for this kind, but the Ruling 8 half it does guarantee holds
    color.value = '';
    fire(color, 'input');
    assert.equal(JSON.stringify(instance.dump()), styled, 'nothing written, nothing deleted');
    assert.equal(node('plain').props.color, '#ff0000');
  });

edit('R-a: a component-own prop in the shared category opens the pane; the badge counts the shared group only',
  ({node, panel}) => {
    const shown = panel('own');
    assert.deepEqual(appearancePane(shown), {title: 'Appearance', expanded: 'true'},
      'expanded for the own prop, which the badge does not count');
    type(shown, 'Appearance', 'padding', '4px');
    assert.deepEqual(appearancePane(shown), {title: 'Appearance (1)', expanded: 'true'});
    assert.equal(node('own').props.accent, 'x', 'the own prop is untouched');

    type(shown, 'Appearance', 'accent', '');
    assert.equal(node('own').props.accent, '',
      'clearing the own prop writes \'\' like any category — delete-on-empty is the shared group\'s only');
  });

edit('the read-only panel hides unassigned shared appearance props and empty sections',
  ({instance, node}) => {
    const handler = new SpecNodeHandler();
    const readOnly = (n) => handler.renderProperties(new SpecNodeRef(instance, node(n)));
    const rows = (shown, title) => {
      const el = section(shown, title);
      return el === null ? null : [...el.querySelectorAll('tr')].map((tr) => tr.children[0].textContent);
    };
    assert.equal(section(readOnly('plain'), 'Appearance'), null,
      'an unstyled node shows no Appearance section at all');
    const styled = readOnly('styled');
    assert.deepEqual(rows(styled, 'Appearance'), ['color', 'padding'], 'only the assigned rows');
    assert.ok(rows(styled, 'Properties').includes('label'), 'the ordinary sections keep every row');
  });

edit('bindsOf hides the unbound appearance rows; a bound one shows and its field goes read-only',
  ({ref, panel}) => {
    assert.deepEqual(Object.keys(bindsOf(ref('plain'))), ['value']);
    const bound = bindsOf(ref('boundColor'));
    assert.deepEqual(Object.keys(bound), ['value', 'color']);
    assert.equal(bound.color, '$.c');

    const shown = panel('boundColor');
    assert.deepEqual([...section(shown, 'Bindings').querySelectorAll('[data-u2-prop]')]
      .map((el) => el.dataset.u2Prop), ['value', 'color', 'add-binding']);
    assert.equal(field(shown, 'Bindings', 'color').value, '$.c');
    assert.equal(field(shown, 'Appearance', 'color').disabled, true,
      'the Bindings row is where a bound prop is edited');

    const select = section(shown, 'Bindings').querySelector('[data-u2-prop="add-binding"] select');
    const offered = [...select.querySelectorAll('option')].map((o) => o.value).filter((v) => v !== '');
    assert.ok(offered.includes('padding') && !offered.includes('value') && !offered.includes('color'),
      `appearance props only, minus the bound one: ${offered}`);
  });

edit('an "Add binding…" pick on a plain node creates the appearance bind through the picker',
  async ({node, panel, patches}) => {
    const shown = panel('plain');
    assert.deepEqual(appearancePane(shown)?.expanded, 'false');
    const select = section(shown, 'Bindings').querySelector('[data-u2-prop="add-binding"] select');
    select.value = 'color';
    fire(select, 'change');
    fire(shown.querySelector('[data-u2-bind-pick="add-binding"]'), 'click');
    const dialog = picker();
    assert.notEqual(dialog, null, 'the … button opens the picker');
    fire(pickerRow(dialog, 'c'), 'click');
    okButton(dialog).click();
    await flush();

    assert.deepEqual(patches.map((p) => [p.op, p.name, p.path]), [['set-bind', 'color', '$.c']]);
    assert.equal(node('plain').bind.color, '$.c');
    const after = panel('plain');
    assert.equal(field(after, 'Bindings', 'color').value, '$.c');
    assert.equal(field(after, 'Appearance', 'color').disabled, true);
  });
