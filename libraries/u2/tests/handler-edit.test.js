/* WO-10/WO-2 — the writable context panel: the property model a node answers with, and what a
   field edit, a refusal, an undo and a binding change do to the document, the canvas and the panel
   itself. The panel is one propertyForm per section (F6/R1): inputs persist for the panel's
   lifetime and refresh in place — the same-element-identity test below is the invariant, and must
   fail if anyone reintroduces a rebuild. `DG.ObjectHandler` and `grok.shell` come from
   tests/dg-stub.mjs; the panel lives in the platform's context pane, so only its DOM is
   exercised here. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Signal} from '../src/core/signals.js';
import {Control} from '../src/core/component.js';
import {TextInput} from '../src/components/inputs/text-input.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {SpecEditor} from '../src/spec/editor.js';
import {backends} from '../src/sources/backends.js';
import {Property} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {SpecNodeRef} = await import('../src/dg/designer/node-ref.js');
const {SpecNodeHandler, propsFor, disposePanel} = await import('../src/dg/designer/handler.js');
const {shell} = await import('datagrok-api/grok');

/** A prop whose value lives nowhere but the component — what the dump-noise guard is about. */
class FakeInput extends TextInput {
  constructor(options) {
    super(options);
    this.hint = options.hint ?? 'Type here';
  }
}

const INPUT = {
  tag: 'u2-e-input',
  create: (props) => new FakeInput({
    label: props.label,
    value: props.value instanceof Signal ? undefined : props.value,
    bind: props.value instanceof Signal ? props.value : undefined,
    hint: props.hint,
  }),
  description: 'Fake input for the panel-editing tests',
  props: [
    {name: 'label', type: 'string', category: 'Appearance'},
    {name: 'value', type: 'string', bindable: true, twoWay: true},
    {name: 'stage', type: 'string', choices: ['Discovery', 'Phase I']},
    {name: 'sizes', type: 'object'},
    {name: 'items', type: 'string_list'},
    {name: 'hint', type: 'string'},
  ],
  events: ['input', 'change'],
  example: {tag: 'u2-e-input'},
};

const BOX = {
  tag: 'u2-e-box',
  create: () => new Control(),
  description: 'Fake container',
  props: [],
  acceptsChildren: true,
  example: {tag: 'u2-e-box'},
};

const TABS = {
  tag: 'u2-e-tabs',
  create: () => new Control(),
  createWithChildren: (_props, children, nodes) => {
    const strip = new Control();
    for (let i = 0; i < children.length; i++) {
      const label = document.createElement('span');
      label.className = 'fake-tab-label';
      label.textContent = nodes[i].props?.title ?? `Tab ${i + 1}`;
      strip.root.append(label, children[i] instanceof Control ? children[i].root : children[i]);
    }
    return strip;
  },
  description: 'Fake tab strip that consumes its children at construction',
  props: [],
  childProps: [{name: 'title', type: 'string', description: 'Tab label'}],
  acceptsChildren: true,
  example: {tag: 'u2-e-tabs'},
};

/** A source whose params are the function's own — what gives the panel a Parameters section. */
const SOURCE = {
  tag: 'u2-e-source',
  visual: false,
  createComponent: () => new Control(),
  description: 'Fake data source',
  props: [
    {name: 'func', type: 'string'},
    {name: 'params', type: 'object', subBindable: true},
  ],
  example: {tag: 'u2-e-source', name: 'orders'},
};

/** What `backends.funcRunner` answers for the source's function — the panel reads its inputs. */
const RUNNER = {
  find: (name) => name !== 'Orders' ? null : {name, kind: 'query', outputs: [],
    inputs: [new Property('days', 'int'), new Property('city', 'string')]},
  run: () => Promise.resolve({}),
};

/** A document of its own per test: patches mutate the spec in place. */
function editable() {
  return {$schema: 'dg-ui/1',
    components: [{tag: 'u2-e-source', name: 'orders', props: {func: 'Orders', params: {days: 30}}}],
    root: {tag: 'u2-e-box', name: 'layout', children: [
    {tag: 'u2-e-input', name: 'nameInput', props: {label: 'Name', value: 'Aspirin', sizes: [1, 2]},
      on: {input: 'cmd:touched'}},
    {tag: 'u2-e-input', name: 'boundInput', props: {label: 'Bound'}, bind: {value: '$.reagent'}},
    {tag: 'u2-e-input', name: 'brokenInput', props: {label: 'Ghost'}, bind: {value: '$.nowhere'}},
    {tag: 'u2-e-tabs', name: 'tabs', children: [
      {tag: 'u2-e-box', name: 'firstPane', props: {title: 'First'}},
    ]},
    {tag: 'div', name: 'block', props: {text: 'Hello', cls: 'greeting'}},
    {tag: 'a', name: 'docs', props: {text: 'Docs', href: '/help'}},
    {tag: 'img', name: 'logo', props: {src: '/logo.png'}},
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
  return [spec.root, ...spec.components ?? []].map(walk).find((node) => node !== null) ?? null;
}

/** The section under an `h3` — the form or the table the handler put there. Bindings and Events
 * are accordion panes instead: they fold, so what sits below them stays reachable. */
function section(panel, title) {
  const kids = [...panel.children];
  const i = kids.findIndex((el) => el.tagName === 'H3' && el.textContent === title);
  if (i >= 0)
    return kids[i + 1];
  const pane = [...panel.querySelectorAll('.u2-accordion-pane')]
    .find((p) => p.querySelector('.u2-accordion-title')?.textContent === title);
  return pane?.querySelector('.u2-accordion-content') ?? null;
}

/** The folding sections, with what each one is showing: `[title, aria-expanded]`. */
const panes = (panel) => [...panel.querySelectorAll('.u2-accordion-pane')]
  .map((p) => [p.querySelector('.u2-accordion-title').textContent,
    p.querySelector('.u2-accordion-header').getAttribute('aria-expanded')]);

function sections(panel) {
  return [...panel.children].filter((el) => el.tagName === 'H3').map((el) => el.textContent);
}

/** The editor of one panel field — every field root carries `data-u2-prop`, the automation id. */
function field(panel, title, name) {
  return section(panel, title).querySelector(`[data-u2-prop="${name}"]`)
    .querySelector('input, textarea, select');
}

/** Types into a panel field the way a user does — the form commits on every keystroke. */
function type(panel, title, name, text) {
  const input = field(panel, title, name);
  input.value = text;
  fire(input, 'input');
  return input;
}

/** Commits a guarded field (Bindings/Events) the way a user does — on blur or Enter. */
function commit(panel, title, name, text) {
  const input = field(panel, title, name);
  input.value = text;
  fire(input, 'change');
  return input;
}

/** The binding picker the `…` button opens, with a viewport handed to its virtual tree — the shim
 * measures nothing, and three overscan rows would hide everything below them. */
function picker() {
  const dialog = document.querySelector('.u2-bind-picker');
  const list = dialog?.querySelector('.u2-list');
  if (list) {
    list.clientHeight = 400;
    fire(list, 'scroll');
  }
  return dialog;
}

const pickerRows = (dialog) => [...dialog.querySelectorAll('.u2-tree-label')].map((el) => el.textContent);

const pickerRow = (dialog, label) => [...dialog.querySelectorAll('.u2-list-row')]
  .find((row) => (row.querySelector('.u2-tree-label')?.textContent ?? '').startsWith(label));

const okButton = (dialog) => [...dialog.querySelectorAll('.u2-dialog-buttons button')]
  .find((b) => b.textContent === 'OK');

function edit(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const warn = console.warn;
    const warnings = [];
    const infos = [];
    const warning = shell.warning;
    const info = shell.info;
    console.warn = () => {};
    shell.warning = (message) => warnings.push(message);
    shell.info = (message) => infos.push(message);
    const reg = new Registry();
    for (const meta of [INPUT, BOX, TABS, SOURCE])
      reg.register(meta);
    const spec = editable();
    const instance = renderSpec(spec, new SpecContext({data: {reagent: 'Ethanol'}}), reg);
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
        instance, editor, handler, patches, warnings, infos,
        node: (n) => find(spec, n),
        ref: (n) => new SpecNodeRef(instance, find(spec, n), editor),
        panel: (n) => handler.renderProperties(new SpecNodeRef(instance, find(spec, n), editor)),
      });
    } finally {
      disposePanel();
      instance.dispose();
      delete backends.funcRunner;
      shell.warning = warning;
      shell.info = info;
      console.warn = warn;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

edit('propsFor(): sections by category, editor routing, a bound prop and a structured value',
  ({ref}) => {
    const model = propsFor(ref('nameInput'));
    assert.deepEqual(model.map((s) => s.title), ['Appearance', 'Properties']);
    const [appearance, props] = model;
    assert.deepEqual(appearance.props.map((p) => p.name), ['label']);
    assert.deepEqual(props.props.map((p) => p.name), ['value', 'stage', 'sizes', 'items', 'hint']);
    const by = (name) => props.props.find((p) => p.name === name);
    assert.equal(typeof by('value').set, 'function');
    assert.deepEqual(by('stage').choices, ['Discovery', 'Phase I'],
      'a string constrained to choices routes to the choice editor');
    assert.equal(by('sizes').set, undefined, 'a structured value has no editor yet');
    assert.equal(by('items').type, 'list', 'string_list edits through the list editor');
    assert.equal(appearance.values.label, 'Name');
    assert.equal(props.values.value, 'Aspirin', 'the live value, read through getProperties()');
    assert.equal(props.values.sizes, '[1,2]', 'and the structured one as the JSON the spec carries');
    assert.equal(props.values.hint, 'Type here', 'a prop only the component knows about still reads');

    const bound = propsFor(ref('boundInput'));
    const value = bound.find((s) => s.title === 'Properties').props.find((p) => p.name === 'value');
    assert.equal(value.set, undefined,
      'a bound prop is the context\'s to change — the Bindings field is where it is edited');
    assert.equal(typeof bound.find((s) => s.title === 'Appearance').props[0].set, 'function');
  });

edit('propsFor(): a plain HTML node, a per-child prop, and a node that failed to build',
  ({ref}) => {
    const block = propsFor(ref('block'));
    assert.deepEqual(block.map((s) => s.title), ['Properties']);
    assert.deepEqual(block[0].props.map((p) => p.name), ['text', 'cls']);
    assert.deepEqual(block[0].values, {text: 'Hello', cls: 'greeting'});
    assert.deepEqual(propsFor(ref('docs'))[0].props.map((p) => p.name), ['text', 'cls', 'href']);
    assert.deepEqual(propsFor(ref('logo'))[0].props.map((p) => p.name), ['text', 'cls', 'src']);

    const pane = propsFor(ref('firstPane'));
    assert.deepEqual(pane.map((s) => s.title), ['Parent (u2-e-tabs)'],
      'what the parent reads off the child, under the parent\'s own section');
    assert.equal(pane[0].values.title, 'First');

    const broken = propsFor(ref('brokenInput'));
    assert.deepEqual(broken.flatMap((s) => s.props.map((p) => p.name)),
      INPUT.props.map((p) => p.name),
      'a node that failed to build still has the props that need fixing');
    assert.equal(broken.find((s) => s.title === 'Appearance').values.label, 'Ghost',
      'read from the node, since nothing was built');
  });

edit('a field edit is one patch: the document, the canvas and the dump all follow',
  ({instance, node, panel, patches}) => {
    const target = node('nameInput');
    const shown = panel('nameInput');
    const before = instance.nodes().get(target);
    type(shown, 'Appearance', 'label', 'Reagent');

    assert.equal(target.props.label, 'Reagent');
    assert.deepEqual(patches.map((p) => [p.op, p.name, p.value]), [['set-prop', 'label', 'Reagent']]);
    assert.notEqual(instance.nodes().get(target), before, 'the node was re-rendered in place');
    assert.equal(instance.nodes().get(target).root.querySelector('label').textContent, 'Reagent');
    assert.equal(instance.dump().root.children[0].props.label, 'Reagent');

    type(shown, 'Properties', 'hint', '');
    assert.equal(patches.length, 1, 'clearing a prop the spec never carried is dump noise, not an edit');
    assert.equal('hint' in target.props, false);
  });

edit('per-keystroke commits coalesce into a single undo entry', ({editor, node, panel}) => {
  const shown = panel('nameInput');
  const input = field(shown, 'Appearance', 'label');
  for (const text of ['R', 'Re', 'Rea', 'Reag']) {
    input.value = text;
    fire(input, 'input');
  }
  assert.equal(node('nameInput').props.label, 'Reag', 'every keystroke reached the document');

  editor.undo();
  assert.equal(node('nameInput').props.label, 'Name', 'one undo reverts the whole typed word');
  assert.equal(editor.canUndo.peek(), false, 'four keystrokes were one entry');
});

edit('a refused edit warns with the full reason and keeps the typed text for correction',
  async ({node, panel, patches, warnings}) => {
    const shown = panel('boundInput');
    const input = commit(shown, 'Bindings', 'value', 'reagent');
    await flush();

    assert.deepEqual(patches, [], 'nothing the engine refuses reaches the document');
    assert.equal(node('boundInput').bind.value, '$.reagent');
    assert.equal(input.value, 'reagent', 'the typed text stays on screen — the document was never touched');
    assert.equal(warnings.length, 1);
    assert.match(warnings[0], /must start with "\$\."/);

    commit(shown, 'Bindings', 'value', '$.other');
    assert.equal(node('boundInput').bind.value, '$.other', 'correcting the kept text commits');
    assert.equal(warnings.length, 1);
  });

edit('a bound prop and a structured value render read-only', ({panel}) => {
  const bound = panel('boundInput');
  const value = field(bound, 'Properties', 'value');
  assert.equal(value.disabled, true, 'the binding field is where a bound prop is edited');
  assert.equal(value.value, 'Ethanol', 'still showing the live bound value');
  assert.equal(field(bound, 'Appearance', 'label').disabled, false);

  const named = panel('nameInput');
  const sizes = field(named, 'Properties', 'sizes');
  assert.equal(sizes.disabled, true);
  assert.equal(sizes.value, '[1,2]', 'a structured value shows as the JSON the spec carries');
});

edit('a choice prop edits through a choice list, a string_list through the list editor',
  ({node, panel, patches}) => {
    const shown = panel('nameInput');
    const stage = field(shown, 'Properties', 'stage');
    assert.equal(stage.tagName, 'SELECT');
    stage.value = 'Phase I';
    fire(stage, 'change');
    assert.equal(node('nameInput').props.stage, 'Phase I');

    const items = field(shown, 'Properties', 'items');
    items.value = 'One,Two';
    fire(items, 'input');
    assert.deepEqual(node('nameInput').props.items, ['One', 'Two'], 'committed as the parsed list');
    assert.deepEqual(patches.map((p) => [p.op, p.name]), [['set-prop', 'stage'], ['set-prop', 'items']]);
  });

edit('without an editor the panel is read-only — and an HTML node still shows and edits its props',
  ({instance, node, panel, patches}) => {
    const readOnly = new SpecNodeHandler().renderProperties(new SpecNodeRef(instance, node('block')));
    assert.deepEqual(sections(readOnly), ['Node', 'Properties']);
    assert.equal(readOnly.querySelectorAll('input').length, 0, 'nothing to type into');
    assert.match(readOnly.textContent, /Hello/);
    assert.match(readOnly.textContent, /greeting/);

    const writable = panel('block');
    assert.deepEqual(sections(writable), ['Node', 'Properties']);
    assert.equal(writable.querySelectorAll('input').length, 2, 'the same model, with editors');
    type(writable, 'Properties', 'text', 'Edited');
    assert.equal(node('block').props.text, 'Edited');
    assert.deepEqual(patches.map((p) => p.op), ['set-prop']);
  });

edit('an undo refreshes the field in place — the input the user is typing in is never rebuilt',
  ({editor, node, panel}) => {
    const shown = panel('nameInput');
    const input = type(shown, 'Appearance', 'label', 'Reagent');
    assert.equal(node('nameInput').props.label, 'Reagent');

    editor.undo();
    assert.equal(node('nameInput').props.label, 'Name');
    assert.equal(field(shown, 'Appearance', 'label'), input, 'the very same input element');
    assert.equal(input.value, 'Name', 'showing what the document says now');

    editor.redo();
    assert.equal(input.value, 'Reagent');
  });

edit('bindings and events edit through the panel, and clearing one drops it from the dump',
  ({instance, node, panel, patches, warnings, infos}) => {
    const bound = panel('boundInput');
    commit(bound, 'Bindings', 'value', '$.other');
    assert.equal(node('boundInput').bind.value, '$.other');
    assert.deepEqual(infos, [], 'a rewire is visible on the canvas — no balloon');

    commit(bound, 'Bindings', 'value', '');
    assert.equal(node('boundInput').bind, undefined);
    assert.equal('bind' in instance.dump().root.children[1], false, 'and the dump carries no empty container');
    assert.deepEqual(infos, ['value: binding removed'], 'clearing wiring is deliberate but never silent');

    const named = panel('nameInput');
    assert.deepEqual([...section(named, 'Events').querySelectorAll('[data-u2-prop]')]
      .map((el) => el.dataset.u2Prop), ['input', 'change'], 'every declared event, wired or not');
    commit(named, 'Events', 'change', 'save');
    assert.equal(warnings[warnings.length - 1],
      'an event must name a command as \'cmd:\' followed by the command name — got \'save\'');
    assert.equal(node('nameInput').on.change, undefined);

    commit(named, 'Events', 'change', 'cmd:save');
    assert.equal(node('nameInput').on.change, 'cmd:save');
    assert.deepEqual(patches.map((p) => p.op), ['set-bind', 'set-bind', 'set-event']);
  });

edit('a guarded field commits on change, never per keystroke — partial values neither warn nor wipe',
  async ({node, panel, patches, warnings}) => {
    const shown = panel('nameInput');
    const input = field(shown, 'Events', 'input');
    for (const text of ['c', 'cm', 'cmd', 'cmd:', 'cmd:p', 'cmd:pi']) {
      input.value = text;
      fire(input, 'input');
    }
    await flush();
    assert.deepEqual(patches, [], 'no keystroke reached the engine');
    assert.equal(warnings.length, 0, 'so nothing was refused mid-word');
    assert.equal(input.value, 'cmd:pi', 'and the in-progress text was never wiped');

    input.value = 'cmd:ping';
    fire(input, 'change');
    assert.equal(node('nameInput').on.input, 'cmd:ping', 'the commit lands on blur/Enter');
    assert.deepEqual(patches.map((p) => [p.op, p.command]), [['set-event', 'cmd:ping']]);

    // select-all-retype (the acceptance corruption): under per-keystroke commits the first
    // keystroke's refusal-revert restored the old text under the caret and the rest APPENDED —
    // 'cmd:saveave' silently rewired the button; keystrokes must not commit and the refused
    // final value must stay put
    for (const text of ['s', 'sa', 'sav', 'save']) {
      input.value = text;
      fire(input, 'input');
    }
    await flush();
    assert.equal(input.value, 'save', 'nothing restored old text under the caret');
    fire(input, 'change');
    assert.equal(node('nameInput').on.input, 'cmd:ping', 'the refused retype never reached the document');
    assert.equal(input.value, 'save', 'and the typed text stays for correction');
    assert.equal(warnings.length, 1);
    assert.deepEqual(patches.map((p) => [p.op, p.command]), [['set-event', 'cmd:ping']]);
  });

edit('a per-child prop is the parent\'s to render: editing a pane title rebuilds the tab strip',
  ({instance, node, panel}) => {
    const tabs = node('tabs');
    const before = instance.nodes().get(tabs);
    const shown = panel('firstPane');
    assert.equal(field(shown, 'Parent (u2-e-tabs)', 'title').value, 'First');

    type(shown, 'Parent (u2-e-tabs)', 'title', 'Renamed');
    assert.equal(node('firstPane').props.title, 'Renamed');
    assert.notEqual(instance.nodes().get(tabs), before, 'the parent that consumes the title re-rendered');
    assert.equal(instance.root.querySelector('.fake-tab-label').textContent, 'Renamed');
  });

/* WO-15 — the Bindings section as the place a binding is ADDED, not only corrected: one row per
   bindable prop, and the picker beside the path field for whoever does not want to write the
   grammar by hand. */

edit('the Bindings section lists every bindable prop, bound or not — an empty row commits nothing',
  ({node, panel, patches}) => {
    const shown = panel('nameInput');
    assert.deepEqual([...section(shown, 'Bindings').querySelectorAll('[data-u2-prop]')]
      .map((el) => el.dataset.u2Prop), ['value'], 'the tag\'s bindable prop, unbound');
    assert.equal(field(shown, 'Bindings', 'value').value, '');

    commit(shown, 'Bindings', 'value', '');
    assert.deepEqual(patches, [], 'an untouched empty row is not an edit');
    assert.equal(node('nameInput').bind, undefined);

    assert.equal(field(panel('boundInput'), 'Bindings', 'value').value, '$.reagent',
      'and a bound row shows the path the document carries');
  });

edit('the picker commits the path it assembles, through the funnel a typed path takes',
  async ({node, panel, patches}) => {
    const shown = panel('nameInput');
    fire(shown.querySelector('[data-u2-bind-pick="value"]'), 'click');
    const dialog = picker();
    assert.notEqual(dialog, null, 'the … button opens the picker');
    assert.deepEqual(pickerRows(dialog).slice(0, 2), ['App data', 'reagent : string ⇄'],
      'the roots sit under the heading that says where they came from');

    fire(pickerRow(dialog, 'reagent'), 'click');
    okButton(dialog).click();
    await flush();

    assert.deepEqual(patches.map((p) => [p.op, p.name, p.path]), [['set-bind', 'value', '$.reagent']]);
    assert.equal(node('nameInput').bind.value, '$.reagent');
    assert.equal(field(shown, 'Bindings', 'value').value, '$.reagent', 'and the field follows the document');
    assert.equal(document.querySelector('.u2-bind-picker'), null, 'the dialog closes with its OK');
  });

/* WO-16 — the Parameters section of a source that names a function: one typed editor per declared
   input, edits coalescing into the one `params` prop, and the picker writing a dotted sub-bind. */

edit('the function\'s inputs edit as one coalescing params prop', ({editor, node, panel, patches}) => {
  backends.funcRunner = RUNNER;
  const shown = panel('orders');
  assert.deepEqual(sections(shown), ['Node', 'Properties', 'Parameters']);
  assert.deepEqual([...section(shown, 'Parameters').querySelectorAll('[data-u2-prop]')]
    .map((el) => el.dataset.u2Prop), ['days', 'city'], 'a row per input the function declares');
  assert.equal(field(shown, 'Parameters', 'days').value, '30', 'showing what the node passes');

  type(shown, 'Parameters', 'days', '45');
  type(shown, 'Parameters', 'city', 'Kyiv');
  assert.deepEqual(node('orders').props.params, {days: 45, city: 'Kyiv'});
  assert.deepEqual(patches.map((p) => [p.op, p.name]), [['set-prop', 'params'], ['set-prop', 'params']],
    'every param edits the one prop the source declares');

  editor.undo();
  assert.deepEqual(node('orders').props.params, {days: 30},
    'and edits to it coalesce: one undo takes the whole pass back');
});

edit('a param bind is a dotted set-bind, and the row it lands on becomes the binding\'s',
  async ({instance, node, panel, patches}) => {
    backends.funcRunner = RUNNER;
    const shown = panel('orders');
    fire(shown.querySelector('[data-u2-bind-pick="params.city"]'), 'click');
    const dialog = picker();
    assert.notEqual(dialog, null, 'every param row carries the picker');

    fire(pickerRow(dialog, 'reagent'), 'click');
    okButton(dialog).click();
    await flush();
    assert.deepEqual(patches.map((p) => [p.op, p.name, p.path]),
      [['set-bind', 'params.city', '$.reagent']]);
    assert.equal(node('orders').bind['params.city'], '$.reagent');
    assert.equal(shell.o?.node, node('orders'), 'the panel is asked for again: its shape changed');

    const after = panel('orders');
    const city = field(after, 'Parameters', 'city');
    assert.equal(city.value, '$.reagent', 'the row shows what it follows');
    assert.equal(city.disabled, true, 'and the Bindings row below is where it is changed');
    assert.equal(field(after, 'Bindings', 'params.city').value, '$.reagent');
    assert.equal(instance.dump().components[0].bind['params.city'], '$.reagent');
  });

edit('a cycle picked in the tree is refused: the loop is named and nothing reaches the document',
  async ({node, panel, patches, warnings}) => {
    const shown = panel('nameInput');
    fire(shown.querySelector('[data-u2-bind-pick="value"]'), 'click');
    const dialog = picker();
    fire(pickerRow(dialog, 'nameInput').querySelector('.u2-tree-twistie'), 'click');
    await flush();
    fire(pickerRow(dialog, 'value'), 'click');
    okButton(dialog).click();
    await flush();

    assert.deepEqual(patches, [], 'a self-reference never becomes a patch');
    assert.equal(node('nameInput').bind, undefined);
    assert.equal(field(shown, 'Bindings', 'value').value, '', 'and the row is still empty');
    assert.equal(warnings.length, 1);
    assert.match(warnings[0], /would create a cycle: nameInput → nameInput/);
  });

/* P4 acceptance — the panel was taller than the platform's pane, which scrolls with no scrollbar
   and no fold shadow: the Events section simply was not there. The wiring sections fold. */
edit('Bindings folds until something is bound; Events, the section it was hiding, stays open',
  ({panel}) => {
    const unbound = panel('nameInput');
    assert.deepEqual(panes(unbound), [['Bindings', 'false'], ['Events', 'true']]);
    assert.notEqual(field(unbound, 'Bindings', 'value'), null,
      'folded, never unbuilt — the fields refresh in place on every patch, opened or not');
    assert.deepEqual(panes(panel('boundInput')), [['Bindings', 'true'], ['Events', 'true']],
      'wiring that exists is shown');
  });
