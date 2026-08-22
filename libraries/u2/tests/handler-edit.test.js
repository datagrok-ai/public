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
import {Signal, signal} from '../src/core/signals.js';
import {Control} from '../src/core/component.js';
import {TextInput} from '../src/components/inputs/text-input.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {SpecEditor} from '../src/spec/editor.js';
import {backends} from '../src/sources/backends.js';
import {dfBindings} from '../src/sources/df-bindings.js';
import {DataFrame, Property, WidgetDescriptor, platform} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {SpecNodeRef} = await import('../src/dg/designer/node-ref.js');
const {SpecNodeHandler, propsFor, disposePanel} = await import('../src/dg/designer/handler.js');
const {TABLE_HINT} = await import('../src/dg/designer/prop-model.js');
const {shell} = await import('datagrok-api/grok');
const {registerPlatformComponents} = await import('../src/dg/viewers/registrations.js');

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
    const writes = shell.dart.writes.length;
    okButton(dialog).click();
    await flush();
    assert.deepEqual(patches.map((p) => [p.op, p.name, p.path]),
      [['set-bind', 'params.city', '$.reagent']]);
    assert.equal(node('orders').bind['params.city'], '$.reagent');
    assert.equal(field(shown, 'Bindings', 'params.city').value, '$.reagent',
      'the Bindings section shows the row at once, in place — the platform renders a shell.o write one gesture behind');
    assert.equal(shell.o?.node, node('orders'), 'the panel is asked for again: its shape changed');
    assert.equal(shell.dart.writes.length - writes, 1, 'one shell.o write per pick');

    const after = panel('orders');
    const city = field(after, 'Parameters', 'city');
    assert.equal(city.value, '$.reagent', 'the row shows what it follows');
    assert.equal(city.disabled, true, 'and the Bindings row below is where it is changed');
    assert.equal(field(after, 'Bindings', 'params.city').value, '$.reagent');
    assert.equal(instance.dump().components[0].bind['params.city'], '$.reagent');
  });

/* U1 of the viewers acceptance pass: OK on a group closed the picker and bound nothing — four
   attempts read as a broken feature. A group warns and the dialog stays open. */
edit('OK on a group warns and leaves the picker open; nothing reaches the document',
  async ({panel, patches, warnings}) => {
    const shown = panel('nameInput');
    fire(shown.querySelector('[data-u2-bind-pick="value"]'), 'click');
    const dialog = picker();
    // the fake source has no default step: its row is a group, as every heading is
    fire(pickerRow(dialog, 'orders'), 'click');
    okButton(dialog).click();
    await flush();
    assert.deepEqual(warnings, ['Pick a value — a group has nothing to bind to on its own']);
    assert.equal(patches.length, 0);
    assert.notEqual(document.querySelector('.u2-bind-picker'), null, 'the dialog stays up for the real pick');

    fire(pickerRow(dialog, 'reagent'), 'click');
    okButton(dialog).click();
    await flush();
    assert.deepEqual(patches.map((p) => [p.op, p.name, p.path]), [['set-bind', 'value', '$.reagent']]);
    assert.equal(document.querySelector('.u2-bind-picker'), null, 'and a leaf closes it');
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

/* WO-V5 — a property-tier node in the panel (VP-13): a set-prop re-creates a platform viewer, so
   its sections commit on blur/Enter, while a u2 input beside it keeps committing per keystroke.
   The viewer is the WO-V3 double, registered as `u2-viewer-scatter-plot` on a private registry. */
function tier(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const warn = console.warn;
    console.warn = () => {};
    // a descriptor with no label of its own answers the name as friendlyName and caption alike
    WidgetDescriptor.registry = [new WidgetDescriptor('Scatter plot', [
      new Property('xColumnName', 'string', {category: 'Data', friendlyName: 'xColumnName', caption: 'xColumnName'}),
      new Property('yColumnName', 'string', {category: 'Data', friendlyName: 'Y'}),
      new Property('showColumnNames', 'bool', {category: 'Data', friendlyName: 'Show Column Names',
        caption: 'Show Column Names'}),
      new Property('markerMinSize', 'num', {category: 'Misc'}),
      new Property('showRegressionLine', 'bool', {category: 'Regression'}),
    ])];
    const reg = new Registry();
    registerPlatformComponents(reg);
    const scope = new Scope();
    const df = new DataFrame([{name: 'city', type: 'string'}, {name: 'total', type: 'double'}],
      [{city: 'Kyiv', total: 1240}], 'orders');
    const spec = {$schema: 'dg-ui/1', root: {tag: 'u2-div-v', name: 'box', children: [
      {tag: 'u2-viewer-scatter-plot', name: 'plot', bind: {table: '$.orders'}, props: {xColumnName: 'total'}},
      {tag: 'u2-text-input', name: 'nameInput', props: {label: 'Name'}},
    ]}};
    const instance = renderSpec(spec,
      new SpecContext({data: {orders: dfBindings(signal(df), scope), reagent: 'total'}}), reg);
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
        instance, patches, editor,
        node: (n) => find(spec, n),
        panel: (n) => handler.renderProperties(new SpecNodeRef(instance, find(spec, n), editor)),
      });
    } finally {
      disposePanel();
      instance.dispose();
      scope.dispose();
      WidgetDescriptor.registry = [];
      platform.reset();
      console.warn = warn;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

tier('a property-tier node commits on change, never per keystroke; a u2 input keeps per-keystroke commits',
  ({instance, node, panel, patches}) => {
    const plot = node('plot');
    const shown = panel('plot');
    const before = instance.nodes().get(plot);
    assert.equal(before.propertyTier, true);
    const x = shown.querySelector('[data-u2-prop="xColumnName"] input');
    assert.equal(x.value, 'total', 'the live look value, read through the tier');
    for (const text of ['c', 'ci', 'cit'])
      fire(Object.assign(x, {value: text}), 'input');
    assert.equal(patches.length, 0, 'a keystroke builds no scatter plot');
    assert.equal(instance.nodes().get(plot), before);
    assert.equal(plot.props.xColumnName, 'total');

    fire(Object.assign(x, {value: 'city'}), 'change');
    assert.deepEqual(patches.map((p) => [p.op, p.name, p.value]), [['set-prop', 'xColumnName', 'city']]);
    assert.notEqual(instance.nodes().get(plot), before, 'one commit, one re-created viewer');
    assert.equal(instance.nodes().get(plot).props.xColumnName, 'city');
    assert.equal(x.value, 'city', 'the field refreshed in place from the new viewer');

    const label = panel('nameInput').querySelector('[data-u2-prop="label"] input');
    fire(Object.assign(label, {value: 'Reagent'}), 'input');
    assert.equal(patches.length, 2, 'a u2 input commits as it is typed');
    assert.equal(node('nameInput').props.label, 'Reagent');
  });

/* UX-B — the viewer's panel as the platform's own shows a viewer: the frame prop is a binding, not a
   field; the categories fold past the first, Misc/Description/Events last; the labels read as words;
   the Bindings section lists what is bound and one row to add the next. */

const label = (panel, title, name) =>
  section(panel, title).querySelector(`[data-u2-prop="${name}"] .u2-input-label`)?.textContent;

tier('a property-tier node: categories fold past the first, Misc last, labels humanized, the frame kept out of the fields',
  ({panel}) => {
    const shown = panel('plot');
    assert.deepEqual(panes(shown), [['Data', 'true'], ['Regression', 'false'], ['Misc', 'false'],
      ['Bindings', 'true'], ['Events', 'true']],
    'one pane per category, the first open, Misc behind the rest, the wiring sections after');
    assert.deepEqual(sections(shown), ['Node'], 'no plain section: every category folds');
    assert.equal(section(shown, 'Data').querySelector('[data-u2-prop="table"]'), null,
      'a dataframe prop is bind-only');
    assert.notEqual(section(shown, 'Misc').querySelector('[data-u2-prop="look"]'), null,
      'the look escape hatch files under Misc, behind the descriptor\'s own categories');
    assert.equal(label(shown, 'Data', 'xColumnName'), 'X Column Name',
      'camelCase read as words — the name echoed back as friendlyName/caption is no label');
    assert.equal(label(shown, 'Data', 'yColumnName'), 'Y', 'a friendlyName the descriptor gives wins');
    assert.equal(label(shown, 'Data', 'showColumnNames'), 'Show Column Names', 'as does a platform-authored caption');
    assert.equal(label(shown, 'Regression', 'showRegressionLine'), 'Show Regression Line');
    assert.equal(label(panel('nameInput'), 'Properties', 'label'), 'label', 'a u2 control keeps its prop names');
  });

tier('a property-tier node binds through one "Add binding…" row: the bound props, then the next one picked',
  async ({node, panel, patches}) => {
    const plot = node('plot');
    const shown = panel('plot');
    const rows = (p) => [...section(p, 'Bindings').querySelectorAll('[data-u2-prop]')].map((el) => el.dataset.u2Prop);
    assert.deepEqual(rows(shown), ['table', 'add-binding'], 'what is bound, and the row that adds');
    assert.equal(field(shown, 'Bindings', 'table').value, '$.orders');
    const select = section(shown, 'Bindings').querySelector('[data-u2-prop="add-binding"] select');
    const offered = select.querySelectorAll('option').map((o) => o.value).filter((v) => v !== '');
    assert.ok(offered.includes('xColumnName') && !offered.includes('table'), `the unbound props: ${offered}`);

    const said = [];
    const warning = shell.warning;
    shell.warning = (text) => said.push(text);
    try {
      fire(shown.querySelector('[data-u2-bind-pick="add-binding"]'), 'click');
      assert.equal(picker(), null, 'no property picked, no picker');
      assert.equal(said.length, 1);
      assert.match(said[0], /pick the property first/);
    } finally {
      shell.warning = warning;
    }

    select.value = 'xColumnName';
    fire(select, 'change');
    fire(shown.querySelector('[data-u2-bind-pick="add-binding"]'), 'click');
    const dialog = picker();
    assert.notEqual(dialog, null, 'the picker opens for the property named');
    fire(pickerRow(dialog, 'reagent'), 'click');
    const writes = shell.dart.writes.length;
    okButton(dialog).click();
    await flush();
    assert.deepEqual(patches.map((p) => [p.op, p.name, p.path]), [['set-bind', 'xColumnName', '$.reagent']]);
    assert.equal(plot.bind.xColumnName, '$.reagent');
    assert.deepEqual(rows(shown), ['table', 'xColumnName', 'add-binding'],
      'the row is there at once, in place — the platform renders a shell.o write one gesture behind');
    assert.equal(field(shown, 'Bindings', 'xColumnName').value, '$.reagent');
    const next = section(shown, 'Bindings').querySelector('[data-u2-prop="add-binding"] select');
    assert.equal([...next.querySelectorAll('option')].some((o) => o.value === 'xColumnName'), false,
      'and "Add binding…" no longer offers what is bound');
    assert.equal(shell.o?.node, plot, 'the panel is asked for again: it has a row more');
    assert.equal(shell.dart.writes.length - writes, 1, 'one shell.o write per pick');

    const after = panel('plot');
    assert.deepEqual(rows(after), ['table', 'xColumnName', 'add-binding']);
    assert.equal(field(after, 'Bindings', 'xColumnName').value, '$.reagent');
    assert.equal(field(after, 'Data', 'xColumnName').disabled, true, 'and the field is the binding\'s now');
  });

/** A viewer as the registry describes one, built by a create that wants its frame — the placeholder
 * every dropped viewer is until `table` is bound. */
const VIEWER = {
  tag: 'u2-e-viewer',
  category: 'Viewers',
  create: (props) => {
    if (props.table === undefined)
      throw new Error('u2-e-viewer needs a table');
    return new Control();
  },
  description: 'Fake viewer that builds only over a frame',
  props: [
    {name: 'table', type: 'dataframe', bindable: true},
    {name: 'xColumnName', type: 'string', bindable: true, category: 'Misc'},
    {name: 'yColumnName', type: 'string', bindable: true, category: 'Axes'},
    {name: 'title', type: 'string', bindable: true, category: 'Description'},
  ],
  example: {tag: 'u2-e-viewer'},
};

test('a broken viewer without its frame: the Error row says what to bind, Bindings opens on it, the frame is no field',
  async () => {
    const live = Scope.liveCount;
    const warn = console.warn;
    console.warn = () => {};
    const reg = new Registry();
    for (const meta of [BOX, VIEWER])
      reg.register(meta);
    const spec = {$schema: 'dg-ui/1', root: {tag: 'u2-e-box', name: 'layout', children: [
      {tag: 'u2-e-viewer', name: 'plot'}]}};
    const instance = renderSpec(spec, new SpecContext({data: {reagent: 'Ethanol'}}), reg);
    const editor = new SpecEditor(instance);
    const handler = new SpecNodeHandler();
    try {
      const ref = new SpecNodeRef(instance, find(spec, 'plot'), editor);
      assert.match(ref.error(), /needs a table/);
      const shown = handler.renderProperties(ref);
      const error = [...section(shown, 'Node').querySelectorAll('tr')]
        .find((tr) => tr.children[0].textContent === 'Error');
      assert.equal(error.querySelector('.u2-designer-hint')?.textContent, TABLE_HINT, 'the way out, under the error');
      assert.deepEqual(sections(shown), ['Node', 'Axes', 'Misc', 'Description'],
        'no frame field; Misc and Description last, the rest first-seen');
      assert.deepEqual(panes(shown), [['Bindings', 'true']], 'open on the binding that is missing');
      assert.deepEqual([...section(shown, 'Bindings').querySelectorAll('[data-u2-prop]')].map((el) => el.dataset.u2Prop),
        ['table', 'xColumnName', 'yColumnName', 'title'], 'a node that did not build is no tier: every bindable row');
    } finally {
      disposePanel();
      instance.dispose();
      console.warn = warn;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
