/* The read-only inspector (P2): the models a designer derives from a rendered spec — breadcrumb
   path, structure tree, broken-node count — and what the node ObjectHandler renders from them.
   `DG.ObjectHandler` comes from tests/dg-stub.mjs. The glass pane hit-test and the shell view need
   a real browser and a running platform; they are the live pass, not this file. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {Control} from '../src/core/component.js';
import {TextInput} from '../src/components/inputs/text-input.js';
import {TabStrip} from '../src/components/containers/tabs.js';
import {Accordion} from '../src/components/containers/accordion.js';
import {Registry, registry as globalRegistry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {SAMPLES} from '../src/dg/designer/samples.js';
import {dfBindings} from '../src/sources/df-bindings.js';
import {DataFrame, Property, WidgetDescriptor, platform} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {SpecNodeRef, SpecNodesRef, specTree, idPath, brokenCount, nodeLabel} =
  await import('../src/dg/designer/node-ref.js');
const {SpecNodeHandler, SpecNodesHandler} = await import('../src/dg/designer/handler.js');
const {SpecDesigner} = await import('../src/dg/designer/view.js');
const {registerPlatformComponents} = await import('../src/dg/viewers/registrations.js');
const {shell} = await import('datagrok-api/grok');

/** Every test runs against a clean document and must leave the live-scope count where it was. */
function designer(name, body) {
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

/** Everything the body assigns to `grok.shell.o` — the platform rebuilds the property panel once
 * per assignment, so the count is what a test asserts. */
async function shellWrites(body) {
  const from = shell.dart.writes.length;
  try {
    await body();
    return shell.dart.writes.slice(from);
  } finally {
    shell.o = null;
  }
}

const INPUT = {
  tag: 'u2-fake-input',
  create: (props) => new TextInput({label: props.label, value: props.value}),
  description: 'Fake input for the designer tests',
  props: [
    {name: 'label', type: 'string', category: 'Appearance'},
    {name: 'value', type: 'string', bindable: true, twoWay: true},
  ],
  events: ['input', 'change'],
  example: {tag: 'u2-fake-input'},
};

const BOX = {
  tag: 'u2-fake-box',
  create: () => new Control(),
  description: 'Fake container',
  props: [],
  acceptsChildren: true,
  example: {tag: 'u2-fake-box'},
};

const LEAF = {
  tag: 'u2-fake-leaf',
  create: () => new Control(),
  description: 'Fake leaf — drops the children a spec gives it',
  props: [],
  example: {tag: 'u2-fake-leaf'},
};

const SPEC = {
  $schema: 'dg-ui/1',
  root: {tag: 'u2-fake-box', name: 'layout', children: [
    {tag: 'u2-fake-box', children: [
      {tag: 'u2-fake-input', name: 'nameInput', props: {label: 'Name', value: 'Aspirin'},
        on: {input: 'cmd:touched'}},
    ]},
    {tag: 'u2-fake-leaf', name: 'leaf', children: [{tag: 'u2-fake-input', name: 'dropped'}]},
    {tag: 'u2-nope', name: 'broken'},
  ]},
};

/** The designer renders against the global registry, so the fakes have to live there too. */
for (const meta of [INPUT, BOX, LEAF])
  globalRegistry.register(meta);

/** Renders SPEC against a private registry, muting the warnings it provokes on purpose. */
function rendered() {
  const reg = new Registry();
  for (const meta of [INPUT, BOX, LEAF])
    reg.register(meta);
  const warn = console.warn;
  console.warn = () => {};
  try {
    return renderSpec(SPEC, new SpecContext(), reg);
  } finally {
    console.warn = warn;
  }
}

function node(name) {
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
  return walk(SPEC.root);
}

designer('path(): the breadcrumb names a node by its spec identity, a nameless one by its tag', () => {
  const instance = rendered();
  assert.equal(new SpecNodeRef(instance, node('nameInput')).path(),
    'layout › u2-fake-box › nameInput');
  assert.equal(new SpecNodeRef(instance, node('layout')).path(), 'layout');
  assert.equal(nodeLabel(node('nameInput')), 'nameInput');
  assert.equal(nodeLabel(SPEC.root.children[0]), 'u2-fake-box');
  instance.dispose();
});

designer('specTree(): the tree is what rendered — a dropped child is in the spec but not in it', () => {
  const instance = rendered();
  const tree = specTree(instance);

  assert.equal(tree.roots.length, 1);
  const root = tree.roots[0];
  assert.equal(root.id, '0');
  assert.equal(root.label, 'layout');
  assert.deepEqual(root.children.map((c) => c.id), ['0.0', '0.1', '0.2']);
  assert.deepEqual(root.children.map((c) => c.label), ['u2-fake-box', 'leaf', 'broken']);
  assert.deepEqual(root.children[0].children.map((c) => c.label), ['nameInput']);
  assert.equal(root.children[1].children, undefined, 'the leaf dropped its child, so the tree has none');
  assert.equal(tree.ids.get(node('nameInput')), '0.0.0');
  instance.dispose();
});

designer('idPath(): every id on the way down, root first', () => {
  assert.deepEqual(idPath('0.0.0'), ['0', '0.0', '0.0.0']);
  assert.deepEqual(idPath('0'), ['0']);
});

designer('a node that could not be built is still selectable, and says why', () => {
  const instance = rendered();
  const ref = new SpecNodeRef(instance, node('broken'));
  assert.equal(brokenCount(instance), 1);
  assert.match(ref.error(), /^u2-nope: /);
  assert.equal(ref.built() instanceof Control, false);
  assert.equal(ref.element(), instance.nodes().get(node('broken')));
  assert.equal(new SpecNodeRef(instance, node('nameInput')).error(), null);
  instance.dispose();
});

designer('the handle reaches the built component and the metadata behind its tag', () => {
  const instance = rendered();
  const ref = new SpecNodeRef(instance, node('nameInput'));
  const built = ref.built();
  assert.ok(built instanceof TextInput);
  assert.equal(ref.element(), built.root);
  assert.equal(ref.meta().tag, 'u2-fake-input');
  assert.equal(ref.meta().events.length, 2);
  instance.dispose();
});

designer('the handler answers for a node ref and captions it by identity and tag', () => {
  const instance = rendered();
  const handler = new SpecNodeHandler();
  assert.equal(handler.type, 'u2-spec-node');
  assert.equal(handler.isApplicable(new SpecNodeRef(instance, node('nameInput'))), true);
  assert.equal(handler.isApplicable({tag: 'u2-fake-input'}), false);
  assert.equal(handler.getCaption(new SpecNodeRef(instance, node('nameInput'))),
    'nameInput (u2-fake-input)');
  assert.equal(handler.getCaption(new SpecNodeRef(instance, SPEC.root.children[0])), 'u2-fake-box');
  instance.dispose();
});

designer('the context panel is the node: identity, properties by category, events and their wiring', () => {
  const instance = rendered();
  const panel = new SpecNodeHandler().renderProperties(new SpecNodeRef(instance, node('nameInput')));
  const caption = panel.children[0];
  assert.equal(caption.classList.contains('u2-designer-caption'), true, 'the panel heads itself');
  assert.equal(caption.textContent, 'nameInput (u2-fake-input)');
  const sections = [...panel.querySelectorAll('h3')].map((h) => h.textContent);
  assert.deepEqual(sections, ['Node', 'Appearance', 'Properties', 'Events'],
    'no bindings section on a node that binds nothing');

  const text = panel.textContent;
  assert.match(text, /u2-fake-input/);
  assert.match(text, /layout › u2-fake-box › nameInput/);
  assert.match(text, /Name/, 'the categorized label prop');
  assert.match(text, /Aspirin/, 'the live value, read through getProperties()');
  assert.match(text, /cmd:touched/, 'the input event, wired');
  instance.dispose();
});

designer('a broken node renders its failure instead of properties it does not have', () => {
  const instance = rendered();
  const ref = new SpecNodeRef(instance, node('broken'));
  const panel = new SpecNodeHandler().renderProperties(ref);
  assert.deepEqual([...panel.querySelectorAll('h3')].map((h) => h.textContent), ['Node']);
  assert.match(panel.textContent, /u2-nope/);
  assert.equal(panel.children[0].textContent, 'broken (u2-nope) — broken',
    'the header says the node did not build');
  instance.dispose();
});

designer('open() selects the new root even when the old selection\'s row id collides', async () => {
  const warn = console.warn;
  console.warn = () => {};
  const view = new SpecDesigner(SPEC);
  const spec = view._instance.spec;
  try {
    view._select(find(spec, 'nameInput'));
    await flush();

    // ids are ordinals: nameInput is '0.0.0' here too, so keyed selection preservation
    // would silently pick otherInput — the collision this test exists for
    const other = {$schema: 'dg-ui/1', root: {tag: 'u2-fake-box', name: 'other', children: [
      {tag: 'u2-fake-box', children: [
        {tag: 'u2-fake-input', name: 'otherInput', props: {label: 'Other'}},
      ]},
    ]}};
    const writes = await shellWrites(async () => {
      view.open(other);
      await flush();
      assert.equal(shell.o.node.name, 'other', 'the new root is what ends up selected');
    });
    assert.match(view.status.value, /^other · /, 'and the status bar agrees');
    assert.equal(writes.some((w) => w?.node?.name === 'otherInput'), false,
      'no transient selection of the colliding node was ever announced');
  } finally {
    console.warn = warn;
    shell.o = null;
    view.dispose();
  }
});

designer('selecting a node always writes the shell object, so a re-click recovers a lost panel', async () => {
  const warn = console.warn;
  console.warn = () => {};
  const view = new SpecDesigner(SPEC);
  console.warn = warn;
  const spec = view._instance.spec;

  assert.ok(shell.o instanceof SpecNodeRef, 'opening selects the root outright');
  assert.equal(shell.o.node, spec.root);

  view._select(find(spec, 'nameInput'));
  await flush();
  assert.equal(shell.o.node, find(spec, 'nameInput'));

  shell.o = null;
  view._select(find(spec, 'nameInput'));
  assert.ok(shell.o instanceof SpecNodeRef, 'the same node again is not a no-op');
  assert.equal(shell.o.node, find(spec, 'nameInput'));
  assert.match(view.status.value, /nameInput/, 'and the status bar still names it');

  shell.o = null;
  fire(view.toolbox.querySelector('.u2-tree'), 'click');
  assert.equal(shell.o?.node, find(spec, 'nameInput'),
    'a click on the row already selected re-asserts it, though no selection signal moved');

  shell.o = null;
  view.dispose();
});

/* WO-9 — the editing gestures: what the palette offers, and what a drop, a context action, a
   keystroke and an undo do to the document, the tree, the selection and the canvas. */

const {Palette} = await import('../src/dg/designer/palette.js');

const ACTIONS = ['Delete', 'Move Up', 'Move Down', 'Duplicate'];

/** A document of its own per test: patches mutate the spec in place. */
function editable() {
  return {$schema: 'dg-ui/1', root: {tag: 'u2-fake-box', name: 'layout', children: [
    {tag: 'u2-fake-box', name: 'group', children: [
      {tag: 'u2-fake-input', name: 'first', props: {label: 'First'}},
      {tag: 'u2-fake-input', name: 'second', props: {label: 'Second'}},
    ]},
    {tag: 'u2-fake-box', name: 'empty'},
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

/** Where the drag layer would land a drop — the rect is the indicator's, not the resolver's input. */
function into(parent, index) {
  return {parent, index, kind: 'into', rect: {x: 0, y: 0, width: 0, height: 0}};
}

function element(view, node) {
  const built = view._instance.nodes().get(node);
  return built instanceof Control ? built.root : built;
}

function stamped(view) {
  return view.root.querySelectorAll('[data-u2-name]').map((el) => el.dataset.u2Name);
}

designer('the palette groups the registry into accordion panes, filters across them, and tags every row', () => {
  const reg = new Registry();
  reg.register({tag: 'u2-p-text', category: 'Inputs', description: 'Single-line text editor',
    create: () => new Control(), props: [], example: {tag: 'u2-p-text'}});
  reg.register({tag: 'u2-p-box', category: 'Containers', description: 'Vertical stack',
    create: () => new Control(), props: [], acceptsChildren: true, example: {tag: 'u2-p-box'}});
  reg.register({tag: 'u2-p-number', category: 'Inputs', description: 'Numeric editor',
    create: () => new Control(), props: [], example: {tag: 'u2-p-number'}});
  const palette = new Palette(reg);
  const panes = () => palette.root.querySelectorAll('.u2-accordion-pane').filter((p) => !p.hidden);
  const headers = () => panes().map((p) => p.querySelector('.u2-accordion-title').textContent);
  const items = () => palette.root.querySelectorAll('.u2-palette-item')
    .filter((el) => !el.hidden && !el.closest('.u2-accordion-pane').hidden)
    .map((el) => el.dataset.u2Tag);

  assert.deepEqual(headers(), ['Containers', 'Inputs'], 'one pane per category, in the fixed pane order');
  assert.deepEqual(items(), ['u2-p-box', 'u2-p-text', 'u2-p-number'], 'registration order within a pane');
  assert.ok(panes().every((p) =>
    p.querySelector('.u2-accordion-header').getAttribute('aria-expanded') === 'true'),
  'every pane starts expanded');
  const text = palette.root.querySelector('[data-u2-tag="u2-p-text"]');
  assert.equal(text.textContent, 'p-text', 'the row is labelled by the tag suffix');
  assert.equal(text.title, 'Single-line text editor');

  palette.filter.value.value = 'numeric';
  assert.deepEqual(items(), ['u2-p-number'], 'the filter reads the description too');
  assert.deepEqual(headers(), ['Inputs'], 'a pane with nothing left in it hides whole');
  palette.filter.value.value = '';
  assert.equal(items().length, 3);
  palette.dispose();
});

designer('a viewer row: its platform label and icon, the Viewers pane sorted and filed after Containers, usage searchable', () => {
  const reg = new Registry();
  const plain = (tag, category, description, extra = {}) => reg.register({tag, category, description,
    create: () => new Control(), props: [], example: {tag}, ...extra});
  plain('u2-p-text', 'Inputs', 'Single-line text editor');
  plain('u2-viewer-scatter-plot', 'Viewers', 'Dots', {label: 'Scatter plot', usage: 'Two numeric columns',
    icon: () => Object.assign(document.createElement('i'), {className: 'fake-glyph'})});
  plain('u2-viewer-bar-chart', 'Viewers', 'Bars', {label: 'Bar chart'});
  plain('u2-p-odd', 'Oddities', 'Unlisted category');
  plain('u2-p-box', 'Containers', 'Vertical stack', {acceptsChildren: true});
  const palette = new Palette(reg);
  const panes = () => palette.root.querySelectorAll('.u2-accordion-pane').filter((p) => !p.hidden);
  const headers = () => panes().map((p) => p.querySelector('.u2-accordion-title').textContent);
  const items = () => palette.root.querySelectorAll('.u2-palette-item')
    .filter((el) => !el.hidden && !el.closest('.u2-accordion-pane').hidden)
    .map((el) => el.dataset.u2Tag);

  assert.deepEqual(headers(), ['Containers', 'Viewers', 'Inputs', 'Oddities'],
    'Viewers right after Containers; a category the order does not name follows, first-seen');
  assert.deepEqual(items(), ['u2-p-box', 'u2-viewer-bar-chart', 'u2-viewer-scatter-plot', 'u2-p-text', 'u2-p-odd'],
    'the Viewers pane sorts by label, the others keep registration order');
  const scatter = palette.root.querySelector('[data-u2-tag="u2-viewer-scatter-plot"]');
  const bar = palette.root.querySelector('[data-u2-tag="u2-viewer-bar-chart"]');
  assert.equal(scatter.textContent, 'Scatter plot', 'the meta label, not the tag suffix');
  assert.notEqual(scatter.querySelector('.u2-palette-item-icon .fake-glyph'), null, 'the platform glyph before it');
  assert.equal(bar.querySelector('.u2-palette-item-icon'), null, 'no glyph, no box');
  assert.equal(scatter.title, 'Two numeric columns', 'usage is the tooltip where there is one');
  assert.equal(bar.title, 'Bars', 'the description otherwise');

  palette.filter.value.value = 'numeric';
  assert.deepEqual(items(), ['u2-viewer-scatter-plot'], 'usage is searchable');
  palette.filter.value.value = 'scatter plot';
  assert.deepEqual(items(), ['u2-viewer-scatter-plot'], 'so is the label');
  palette.dispose();
});

designer('a collapsed palette pane is remembered across rebuilds, and the filter force-opens its matches', () => {
  const reg = new Registry();
  reg.register({tag: 'u2-m-alpha', category: 'AlphaCat', description: 'Alpha thing',
    create: () => new Control(), props: [], example: {tag: 'u2-m-alpha'}});
  reg.register({tag: 'u2-m-beta', category: 'BetaCat', description: 'Beta thing',
    create: () => new Control(), props: [], example: {tag: 'u2-m-beta'}});
  const pane = (p, title) => p.root.querySelectorAll('.u2-accordion-pane')
    .find((el) => el.querySelector('.u2-accordion-title').textContent === title);
  const expanded = (p, title) =>
    pane(p, title).querySelector('.u2-accordion-header').getAttribute('aria-expanded') === 'true';

  const palette = new Palette(reg);
  fire(pane(palette, 'BetaCat').querySelector('.u2-accordion-header'), 'click');
  assert.equal(expanded(palette, 'BetaCat'), false, 'the header click collapsed the pane');

  palette.filter.value.value = 'beta';
  assert.equal(expanded(palette, 'BetaCat'), true, 'a filter match force-opens its pane');
  assert.equal(pane(palette, 'AlphaCat').hidden, true, 'the pane with no match hides');
  palette.filter.value.value = '';
  assert.equal(expanded(palette, 'BetaCat'), false, 'clearing the query restores what the user chose');
  assert.equal(pane(palette, 'AlphaCat').hidden, false);
  palette.dispose();

  const rebuilt = new Palette(reg);
  assert.equal(expanded(rebuilt, 'BetaCat'), false, 'the collapse survives a palette rebuild');
  assert.equal(expanded(rebuilt, 'AlphaCat'), true, 'an untouched pane stays expanded');
  fire(pane(rebuilt, 'BetaCat').querySelector('.u2-accordion-header'), 'click');
  assert.equal(expanded(rebuilt, 'BetaCat'), true, 'and a re-expand is remembered too');
  rebuilt.dispose();
});

designer('a palette drop inserts an auto-named node, selects it, and the status says modified', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const group = find(spec, 'group');
  assert.ok(view.toolbox.querySelector('.u2-palette'), 'the palette sits in the toolbox');
  assert.match(view.status.value, /^layout . 5 nodes$/, view.status.value);

  view._dragLayer.drop({tag: 'u2-fake-input', label: 'fake-input'}, into(group, 2));
  const added = group.children[2];
  assert.equal(added.name, 'fakeInput1', 'the Delphi rule: the camelCased tag suffix plus a free counter');
  assert.equal(shell.o.node, added, 'the dropped node is what ends up selected');
  assert.deepEqual(view._model.roots[0].children[0].children.map((c) => c.label),
    ['first', 'second', 'fakeInput1'], 'the tree followed');
  assert.match(view.status.value, /fakeInput1 . 6 nodes . modified$/, view.status.value);

  view._dragLayer.drop({tag: 'u2-fake-input', label: 'fake-input'}, into(group, 3));
  assert.equal(group.children[3].name, 'fakeInput2', 'the counter skips what the spec already carries');

  shell.o = null;
  view.dispose();
});

designer('a multi-host drops with its seeded pane: one insert, one undo, nothing shared with the next', () => {
  const reg = new Registry();
  reg.register({tag: 'u2-host-box', category: 'Containers', description: 'Private container',
    create: () => new Control(), props: [], acceptsChildren: true, example: {tag: 'u2-host-box'}});
  reg.register({tag: 'u2-host-tabs', category: 'Containers', description: 'Private multi-host',
    create: () => new Control(),
    createWithChildren: (_props, children) => {
      const host = new Control();
      for (const child of children)
        host.root.append(child instanceof Control ? child.root : child);
      return host;
    },
    props: [], childProps: [{name: 'title', type: 'string'}], acceptsChildren: true,
    defaultChildren: [{tag: 'u2-host-box', props: {title: 'Tab 1'}}],
    example: {tag: 'u2-host-tabs'}});
  const view = new SpecDesigner({$schema: 'dg-ui/1', root: {tag: 'u2-host-box', name: 'form1'}},
    undefined, reg);
  const spec = view._instance.spec;

  view._dragLayer.drop({tag: 'u2-host-tabs', label: 'host-tabs'}, into(spec.root, 0));
  const host = spec.root.children[0];
  assert.equal(host.name, 'hostTabs1');
  assert.deepEqual(host.children.map((c) => [c.name, c.props.title]), [['hostBox1', 'Tab 1']],
    'the pane the registration declares, named like anything else the designer inserts');
  assert.equal(shell.o.node, host, 'the container is still what the drop selects');
  assert.ok(stamped(view).includes('hostBox1'), 'and the pane is on the canvas under its own name');
  assert.match(view.status.value, /hostTabs1 . 3 nodes . modified$/, view.status.value);

  view._dragLayer.drop({tag: 'u2-host-tabs', label: 'host-tabs'}, into(spec.root, 1));
  const second = spec.root.children[1];
  assert.deepEqual(second.children.map((c) => c.name), ['hostBox2'], 'the next insert names its own');
  assert.notEqual(second.children[0].props, host.children[0].props, 'and shares nothing with it');

  fire(view.root, 'keydown', {key: 'z', ctrlKey: true});
  fire(view.root, 'keydown', {key: 'z', ctrlKey: true});
  assert.equal(spec.root.children, undefined, 'one undo per insert takes the pane with the host');

  shell.o = null;
  view.dispose();
});

designer('Delete removes the node and selects its parent; the root refuses every action', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const group = find(spec, 'group');
  view._select(find(spec, 'second'));
  view._run('Delete');

  assert.deepEqual(group.children.map((c) => c.name), ['first']);
  assert.equal(shell.o.node, group, 'the selection falls back to the parent the gesture put aside');
  assert.deepEqual(stamped(view), ['layout', 'group', 'first', 'empty'], 'and the canvas re-rendered');

  view._select(spec.root);
  const actions = view._actions();
  assert.deepEqual(actions.map((a) => a.name), ACTIONS);
  assert.deepEqual(actions.filter((a) => a.enabled).map((a) => a.name), [],
    'the root can be neither removed, reordered nor duplicated');

  shell.o = null;
  view.dispose();
});

designer('Move Up and Move Down reorder the document, the canvas and the tree, and stop at the edges', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const group = find(spec, 'group');
  const second = find(spec, 'second');
  view._select(second);
  assert.deepEqual(view._actions().filter((a) => a.enabled).map((a) => a.name),
    ['Delete', 'Move Up', 'Duplicate'], 'the last child cannot move down');

  view._run('Move Up');
  assert.deepEqual(group.children.map((c) => c.name), ['second', 'first']);
  assert.deepEqual(stamped(view), ['layout', 'group', 'second', 'first', 'empty']);
  assert.deepEqual(view._model.roots[0].children[0].children.map((c) => c.label), ['second', 'first']);
  assert.equal(shell.o.node, second, 'the moved node keeps the selection');
  assert.deepEqual(view._actions().filter((a) => a.enabled).map((a) => a.name),
    ['Delete', 'Move Down', 'Duplicate'], 'and now it is the first child');

  view._run('Move Down');
  assert.deepEqual(group.children.map((c) => c.name), ['first', 'second']);

  shell.o = null;
  view.dispose();
});

designer('Duplicate clones the subtree under free names and selects the clone', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const group = find(spec, 'group');
  view._select(group);
  view._run('Duplicate');

  const clone = spec.root.children[1];
  assert.equal(clone.name, 'group1', 'the base is the name stripped of its trailing digits');
  assert.deepEqual(clone.children.map((c) => c.name), ['first1', 'second1'], 'every name in the clone is fresh');
  assert.notEqual(clone.children[0], group.children[0], 'a deep copy, not the same nodes');
  assert.equal(shell.o.node, clone);
  assert.deepEqual(view._model.roots[0].children.map((c) => c.label), ['group', 'group1', 'empty'],
    'inserted right after the original');

  shell.o = null;
  view.dispose();
});

designer('Ctrl+Z on the designer root undoes; a Delete keystroke inside an editor is not a gesture', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const group = find(spec, 'group');
  const second = find(spec, 'second');
  view._select(second);
  view._run('Delete');
  assert.deepEqual(group.children.map((c) => c.name), ['first']);

  fire(view.root, 'keydown', {key: 'z', ctrlKey: true});
  assert.deepEqual(group.children.map((c) => c.name), ['first', 'second']);
  assert.equal(group.children[1], second, 'the very same node object comes back');
  fire(view.root, 'keydown', {key: 'y', ctrlKey: true});
  assert.deepEqual(group.children.map((c) => c.name), ['first'], 'and Ctrl+Y takes it away again');
  fire(view.root, 'keydown', {key: 'z', ctrlKey: true});
  fire(view.root, 'keydown', {key: 'Z', ctrlKey: true, shiftKey: true});
  assert.deepEqual(group.children.map((c) => c.name), ['first'], 'Ctrl+Shift+Z redoes too');
  fire(view.root, 'keydown', {key: 'z', ctrlKey: true});

  view._select(find(spec, 'first'));
  fire(view.root.querySelector('input'), 'keydown', {key: 'Delete'});
  assert.deepEqual(group.children.map((c) => c.name), ['first', 'second'],
    'the Delete key belongs to the text under the cursor');
  fire(view.root, 'keydown', {key: 'Delete'});
  assert.deepEqual(group.children.map((c) => c.name), ['second']);

  shell.o = null;
  view.dispose();
});

designer('a right-click on a tree row selects it and answers the same action list', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const first = find(spec, 'first');
  const actions = view._rowActions({id: '0.0.0', label: 'first', data: first});
  assert.equal(shell.o.node, first, 'the row is selected before its menu is built');
  assert.deepEqual(actions.map((a) => a.name), ACTIONS);
  assert.deepEqual(actions.filter((a) => a.enabled).map((a) => a.name),
    ['Delete', 'Move Down', 'Duplicate']);

  shell.o = null;
  view.dispose();
});

designer('an empty container advertises itself, until something is dropped into it', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const empty = find(spec, 'empty');
  const group = find(spec, 'group');
  assert.equal(element(view, empty).classList.contains('u2-designer-empty'), true);
  assert.equal(element(view, group).classList.contains('u2-designer-empty'), false, 'it holds children');
  assert.equal(element(view, find(spec, 'first')).classList.contains('u2-designer-empty'), false,
    'a leaf is no drop target');

  view._dragLayer.drop({tag: 'u2-fake-input', label: 'fake-input'}, into(empty, 0));
  // the element identity changed with the re-render: the class is stamped on what is rendered now
  assert.equal(element(view, empty).classList.contains('u2-designer-empty'), false);

  view._mode.selected.value = ['run'];
  assert.equal(view.root.classList.contains('u2-designer-running'), true,
    'run mode is what scopes the affordance out (css/designer.css)');

  shell.o = null;
  view.dispose();
});

designer('a structural patch assigns the shell object once, and a re-click still re-asserts it', async () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const group = find(spec, 'group');
  view._select(find(spec, 'second'));
  await flush();

  // the tree rebuild re-selects the survivor programmatically, and the row selection that follows
  // must not announce it a second time — the panel would be rebuilt twice per delete
  const deleted = await shellWrites(async () => {
    view._run('Delete');
    await flush();
  });
  assert.equal(deleted.length, 1, `one assignment per patch, got ${deleted.length}`);
  assert.equal(deleted[0].node, group, 'and it is the survivor the gesture put aside');

  const dropped = await shellWrites(async () => {
    view._dragLayer.drop({tag: 'u2-fake-input', label: 'fake-input'}, into(group, 1));
    await flush();
  });
  assert.equal(dropped.length, 1, `one assignment per patch, got ${dropped.length}`);

  // the load-bearing re-click: the row is already selected, so no signal moves — the click alone
  // has to put the node back on the shell (P2 defect 2b)
  const reclicked = await shellWrites(async () => {
    fire(view.toolbox.querySelector('.u2-tree'), 'click');
    await flush();
  });
  assert.equal(reclicked.length, 1, 'a click on the row already selected re-asserts it');
  assert.equal(reclicked[0].node, view._selected);

  view.dispose();
});

designer('only the left button edits: the other two never arm a drag', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const item = view.toolbox.querySelector('.u2-palette-item');

  fire(view._pane, 'mousedown', {button: 1});
  assert.equal(view._drag, undefined, 'a middle-button press would commit a move on release');
  fire(item, 'mousedown', {button: 2});
  assert.equal(view._drag, undefined, 'and a right-button press a ghost that outlives the menu');

  fire(item, 'mousedown', {button: 0});
  assert.equal(view._drag?.tag, item.dataset.u2Tag);

  view._endDrag(false);
  shell.o = null;
  view.dispose();
});

designer('a designer over a private registry builds and offers that registry alone', () => {
  const reg = new Registry();
  reg.register({tag: 'u2-solo-box', category: 'Containers', description: 'Private container',
    create: () => new Control(), props: [], acceptsChildren: true, example: {tag: 'u2-solo-box'}});
  const view = new SpecDesigner({$schema: 'dg-ui/1', root: {tag: 'u2-solo-box', name: 'only'}},
    undefined, reg);

  assert.equal(view._instance.registry, reg);
  assert.equal(brokenCount(view._instance), 0, 'a tag only this registry knows still builds');
  assert.deepEqual(view.toolbox.querySelectorAll('.u2-palette-item').map((el) => el.dataset.u2Tag),
    ['u2-solo-box'], 'and the palette offers what the canvas can build');

  shell.o = null;
  view.dispose();
});

designer('a refused rename warns and the tree goes back to what the node is called', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const warnings = [];
  const warning = shell.warning;
  shell.warning = (message) => warnings.push(message);
  try {
    const row = {id: '0.0.1', label: 'second', data: find(spec, 'second')};
    view._rename(row, 'renamed');
    assert.equal(find(spec, 'renamed').name, 'renamed');
    assert.deepEqual(view._model.roots[0].children[0].children.map((c) => c.label), ['first', 'renamed']);

    view._rename(row, 'first');
    assert.equal(row.data.name, 'renamed', 'a name another node already carries is refused');
    assert.deepEqual(warnings.map((w) => /already taken/.test(w)), [true]);
    assert.deepEqual(view._model.roots[0].children[0].children.map((c) => c.label), ['first', 'renamed'],
      'and the tree is rebuilt, so the label it wrote itself is gone');
  } finally {
    shell.warning = warning;
    shell.o = null;
    view.dispose();
  }
});

designer('a ribbon click hands focus back to the root, so the next Ctrl+Z is the designer\'s', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  document.body.append(view.root);
  const buttons = view.ribbon();
  // the ribbon verbs are icon buttons now: the tooltip is the accessible name (F15)
  const named = (label) => buttons.find((b) => b.getAttribute('aria-label') === label);
  const group = find(spec, 'group');
  view._select(find(spec, 'second'));

  fire(named('Delete'), 'click');
  assert.deepEqual(group.children.map((c) => c.name), ['first']);
  assert.equal(document.activeElement, view.root, 'the action handed focus back');
  fire(view.root, 'keydown', {key: 'z', ctrlKey: true});
  assert.deepEqual(group.children.map((c) => c.name), ['first', 'second'], 'so Ctrl+Z undoes');

  fire(named('Redo'), 'click');
  assert.deepEqual(group.children.map((c) => c.name), ['first']);
  assert.equal(document.activeElement, view.root, 'Redo hands it back too');
  fire(named('Undo'), 'click');
  assert.deepEqual(group.children.map((c) => c.name), ['first', 'second']);
  assert.equal(document.activeElement, view.root);

  const wrapped = view._rowActions({id: '0.0.1', label: 'second', data: find(spec, 'second')});
  document.activeElement = document.body;
  wrapped.find((a) => a.name === 'Delete').run();
  assert.deepEqual(group.children.map((c) => c.name), ['first']);
  assert.equal(document.activeElement, view.root, 'a context-menu action hands focus back as well');

  shell.o = null;
  view.dispose();
});

/* P3.5 WO-3 — canvas & navigation: Escape to the parent, the design-time class stamps, and tabs /
   accordion headers that stay live under the glass pane (plus reveal-on-select from the tree). */

const rootOf = (built) => built instanceof Control ? built.root : built;

const STRIP = {
  tag: 'u2-fake-strip',
  create: () => new TabStrip(),
  createWithChildren: (_props, children, nodes) => {
    const strip = new TabStrip();
    for (let i = 0; i < children.length; i++)
      strip.addTab({id: `tab-${i}`, label: nodes[i].props?.title ?? `Tab ${i + 1}`,
        content: rootOf(children[i])});
    return strip;
  },
  description: 'The real TabStrip behind a fake tag',
  props: [],
  childProps: [{name: 'title', type: 'string'}],
  acceptsChildren: true,
  example: {tag: 'u2-fake-strip'},
};

const FOLD = {
  tag: 'u2-fake-fold',
  create: () => new Accordion(),
  createWithChildren: (_props, children, nodes) => {
    const fold = new Accordion();
    for (let i = 0; i < children.length; i++)
      fold.addPane(nodes[i].props?.title ?? `Pane ${i + 1}`, rootOf(children[i]));
    return fold;
  },
  description: 'The real Accordion behind a fake tag',
  props: [],
  childProps: [{name: 'title', type: 'string'}],
  acceptsChildren: true,
  example: {tag: 'u2-fake-fold'},
};

for (const meta of [STRIP, FOLD])
  globalRegistry.register(meta);

function paned() {
  return {$schema: 'dg-ui/1', root: {tag: 'u2-fake-box', name: 'layout', children: [
    {tag: 'u2-fake-strip', name: 'tabs', children: [
      {tag: 'u2-fake-box', name: 'pane1', props: {title: 'First'}},
      {tag: 'u2-fake-box', name: 'pane2', props: {title: 'Second'}, children: [
        {tag: 'u2-fake-input', name: 'hiddenInput', props: {label: 'Hidden'}},
      ]},
    ]},
    {tag: 'u2-fake-fold', name: 'fold', children: [
      {tag: 'u2-fake-box', name: 'fp1', props: {title: 'One'}, children: [
        {tag: 'u2-fake-input', name: 'foldedInput', props: {label: 'Folded'}},
      ]},
      {tag: 'u2-fake-box', name: 'fp2', props: {title: 'Two'}},
    ]},
  ]}};
}

designer('Escape selects the parent, stays put on the root, and yields to a drag or an editor', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const second = find(spec, 'second');
  view._select(second);
  const consumed = !fire(view.root, 'keydown', {key: 'Escape'});
  assert.equal(shell.o.node, find(spec, 'group'), 'Esc walks to the parent');
  assert.equal(consumed, true, 'a handled Esc is consumed');
  fire(view.root, 'keydown', {key: 'Escape'});
  assert.equal(shell.o.node, spec.root);
  const atRoot = fire(view.root, 'keydown', {key: 'Escape'});
  assert.equal(shell.o.node, spec.root, 'the root has no parent — a no-op');
  assert.equal(atRoot, true, 'and the keystroke is left alone');

  view._select(second);
  view._drag = {label: '', x: 0, y: 0, active: true, target: null};
  fire(view.root, 'keydown', {key: 'Escape'});
  assert.equal(shell.o.node, second, 'a drag in flight owns Escape');
  view._drag = undefined;

  fire(view.root.querySelector('input'), 'keydown', {key: 'Escape'});
  assert.equal(shell.o.node, second, 'so does the rename box or any other text editor');

  shell.o = null;
  view.dispose();
});

designer('design-time stamps: bounds on the root, the outline class on every node, run mode CSS-scoped', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const rootEl = element(view, spec.root);
  const childEl = element(view, find(spec, 'first'));
  assert.equal(rootEl.classList.contains('u2-designer-root'), true);
  assert.equal(rootEl.classList.contains('u2-designer-node'), true);
  assert.equal(childEl.classList.contains('u2-designer-node'), true);
  assert.equal(childEl.classList.contains('u2-designer-root'), false);

  assert.equal(view.root.classList.contains('u2-designer-outlines'), true, 'outlines default on');
  view._outlines.value = false;
  assert.equal(view.root.classList.contains('u2-designer-outlines'), false);
  view._outlines.value = true;

  // the classes survive run mode — css/designer.css scopes the styling under
  // `:not(.u2-designer-running)`, exactly like the empty-container affordance
  view._mode.selected.value = ['run'];
  assert.equal(view.root.classList.contains('u2-designer-running'), true);
  assert.equal(rootEl.classList.contains('u2-designer-root'), true);
  assert.equal(childEl.classList.contains('u2-designer-node'), true);

  shell.o = null;
  view.dispose();
});

designer('a tab header click at design time activates the pane and selects its node', () => {
  const view = new SpecDesigner(paned());
  const spec = view._instance.spec;
  const strip = view._instance.nodes().get(find(spec, 'tabs'));
  assert.equal(strip.activeTab.peek(), 'tab-0');

  const header = strip.root.querySelectorAll('.u2-tabs-tab')[1];
  document.elementsFromPoint = () => [view._pane, header.querySelector('.u2-tabs-label'), header];
  try {
    fire(view._pane, 'click');
    assert.equal(strip.activeTab.peek(), 'tab-1', 'the header stayed live under the glass pane');
    assert.equal(shell.o.node, find(spec, 'pane2'), 'and the pane node is what got selected');
  } finally {
    delete document.elementsFromPoint;
  }

  shell.o = null;
  view.dispose();
});

designer('an accordion header click toggles its pane, selecting the node only when it expands', () => {
  const view = new SpecDesigner(paned());
  const spec = view._instance.spec;
  const fold = view._instance.nodes().get(find(spec, 'fold'));
  assert.equal(fold.paneAt(1).expanded.peek(), false);

  const header = fold.root.children[1].querySelector('.u2-accordion-header');
  document.elementsFromPoint = () => [view._pane, header];
  try {
    fire(view._pane, 'click');
    assert.equal(fold.paneAt(1).expanded.peek(), true, 'the header stayed live under the glass pane');
    assert.equal(shell.o.node, find(spec, 'fp2'), 'expanding selects the pane node');

    view._select(spec.root);
    fire(view._pane, 'click');
    assert.equal(fold.paneAt(1).expanded.peek(), false, 'the second click collapses');
    assert.equal(shell.o.node, spec.root, 'and collapsing moves no selection');
  } finally {
    delete document.elementsFromPoint;
  }

  shell.o = null;
  view.dispose();
});

designer('selecting a hidden-pane node reveals it: the tab flips, the accordion pane expands', () => {
  const view = new SpecDesigner(paned());
  const spec = view._instance.spec;
  const strip = view._instance.nodes().get(find(spec, 'tabs'));
  const fold = view._instance.nodes().get(find(spec, 'fold'));
  assert.equal(strip.activeTab.peek(), 'tab-0');
  assert.equal(fold.paneAt(0).expanded.peek(), false);

  view._select(find(spec, 'hiddenInput'));
  assert.equal(strip.activeTab.peek(), 'tab-1', 'every tabs ancestor activates the on-path pane');
  view._select(find(spec, 'foldedInput'));
  assert.equal(fold.paneAt(0).expanded.peek(), true, 'every accordion ancestor expands it');
  assert.equal(fold.paneAt(1).expanded.peek(), false, 'the off-path pane is left alone');

  shell.o = null;
  view.dispose();
});

/* P3.5 WO-4 — New form (F3) over the sample gallery (F12): confirm-if-modified through the same
   dirty proxy the status bar reads, and the opened sample is a clone — never the module constant. */

registerAll();

designer('New form asks only when the document is modified; OK discards, Cancel keeps', () => {
  const view = new SpecDesigner(editable());
  view._ribbon.newForm();
  assert.equal(document.querySelector('.u2-dialog'), null, 'a clean document needs no confirmation');
  assert.match(view.status.value, /^form1 . 1 node$/, 'the blank sample is open');
  assert.equal(view._instance.spec.root.tag, 'u2-form');
  assert.notEqual(view._instance.spec.root, SAMPLES[0].spec.root, 'the sample was cloned, not adopted');

  view._dragLayer.drop({tag: 'u2-fake-input', label: 'fake-input'}, into(view._instance.spec.root, 0));
  assert.match(view.status.value, /modified/);
  assert.equal(SAMPLES[0].spec.root.children, undefined, 'the module constant never sees the edit');

  view._ribbon.newForm();
  const dialog = document.querySelector('.u2-dialog');
  assert.ok(dialog, 'a modified document asks first');
  assert.match(dialog.textContent, /Discard changes\?/);
  fire(dialog.querySelectorAll('button').find((b) => b.textContent === 'CANCEL'), 'click');
  assert.match(view.status.value, /modified/, 'Cancel keeps the edits');

  view._ribbon.newForm();
  fire(document.querySelectorAll('.u2-dialog button').find((b) => b.textContent === 'OK'), 'click');
  assert.match(view.status.value, /^form1 . 1 node$/, 'OK discards and opens the blank sample');
  assert.equal(view._editor.peek().canUndo.peek(), false, 'a fresh editor — nothing to undo');

  shell.o = null;
  view.dispose();
});

/* P3.5 WO-5 — multi-selection (F5): Ctrl toggles, Shift ranges, one `SpecNodesRef` per gesture,
   the ancestor-filtered multi-delete as one compound patch, and Esc collapsing to the lead. */

/** Renders the tree rows headless: give the virtual list a height, expand `ids`, refresh. */
function showTree(view, ids) {
  view.toolbox.querySelector('.u2-list').clientHeight = 220;
  view._tree.expanded.value = new Set(ids);
  view._tree.setRoots([...view._model.roots]);
}

function treeRow(view, label) {
  return view.toolbox.querySelectorAll('.u2-list-row')
    .find((row) => row.querySelector('.u2-tree-label')?.textContent === label);
}

designer('Ctrl and Shift clicks build a multi-selection: one SpecNodesRef per gesture', async () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  // in the document, so the glass pane counts as connected and the adorners show
  document.body.append(view.root);
  showTree(view, ['0', '0.0']);
  fire(treeRow(view, 'first').querySelector('.u2-tree-label'), 'click');
  await flush();
  assert.equal(shell.o.node, find(spec, 'first'));

  const toggled = await shellWrites(async () => {
    fire(treeRow(view, 'second').querySelector('.u2-tree-label'), 'click', {ctrlKey: true});
    await flush();
  });
  assert.equal(toggled.length, 1, 'a Ctrl-click is exactly one shell.o assignment');
  assert.ok(toggled[0] instanceof SpecNodesRef);
  assert.deepEqual(toggled[0].refs.map((r) => r.node.name), ['first', 'second']);
  assert.equal(view._selected, find(spec, 'second'), 'the lead follows the toggle');
  assert.match(view.status.value, /^2 selected · 5 nodes/, view.status.value);
  const shown = () => view.root.querySelectorAll('.u2-designer-selected')
    .filter((box) => box.style.display === 'block').length;
  assert.equal(shown(), 2, 'one adorner per member');

  // the anchor is the last plain click — `first` — so the range replaces the set
  const ranged = await shellWrites(async () => {
    fire(treeRow(view, 'empty').querySelector('.u2-tree-label'), 'click', {shiftKey: true});
    await flush();
  });
  assert.equal(ranged.length, 1, 'a Shift-click is exactly one assignment too');
  assert.deepEqual(ranged[0].refs.map((r) => r.node.name), ['first', 'second', 'empty']);
  assert.match(view.status.value, /^3 selected/, view.status.value);
  assert.equal(shown(), 3);

  // the 2b re-click, with a multi-selection up: a plain click on the lead row still re-asserts —
  // and collapses the set, so the tree-root listener's modifier guard is what kept it alive above
  const reclicked = await shellWrites(async () => {
    fire(treeRow(view, 'empty').querySelector('.u2-tree-label'), 'click');
    await flush();
  });
  assert.equal(reclicked.length, 1, 'a plain re-click on the lead re-asserts once');
  assert.ok(reclicked[0] instanceof SpecNodeRef);
  assert.equal(reclicked[0].node, find(spec, 'empty'));
  assert.equal(view._multi.length, 1);
  assert.equal(shown(), 1);

  shell.o = null;
  view.dispose();
});

designer('a canvas Ctrl-click toggles membership, and toggling the last member off is a no-op', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const first = find(spec, 'first');
  const second = find(spec, 'second');
  view._select(first);
  document.elementsFromPoint = () => [view._pane, element(view, second)];
  try {
    fire(view._pane, 'click', {ctrlKey: true});
    assert.deepEqual(view._multi.map((n) => n.name), ['first', 'second']);
    assert.ok(shell.o instanceof SpecNodesRef);
    assert.match(view.status.value, /^2 selected/);

    fire(view._pane, 'click', {ctrlKey: true});
    assert.deepEqual(view._multi.map((n) => n.name), ['first'], 'the second toggle removes it');
    assert.equal(shell.o.node, first, 'back to a single selection');

    document.elementsFromPoint = () => [view._pane, element(view, first)];
    fire(view._pane, 'click', {ctrlKey: true});
    assert.deepEqual(view._multi.map((n) => n.name), ['first'], 'the sole member stays selected');
  } finally {
    delete document.elementsFromPoint;
  }
  shell.o = null;
  view.dispose();
});

designer('multi-delete: the ancestor cover, one compound entry, the survivor is the lead\'s parent',
  async () => {
    const view = new SpecDesigner(editable());
    const spec = view._instance.spec;
    const before = view._instance.dump();
    view._select(find(spec, 'first'));
    view._selection.toggleMulti(find(spec, 'group'));
    view._selection.toggleMulti(find(spec, 'empty'));
    assert.match(view.status.value, /^3 selected/);
    assert.deepEqual(view._verbs.cover().map((n) => n.name), ['group', 'empty'],
      'a member under a selected ancestor rides that ancestor\'s remove');
    assert.deepEqual(view._actions().map((a) => [a.name, a.enabled !== false]),
      [['Delete', true], ['Move Up', false], ['Move Down', false], ['Duplicate', false]],
      'the reorder verbs wait for a single node');

    const writes = await shellWrites(async () => {
      fire(view.root, 'keydown', {key: 'Delete'});
      await flush();
    });
    assert.equal(writes.length, 1, 'a multi-delete assigns the shell object exactly once');
    assert.equal(writes[0].node, spec.root, 'the survivor is the lead\'s former parent');
    assert.equal(spec.root.children, undefined);
    assert.match(view.status.value, /^layout · 1 node · modified/, view.status.value);

    fire(view.root, 'keydown', {key: 'z', ctrlKey: true});
    assert.deepEqual(view._instance.dump(), before, 'ONE undo restores the whole compound');
    assert.equal(view._editor.peek().canUndo.peek(), false);

    shell.o = null;
    view.dispose();
  });

designer('Escape collapses a multi-selection to the lead before walking to the parent', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  view._select(find(spec, 'first'));
  view._selection.toggleMulti(find(spec, 'second'));
  assert.equal(view._multi.length, 2);

  fire(view.root, 'keydown', {key: 'Escape'});
  assert.equal(view._multi.length, 1);
  assert.equal(shell.o.node, find(spec, 'second'), 'Esc collapses to the lead first');
  assert.match(view.status.value, /second · 5 nodes/, view.status.value);

  fire(view.root, 'keydown', {key: 'Escape'});
  assert.equal(shell.o.node, find(spec, 'group'), 'and only then walks to the parent');

  shell.o = null;
  view.dispose();
});

designer('the multi handler: caption, tag counts and member paths — read-only', () => {
  const instance = rendered();
  const handler = new SpecNodesHandler();
  const refs = new SpecNodesRef([
    new SpecNodeRef(instance, node('nameInput')),
    new SpecNodeRef(instance, node('leaf')),
    new SpecNodeRef(instance, SPEC.root.children[0]),
    new SpecNodeRef(instance, SPEC.root),
  ]);
  assert.equal(handler.type, 'u2-spec-nodes');
  assert.equal(handler.isApplicable(refs), true);
  assert.equal(handler.isApplicable(new SpecNodeRef(instance, node('leaf'))), false);
  assert.equal(handler.getCaption(refs), '4 nodes');
  assert.equal(handler.renderMarkup(refs).textContent, '4 nodes');

  const panel = handler.renderProperties(refs);
  assert.deepEqual([...panel.querySelectorAll('h3')].map((h) => h.textContent), ['Tags', 'Nodes']);
  assert.equal(panel.querySelectorAll('input').length, 0, 'nothing to type into');
  const tags = Object.fromEntries([...panel.querySelectorAll('tr')]
    .map((tr) => [tr.children[0].textContent, tr.children[1].textContent]));
  assert.deepEqual(tags, {'u2-fake-input': '1', 'u2-fake-leaf': '1', 'u2-fake-box': '2'});
  assert.match(panel.textContent, /layout › u2-fake-box › nameInput/);
  instance.dispose();
});

designer('a tree rebuild collapses a multi-selection to the lead — everywhere, not just in the tree', async () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  document.body.append(view.root);
  showTree(view, ['0', '0.0']);
  fire(treeRow(view, 'first').querySelector('.u2-tree-label'), 'click');
  await flush();
  fire(treeRow(view, 'second').querySelector('.u2-tree-label'), 'click', {ctrlKey: true});
  await flush();
  assert.equal(view._multi.length, 2);
  assert.match(view.status.value, /^2 selected/);

  // an expand/collapse rebuilds the flat rows (setItems): the list keeps the lead by key and
  // collapses the set — the designer must follow, or Delete would remove invisible members
  showTree(view, ['0', '0.0']);
  await flush();
  assert.equal(view._multi.length, 1, 'the designer collapsed with the tree');
  assert.equal(view._selected, find(spec, 'second'), 'to the lead');
  assert.ok(shell.o instanceof SpecNodeRef, 'the panel shows the single selection');
  assert.equal(shell.o.node, find(spec, 'second'));
  assert.match(view.status.value, /second · 5 nodes/, view.status.value);

  shell.o = null;
  view.dispose();
});

designer('multi-delete with the lead under another selected member: the survivor is the cover\'s parent',
  async () => {
    const view = new SpecDesigner({$schema: 'dg-ui/1', root: {tag: 'u2-fake-box', name: 'layout',
      children: [
        {tag: 'u2-fake-box', name: 'outer', children: [
          {tag: 'u2-fake-box', name: 'inner', children: [
            {tag: 'u2-fake-input', name: 'deep'},
          ]},
        ]},
        {tag: 'u2-fake-box', name: 'empty'},
      ]}});
    const spec = view._instance.spec;
    view._select(find(spec, 'inner'));
    view._selection.toggleMulti(find(spec, 'deep'));
    assert.equal(view._selected, find(spec, 'deep'), 'the lead is a descendant of a selected member');
    assert.deepEqual(view._verbs.cover().map((n) => n.name), ['inner']);

    fire(view.root, 'keydown', {key: 'Delete'});
    await flush();
    assert.equal(shell.o.node, find(spec, 'outer'),
      'the survivor is the parent of the cover member holding the lead — not the root');
    assert.deepEqual(spec.root.children.map((c) => c.name), ['outer', 'empty']);
    assert.equal(find(spec, 'inner'), null);

    shell.o = null;
    view.dispose();
  });

designer('a drop over the canvas but over no node lands in the root — the empty middle is a target', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  document.body.append(view.root);
  view._drag = {tag: 'u2-fake-input', label: 'fake-input', x: 0, y: 0, active: true, target: null};
  document.elementsFromPoint = () => [view._pane];
  try {
    fire(document, 'mousemove', {clientX: 5, clientY: 5});
    assert.equal(view._drag.target?.parent, spec.root, 'a null hit over the pane falls back to the root');
    assert.equal(view._drag.target?.index, 2, 'appended at the end');
    assert.equal(view._drag.target?.kind, 'into');
    view._endDrag(true);
    assert.equal(spec.root.children[2].tag, 'u2-fake-input');

    // over the toolbox the pane is not under the pointer: no target, nothing inserted
    view._drag = {tag: 'u2-fake-input', label: 'fake-input', x: 0, y: 0, active: true, target: null};
    document.elementsFromPoint = () => [view.toolbox];
    fire(document, 'mousemove', {clientX: 5, clientY: 5});
    assert.equal(view._drag.target, null);
    view._endDrag(false);
  } finally {
    delete document.elementsFromPoint;
  }
  shell.o = null;
  view.dispose();
});

designer('a finished palette drop hands focus back, so the next Ctrl+Z is the designer\'s', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  document.body.append(view.root, view.toolbox);
  const group = find(spec, 'group');
  const filter = view.toolbox.querySelector('.u2-palette input');
  filter.focus();
  assert.equal(document.activeElement === filter, true);

  fire(view.toolbox.querySelector('.u2-palette-item'), 'mousedown', {button: 0});
  view._drag.active = true;
  view._drag.target = into(group, 2);
  view._endDrag(true);
  assert.equal(group.children.length, 3, 'the drop landed');
  assert.equal(document.activeElement === view.root, true, 'the drop handed focus off the filter box');

  fire(view.root, 'keydown', {key: 'z', ctrlKey: true});
  assert.equal(group.children.length, 2, 'so Ctrl+Z undoes the drop');

  shell.o = null;
  view.dispose();
});

designer('a click then a dblclick anywhere on a fresh-insert row opens the inline rename', async () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  document.body.append(view.root, view.toolbox);
  showTree(view, ['0', '0.0']);
  view._dragLayer.drop({tag: 'u2-fake-input', label: 'fake-input'}, into(find(spec, 'group'), 2));
  await flush();

  const row = treeRow(view, 'fakeInput1');
  assert.ok(row, 'the fresh row is rendered');
  fire(row.querySelector('.u2-tree-label'), 'click');
  await flush();
  assert.equal(treeRow(view, 'fakeInput1'), row,
    'the selecting click recycled no rows — the dblclick needs a stable target');

  // the row, not the label: a short auto-name leaves most of the row bare, and the gesture
  // aimed at "the row" must still open the editor
  fire(row.querySelector('.u2-tree-row'), 'dblclick');
  const editor = view.toolbox.querySelector('.u2-tree-rename');
  assert.equal(editor?.value, 'fakeInput1', 'the inline rename opened');
  fire(editor, 'keydown', {key: 'Escape'});

  shell.o = null;
  view.dispose();
});

designer('saving to the gallery marks the document clean: the status and the New form guard agree', async () => {
  const store = {};
  globalThis.localStorage = {
    getItem: (k) => store[k] ?? null,
    setItem: (k, v) => store[k] = String(v),
    removeItem: (k) => delete store[k],
  };
  const view = new SpecDesigner(editable());
  try {
    view._dragLayer.drop({tag: 'u2-fake-input', label: 'fake-input'}, into(view._instance.spec.root, 0));
    assert.match(view.status.value, /modified/);

    view._ribbon.saveToGallery();
    const dialog = document.querySelector('.u2-dialog');
    const name = dialog.querySelector('input');
    name.value = 'My form';
    fire(name, 'input');
    fire(dialog.querySelectorAll('button').find((b) => b.textContent === 'OK'), 'click');

    assert.ok(store['u2.designer.gallery']?.includes('My form'), 'the entry was stored');
    assert.equal(view._editor.peek().dirty.peek(), false, 'a saved document is clean');
    assert.equal(view._editor.peek().canUndo.peek(), true, 'without losing its history');
    assert.equal(/modified/.test(view.status.value), false, view.status.value);

    view._ribbon.newForm();
    assert.equal(document.querySelector('.u2-dialog'), null, 'a saved document needs no confirmation');
    assert.match(view.status.value, /^form1 . 1 node$/, 'the blank sample is open');
  } finally {
    delete globalThis.localStorage;
  }
  shell.o = null;
  view.dispose();
});

designer('a tree-row hover highlights the node on the canvas, and a drag suppresses it', () => {
  const view = new SpecDesigner(editable());
  const spec = view._instance.spec;
  const tree = view.toolbox.querySelector('.u2-tree');
  tree.querySelector('.u2-list').clientHeight = 220;
  // a fresh array: the signal skips a same-reference write, and the rows need a re-render
  // now that the list has a height
  view._tree.setRoots([...view._model.roots]);

  const row = tree.querySelectorAll('.u2-list-row')[0];
  fire(row.querySelector('.u2-tree-label'), 'mouseover');
  assert.equal(view._hovered, spec.root, 'the hovered row is the hovered node');
  fire(tree, 'mouseleave');
  assert.equal(view._hovered, null);

  view._drag = {label: '', x: 0, y: 0, active: true, target: null};
  fire(row.querySelector('.u2-tree-label'), 'mouseover');
  assert.equal(view._hovered, null, 'no hover while a drag is in flight');
  view._drag = undefined;

  shell.o = null;
  view.dispose();
});

/* P4 WO-12 — the components tray: a chip is where a non-visual component is selected and acted
   on, and the automation handle a component has no canvas element to carry. */

const TRAY_SPEC = {
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-state', name: 'draft', props: {type: 'string'}}, {tag: 'u2-nope', name: 'bad'}],
  root: {tag: 'u2-form', name: 'form', children: [
    {tag: 'u2-text-input', name: 'nameInput', props: {label: 'Name'}, bind: {value: '$.draft'}},
  ]},
};

const chips = (view) => view.root.querySelectorAll('.u2-designer-chip');

designer('the tray shows a chip per component, stamped, broken-marked and selectable', async () => {
  const warn = console.warn;
  console.warn = () => {};
  const view = new SpecDesigner(JSON.parse(JSON.stringify(TRAY_SPEC)));
  console.warn = warn;
  const [draft] = view._instance.spec.components;

  assert.deepEqual(chips(view).map((c) => c.dataset.u2Name), ['draft', 'bad']);
  assert.equal(chips(view)[1].classList.contains('u2-designer-chip-broken'), true);
  assert.equal(chips(view)[0].classList.contains('u2-designer-chip-broken'), false);

  fire(chips(view)[0], 'click');
  await flush();
  assert.equal(shell.o.node, draft, 'a chip click is the ordinary selection pipeline');
  assert.equal(chips(view)[0].classList.contains('u2-designer-chip-selected'), true);
  assert.match(view.status.value, /^draft · 4 nodes/, 'components count as nodes');

  // the reorder verbs have no parent to work against; Delete knows which removal it is
  const actions = view._actions();
  assert.deepEqual(actions.slice(0, 4).map((a) => [a.name, a.enabled !== false]),
    [['Delete', true], ['Move Up', false], ['Move Down', false], ['Duplicate', false]]);
  actions[0].run();
  assert.equal(view._instance.spec.components.includes(draft), false);
  assert.deepEqual(chips(view).map((c) => c.dataset.u2Name), ['bad']);
  assert.equal(view._instance.root.querySelectorAll('.u2-spec-error').length, 1,
    'and the input that bound to it renders broken, contained');

  view._editor.peek().undo();
  assert.equal(view._instance.spec.components[0] === draft, true);
  assert.deepEqual(chips(view).map((c) => c.dataset.u2Name), ['draft', 'bad']);
  assert.equal(view._instance.root.querySelectorAll('.u2-spec-error').length, 0);

  shell.o = null;
  view.dispose();
});

designer('a non-visual palette drag lands on the tray, wherever over the designer it is dropped', () => {
  const warn = console.warn;
  console.warn = () => {};
  const view = new SpecDesigner({$schema: 'dg-ui/1', root: {tag: 'u2-form', name: 'form'}});
  console.warn = warn;

  view._drag = {tag: 'u2-state', label: 'state', x: 0, y: 0, active: true, target: null,
    tray: true, over: true};
  view._endDrag(true);
  const [added] = view._instance.spec.components;
  assert.equal(added.name, 'state1', 'named off the tag like every palette insert');
  assert.equal(added.props.type, 'string', 'and seeded from the registry defaults');
  assert.deepEqual(chips(view).map((c) => c.textContent), ['state1']);
  assert.equal(shell.o.node, added, 'the inserted component is selected');
  assert.match(view.status.value, /^state1 · 2 nodes · modified/);

  view._editor.peek().undo();
  assert.equal(view._instance.spec.components, undefined);
  assert.equal(chips(view).length, 0);

  shell.o = null;
  view.dispose();
});

designer('the Design/Run toggle rebuilds the tray for the mode it switches to', () => {
  const warn = console.warn;
  console.warn = () => {};
  const view = new SpecDesigner({$schema: 'dg-ui/1',
    components: [{tag: 'u2-state', name: 'draft'}], root: {tag: 'u2-form', name: 'form'}});
  console.warn = warn;

  assert.equal(view._instance.designTime, true, 'the designer opens in design mode');
  const before = view._instance.node('draft');
  view._mode.selected.value = ['run'];
  assert.equal(view._instance.designTime, false);
  assert.equal(view._instance.node('draft') === before, false, 'every component is built again');

  shell.o = null;
  view.dispose();
});

/* P4.5 WO-2 — a file or a query dragged out of Browse lands on the tray. The platform half is the
   real one: `ui.makeDroppable` registers dart-side and the stub records the registration, so these
   drive the very callbacks the platform's own mouseup would, with a fake func runner behind them. */

const DG = await import('datagrok-api/dg');
const ui = await import('datagrok-api/ui');
const {backends} = await import('../src/sources/backends.js');
const {sourceStatus} = await import('../src/dg/designer/source-status.js');

/** `active()` reads `_pane.isConnected`, so a designer nobody attached is a designer in run mode. */
function mounted(view) {
  document.body.append(view.root);
  return view;
}

const zone = () => ui.drops[ui.drops.length - 1];
const blank = () => new SpecDesigner({$schema: 'dg-ui/1', root: {tag: 'u2-form', name: 'form'}});
const csv = (path) => new DG.FileInfo(path);
const FUNCS = [
  {name: 'OpenFile', kind: 'function', inputs: [new Property('fullPath', 'string')],
    outputs: [new Property('result', 'int')]},
  // hand-written SQL: the class that must never run without being asked
  {name: 'orders', kind: 'query', inputs: [new Property('days', 'int')],
    outputs: [new Property('rows', 'int')]},
];

/** The designer over a fake platform: what ran, and what the shell was told about it. */
function withBackend(body) {
  const saved = backends.funcRunner;
  const ran = [];
  const said = [];
  const {info, warning} = shell;
  backends.funcRunner = {
    find: (name) => FUNCS.find((f) => f.name === name) ?? null,
    run: (name) => {
      ran.push(name);
      return Promise.resolve({result: 4, rows: 4});
    },
  };
  shell.info = (text) => said.push(text);
  shell.warning = (text) => said.push(text);
  return Promise.resolve(body(ran, said)).finally(() => {
    backends.funcRunner = saved;
    shell.info = info;
    shell.warning = warning;
    shell.o = null;
  });
}

/** Past `AsyncValue.refresh()`, which bypasses the debounce and answers when the run lands. */
const landed = () => new Promise((resolve) => setTimeout(resolve, 100));

designer('a dropped file lands on the tray as an OpenFile source — selected, and run once', () =>
  withBackend(async (ran, said) => {
    const view = mounted(blank());
    zone().doDrop({dragObject: csv('System:Demo/demog.csv'), handled: false});
    await landed();

    assert.deepEqual(view.dump().components, [{tag: 'u2-func-source', name: 'demog1',
      props: {func: 'OpenFile', params: {fullPath: 'System:Demo/demog.csv'}}}]);
    assert.equal(shell.o.node, view._instance.spec.components[0], 'the new source is selected');
    assert.deepEqual(ran, ['OpenFile'], 'the drop IS the explicit act DD9 permits (OR-D1)');
    assert.deepEqual(said, ['demog1: ready'], 'and it says what came of it');

    view.dispose();
  }));

designer('two files are one gesture: two sources, one undo entry', () => withBackend(async (ran) => {
  const view = mounted(blank());
  zone().doDrop({dragObject: [csv('a/one.csv'), csv('a/two.csv')], handled: false});
  await landed();

  assert.deepEqual(view.dump().components.map((c) => c.name), ['one1', 'two1']);
  assert.equal(ran.length, 2, 'each file is read once');
  view._editor.peek().undo();
  assert.equal(view.dump().components, undefined, 'one Ctrl+Z takes the whole drop back');

  view.dispose();
}));

designer('a file nothing can open is refused by name, and nothing lands', () =>
  withBackend(async (ran, said) => {
    const view = mounted(blank());
    zone().doDrop({dragObject: [csv('a/report.pdf'), new DG.FileInfo('a/tables', {directory: true})],
      handled: false});
    await landed();

    assert.equal(view.dump().components, undefined);
    assert.equal(said.length, 1, 'one message per drop, however many items it refused');
    assert.match(said[0], /report\.pdf.*tables/, `"${said[0]}"`);
    assert.deepEqual(ran, []);

    view.dispose();
  }));

designer('a drop in Run mode, and a drop into a disposed designer, do nothing at all', () =>
  withBackend(async () => {
    const view = mounted(blank());
    const drop = zone();
    view._mode.selected.value = ['run'];
    drop.doDrop({dragObject: csv('a/one.csv'), handled: false});
    await landed();
    assert.equal(view.dump().components, undefined, 'there is no tray on screen to point at');

    view._mode.selected.value = ['design'];
    view.dispose();
    // the dart-side registration outlives the designer — `active()` is the whole defence
    drop.doDrop({dragObject: csv('a/one.csv'), handled: false});
    assert.equal(drop.acceptDrag({dragObject: csv('a/one.csv')}), false);
  }));

designer('the tray says where the drop will land, for the whole drag, and hover stays out of it', () =>
  withBackend(async () => {
    const view = mounted(blank());
    const drop = zone();
    const target = () => view._tray.root.classList.contains('u2-designer-tray-drop');

    view._hovered = view._instance.spec.root;
    drop.onBeginDrag();
    assert.equal(target(), true);
    assert.equal(view._hovered, null, 'the adorner of what the pointer left behind goes with it');
    assert.equal(view._dragging, true, 'so the canvas draws no hover adorner under a platform drag');
    fire(view._pane, 'mousemove', {clientX: 5, clientY: 5});
    assert.equal(view._hovered, null);
    drop.onEndDrag();
    assert.equal(target(), false);

    view.dispose();
  }));

designer('a dropped query inserts with no params and says out loud that it has not run', () =>
  withBackend(async (ran) => {
    const view = mounted(blank());
    const query = new DG.DataQuery('orders');
    DG.Func.registry = [query];
    zone().doDrop({dragObject: query, handled: false});
    await landed();
    DG.Func.registry = [];

    assert.deepEqual(view.dump().components,
      [{tag: 'u2-func-source', name: 'orders1', props: {func: 'orders'}}],
      'no params, no designData — the source computes its policy from the kind (OR-D2)');
    assert.deepEqual(ran, [], 'hand-written SQL is never run by a gesture that only dropped it');
    // the chip and the panel must read as "not run yet" — never as an error, and never as silence
    const chip = view.root.querySelector('.u2-designer-chip');
    assert.match(chip.title, /not run yet — use Refresh/);
    assert.equal(chip.classList.contains('u2-designer-chip-error'), false);
    assert.equal(sourceStatus(view._instance.node('orders1')).state, 'idle');

    view.dispose();
  }));

/* P4.5 — the document is the designer's own copy. The editor patches it in place (DD10), so a spec
   a host hands over has to stay the constant it wrote: one host's `DESIGNER_SPEC` was arriving at
   the next designer of the session already edited. */

designer('the spec a host passes in is never the document: edits, undo and a second designer leave it whole', () => {
  const host = editable();
  const before = JSON.stringify(host);
  const view = new SpecDesigner(host);
  const doc = view._instance.spec;
  assert.notEqual(doc, host, 'the canvas holds a copy from the first frame');

  view._dragLayer.drop({tag: 'u2-fake-input', label: 'fake-input'}, into(find(doc, 'group'), 2));
  view._select(find(doc, 'first'));
  view._run('Delete');
  view._rename(find(doc, 'second'), 'renamed');
  view._editor.peek().undo();

  assert.notDeepEqual(view.dump(), JSON.parse(before), 'the document took the edits');
  assert.equal(JSON.stringify(host), before, 'and the host object is byte-identical to what it was');

  const second = new SpecDesigner(host);
  assert.deepEqual(second.dump(), JSON.parse(before),
    'so the next designer of the session opens the spec, not the first one\'s work');

  shell.o = null;
  second.dispose();
  view.dispose();
});

designer('a JSON string opens the same document, and reopening it is a fresh one', () => {
  const text = JSON.stringify(editable());
  const view = new SpecDesigner(text);
  assert.deepEqual(view.dump(), JSON.parse(text));

  view._dragLayer.drop({tag: 'u2-fake-input', label: 'fake-input'}, into(find(view._instance.spec, 'group'), 2));
  assert.match(view.status.value, /6 nodes . modified$/, view.status.value);

  view.open(text);
  assert.deepEqual(view.dump(), JSON.parse(text), 'the string parses fresh every time');
  assert.equal(view._editor.peek().canUndo.peek(), false);

  shell.o = null;
  view.dispose();
});

/* UX-B — the panel after the canvas moved under it: the Design/Run toggle rebuilds every source and
   viewer, and a re-render can turn a placeholder into a component or back; the platform renders
   the panel from the ref it was handed, so each of those is one fresh `shell.o`. */

designer('the Design/Run toggle re-issues the panel for the selected node — once per toggle, from what the new mode built',
  async () => {
    const warn = console.warn;
    console.warn = () => {};
    const view = new SpecDesigner({$schema: 'dg-ui/1',
      components: [{tag: 'u2-state', name: 'draft'}], root: {tag: 'u2-form', name: 'form'}});
    console.warn = warn;
    const draft = view._instance.spec.components[0];
    view._select(draft);
    await flush();

    const toRun = await shellWrites(async () => {
      view._mode.selected.value = ['run'];
      await flush();
    });
    assert.equal(toRun.length, 1, `one assignment per toggle, got ${toRun.length}`);
    assert.equal(toRun[0].node, draft);
    assert.equal(toRun[0].built(), view._instance.node('draft'),
      'the ref reads what Run mode built — the Status section the old panel kept reading is gone with it');

    // the platform's own gear, clicked on a viewer in Run mode, puts that viewer on the panel
    const toDesign = await shellWrites(async () => {
      shell.o = new DG.Viewer({type: 'Grid'});
      view._mode.selected.value = ['design'];
      await flush();
    });
    assert.equal(toDesign.length, 2, 'the gear\'s own write, then the one the toggle makes');
    assert.ok(toDesign[1] instanceof SpecNodeRef, 'leaving Run mode puts the selection back on the panel');
    assert.equal(toDesign[1].node, draft);

    const same = await shellWrites(async () => {
      view._mode.selected.value = ['design'];
      await flush();
    });
    assert.equal(same.length, 0, 'the mode it is already in is no toggle');
    view.dispose();
  });

/** The rows of the panel's Status block — the table under its h3. */
function statusRows(panel) {
  const kids = [...panel.children];
  const i = kids.findIndex((el) => el.tagName === 'H3' && el.textContent === 'Status');
  const rows = {};
  for (const tr of kids[i + 1]?.querySelectorAll('tr') ?? [])
    rows[tr.children[0].textContent] = tr.children[1].textContent;
  return rows;
}

const errorRow = (panel) => [...panel.querySelectorAll('tr')].find((tr) => tr.children[0].textContent === 'Error');

/* The platform renders the panel one `shell.o` write behind, so the one re-issue above would leave
   the old state on screen until the next gesture: the selection first redraws, in place, the parts
   of the panel it rendered that read the node's state — the caption, the Node table, the Status
   block — and then asks once. */

designer('the Design/Run toggle redraws the Status block of the panel on screen from what Run mode built — still one write',
  async () => {
    const saved = backends.funcRunner;
    const frames = {name: 'Frames', kind: 'function', inputs: [], outputs: [new Property('orders', 'dataframe')]};
    backends.funcRunner = {find: (name) => name === 'Frames' ? frames : null,
      run: () => Promise.resolve({orders: new DataFrame([{name: 'customer', type: 'string'}],
        [{customer: 'Bayer'}, {customer: 'Roche'}])})};
    const warn = console.warn;
    console.warn = () => {};
    const view = new SpecDesigner({$schema: 'dg-ui/1',
      components: [{tag: 'u2-func-source', name: 'orders', props: {func: 'Frames', debounce: 0}}],
      root: {tag: 'u2-form', name: 'form'}});
    console.warn = warn;
    try {
      view._select(view._instance.spec.components[0]);
      const shown = new SpecNodeHandler().renderProperties(shell.o);
      document.body.append(shown);
      assert.deepEqual(statusRows(shown), {State: 'idle — not run yet; use Refresh', Rows: '0'});

      const writes = await shellWrites(async () => {
        view._mode.selected.value = ['run'];
        for (let i = 0; i < 5; i++)
          await flush();
      });
      assert.equal(writes.length, 1, 'one assignment per toggle');
      assert.deepEqual(statusRows(shown), {State: 'ready', Rows: '2'},
        'the same panel element reads the source Run mode built and ran');
    } finally {
      backends.funcRunner = saved;
      view.dispose();
    }
  });

designer('a re-render that builds a broken node, or breaks a built one, re-issues the panel; one that does neither is silent',
  async () => {
    const warn = console.warn;
    console.warn = () => {};
    const view = new SpecDesigner({$schema: 'dg-ui/1', root: {tag: 'u2-fake-box', name: 'layout', children: [
      {tag: 'u2-fake-input', name: 'ghost', props: {label: 'Ghost'}, bind: {value: '$.nowhere'}},
    ]}}, new SpecContext({data: {reagent: signal('Aspirin')}}));
    console.warn = warn;
    const ghost = find(view._instance.spec, 'ghost');
    view._select(ghost);
    await flush();
    assert.notEqual(shell.o.error(), null, 'selected while broken');
    const shown = new SpecNodeHandler().renderProperties(shell.o);
    document.body.append(shown);
    assert.equal(shown.querySelector('.u2-designer-caption').textContent, 'ghost (u2-fake-input) — broken');
    assert.notEqual(errorRow(shown), undefined);

    const built = await shellWrites(async () => {
      view._apply({op: 'set-bind', node: ghost, name: 'value', path: '$.reagent'});
      await flush();
    });
    assert.equal(built.length, 1, `the flip is one assignment, got ${built.length}`);
    assert.equal(built[0].node, ghost);
    assert.equal(built[0].error(), null, 'and the new ref says the node is fine');
    assert.equal(shown.querySelector('.u2-designer-caption').textContent, 'ghost (u2-fake-input)',
      'the panel on screen lost its "— broken" without waiting for the platform');
    assert.equal(errorRow(shown), undefined, 'and its Error row');

    const quiet = await shellWrites(async () => {
      view._apply({op: 'set-prop', node: ghost, name: 'label', value: 'Spirit'});
      await flush();
    });
    assert.equal(quiet.length, 0, 'a re-render that keeps the node built refreshes the panel in place (prop-editors)');

    console.warn = () => {};
    const broken = await shellWrites(async () => {
      view._apply({op: 'set-bind', node: ghost, name: 'value', path: '$.nowhere'});
      await flush();
    });
    console.warn = warn;
    assert.equal(broken.length, 1);
    assert.notEqual(broken[0].error(), null, 'the panel is told the node broke');
    assert.equal(shown.querySelector('.u2-designer-caption').textContent, 'ghost (u2-fake-input) — broken');
    assert.notEqual(errorRow(shown), undefined, 'the Error row is back in place');
    view.dispose();
  });

designer('an "Add binding…" pick that breaks the viewer is one shell.o write — the selection\'s re-issue, not the panel\'s',
  async () => {
    WidgetDescriptor.registry = [new WidgetDescriptor('Scatter plot',
      [new Property('xColumnName', 'string', {category: 'Data'})])];
    WidgetDescriptor.prototype.createIcon = () => null;
    const reg = new Registry();
    registerPlatformComponents(reg);
    const scope = new Scope();
    const df = new DataFrame([{name: 'city', type: 'string'}], [{city: 'Kyiv'}], 'orders');
    // the platform refuses a look it cannot build: the bound value is what breaks the viewer
    const fromType = DG.Viewer.fromType;
    DG.Viewer.fromType = (type, table, options) => {
      if (options?.xColumnName === 'boom')
        throw new Error('no such column: boom');
      return fromType(type, table, options);
    };
    const warn = console.warn;
    console.warn = () => {};
    const view = new SpecDesigner({$schema: 'dg-ui/1', root: {tag: 'u2-div-v', name: 'box', children: [
      {tag: 'u2-viewer-scatter-plot', name: 'plot', bind: {table: '$.orders'}}]}},
    new SpecContext({data: {orders: dfBindings(signal(df), scope), boom: signal('boom')}}), reg);
    try {
      const plot = find(view._instance.spec, 'plot');
      view._select(plot);
      await flush();
      const shown = new SpecNodeHandler().renderProperties(shell.o);
      document.body.append(shown);
      const select = shown.querySelector('[data-u2-prop="add-binding"] select');
      select.value = 'xColumnName';
      fire(select, 'change');
      fire(shown.querySelector('[data-u2-bind-pick="add-binding"]'), 'click');
      const dialog = document.querySelector('.u2-bind-picker');
      const list = dialog.querySelector('.u2-list');
      list.clientHeight = 400;
      fire(list, 'scroll');
      fire([...dialog.querySelectorAll('.u2-list-row')]
        .find((row) => (row.querySelector('.u2-tree-label')?.textContent ?? '').startsWith('boom')), 'click');
      const writes = await shellWrites(async () => {
        [...dialog.querySelectorAll('.u2-dialog-buttons button')].find((b) => b.textContent === 'OK').click();
        await flush();
      });
      assert.equal(plot.bind.xColumnName, '$.boom');
      assert.equal(writes.length, 1, `the flip re-issued the panel and the pick did not ask again, got ${writes.length}`);
      assert.notEqual(writes[0].error(), null, 'the one ref says the node broke');
      assert.equal(shown.querySelector('.u2-designer-caption').textContent, 'plot (u2-viewer-scatter-plot) — broken',
        'and the panel on screen says so already');
    } finally {
      console.warn = warn;
      DG.Viewer.fromType = fromType;
      view.dispose();
      scope.dispose();
      WidgetDescriptor.registry = [];
      delete WidgetDescriptor.prototype.createIcon;
      platform.reset();
    }
  });

designer('the selection adorner follows the lead element\'s size: one ResizeObserver, moved with the selection, gone with the designer',
  async () => {
    const Real = globalThis.ResizeObserver;
    const observers = [];
    globalThis.ResizeObserver = class {
      constructor(callback) {
        this.callback = callback;
        this.targets = [];
        this.disconnected = false;
        observers.push(this);
      }
      observe(el) { this.targets.push(el); }
      unobserve(el) { this.targets = this.targets.filter((t) => t !== el); }
      disconnect() {
        this.disconnected = true;
        this.targets = [];
      }
    };
    try {
      const view = mounted(new SpecDesigner(editable()));
      const spec = view._instance.spec;
      // the pane reads as connected only once the designer is: the first reposition after mounting
      // (the first hover, in the app) is what attaches the watch
      view._selection.reposition();
      // the tree's virtual list measures its scroller through an observer of its own
      const watching = observers.filter((o) => o.targets.includes(element(view, spec.root)));
      assert.equal(watching.length, 1, 'the root is selected on open — and watched, once');
      const [observer] = watching;

      view._select(find(spec, 'first'));
      await flush();
      assert.deepEqual(observer.targets, [element(view, find(spec, 'first'))], 'the watch moves with the lead');

      const box = view._selection.box;
      box.style.display = 'none';
      observer.callback();
      assert.equal(box.style.display, 'block', 'a size change re-outlines the lead');

      view._mode.selected.value = ['run'];
      assert.deepEqual(observer.targets, [], 'run mode draws no adorner, so nothing is watched');
      view._mode.selected.value = ['design'];
      assert.equal(observer.targets.length, 1);

      shell.o = null;
      view.dispose();
      assert.equal(observer.disconnected, true, 'disposed with the designer');
    } finally {
      globalThis.ResizeObserver = Real;
    }
  });
