/* The read-only inspector (P2): the models a designer derives from a rendered spec — breadcrumb
   path, structure tree, broken-node count — and what the node ObjectHandler renders from them.
   `DG.ObjectHandler` comes from tests/dg-stub.mjs. The glass pane hit-test and the shell view need
   a real browser and a running platform; they are the live pass, not this file. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {TextInput} from '../src/components/text-input.js';
import {Registry, registry as globalRegistry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';

register('./dg-stub.mjs', import.meta.url);
const {SpecNodeRef, specTree, idPath, brokenCount, nodeLabel} =
  await import('../src/dg/designer/node-ref.js');
const {SpecNodeHandler} = await import('../src/dg/designer/handler.js');
const {SpecDesigner} = await import('../src/dg/designer/view.js');
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
  try {
    view._select(node('nameInput'));
    await flush();

    // ids are ordinals: nameInput is '0.0.0' here too, so keyed selection preservation
    // would silently pick otherInput — the collision this test exists for
    const writes = [];
    const orig = Object.getOwnPropertyDescriptor(Object.getPrototypeOf(shell), 'o') ??
      {get: () => shell._o, set: (v) => shell._o = v};
    let current = shell.o;
    Object.defineProperty(shell, 'o', {
      configurable: true,
      get: () => current,
      set: (v) => {
        current = v;
        writes.push(v);
      },
    });
    const other = {$schema: 'dg-ui/1', root: {tag: 'u2-fake-box', name: 'other', children: [
      {tag: 'u2-fake-box', children: [
        {tag: 'u2-fake-input', name: 'otherInput', props: {label: 'Other'}},
      ]},
    ]}};
    view.open(other);
    await flush();

    assert.equal(shell.o.node.name, 'other', 'the new root is what ends up selected');
    assert.match(view.status.value, /^other · /, 'and the status bar agrees');
    assert.equal(writes.some((w) => w?.node?.name === 'otherInput'), false,
      'no transient selection of the colliding node was ever announced');
    delete shell.o;
    shell.o = current;
    void orig;
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

  assert.ok(shell.o instanceof SpecNodeRef, 'opening selects the root outright');
  assert.equal(shell.o.node, SPEC.root);

  view._select(node('nameInput'));
  await flush();
  assert.equal(shell.o.node, node('nameInput'));

  shell.o = null;
  view._select(node('nameInput'));
  assert.ok(shell.o instanceof SpecNodeRef, 'the same node again is not a no-op');
  assert.equal(shell.o.node, node('nameInput'));
  assert.match(view.status.value, /nameInput/, 'and the status bar still names it');

  shell.o = null;
  fire(view.toolbox.querySelector('.u2-tree'), 'click');
  assert.equal(shell.o?.node, node('nameInput'),
    'a click on the row already selected re-asserts it, though no selection signal moved');

  shell.o = null;
  view.dispose();
});
