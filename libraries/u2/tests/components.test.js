/* The components tray (WO-12): the envelope's second tree, built before the form and started
   after it, edited through its own two patches, and rebuilt whole when the designer flips the
   Design/Run toggle. Everything bound to a component comes back with it — that is the crux, so
   most of what follows watches component identity on both sides of a change. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Component} from '../src/core/component.js';
import {signal} from '../src/core/signals.js';
import {TextInput} from '../src/components/text-input.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {SpecEditor} from '../src/spec/editor.js';
import {registerAll} from '../src/spec/registrations.js';
import {StateSource} from '../src/sources/state.js';
import {brokenCount, specTree} from '../src/dg/designer/node-ref.js';

/** Every test runs against a clean document and must leave the live-scope count where it was. */
function tray(name, body) {
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

/** Every source ever built, in construction order — how a test sees what was re-mounted. */
const built = [];

/** A source that records what it was built with and answers one step: enough to prove the
   construction order, the design-time seam and the dependent re-render. */
class FakeSource extends Component {
  constructor(props, env) {
    super();
    this.value = signal(props.text ?? '');
    this.designTime = env.designTime;
    this.env = env;
    this.started = null;
    this.starts = 0;
    built.push(this);
  }

  start() {
    // the second phase: by now the whole form exists, so a param bound to an input resolves
    this.starts++;
    const path = this.env.subBinds['params.from'];
    this.bound = path === undefined ? null : this.env.resolve(path);
    this.started = this.bound?.peek() ?? null;
  }

  bindStep(name) {
    return name === 'value' || name === '' ? this.value : null;
  }

  bindProps() {
    return [{name: 'value', type: 'string', writable: true}];
  }
}

const SOURCE = {
  tag: 'u2-fake-source',
  category: 'Data',
  visual: false,
  createComponent: (props, env) => new FakeSource(props, env),
  description: 'Fake non-visual source',
  props: [{name: 'text', type: 'string'}, {name: 'params', type: 'object', subBindable: true}],
  example: {tag: 'u2-fake-source', name: 'src'},
};

function registry() {
  const reg = new Registry();
  registerAll(reg);
  reg.register(SOURCE);
  return reg;
}

function open(spec, ctx = new SpecContext()) {
  const reg = registry();
  built.length = 0;
  const instance = renderSpec(JSON.parse(JSON.stringify(spec)), ctx, reg, {designTime: false});
  return {instance, editor: new SpecEditor(instance), reg};
}

const editors = (instance) => instance.root.querySelectorAll('input');

function type(editor, text) {
  editor.value = text;
  fire(editor, 'input');
}

function captureWarnings(body) {
  const warnings = [];
  const original = console.warn;
  console.warn = (message) => warnings.push(message);
  try {
    body();
  } finally {
    console.warn = original;
  }
  return warnings;
}

const STATE_SPEC = {
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-state', name: 'draft', props: {type: 'string', initial: 'seed'}}],
  root: {tag: 'u2-form', name: 'form', children: [
    {tag: 'u2-text-input', name: 'nameInput', props: {label: 'Name'}, bind: {value: '$.draft'}},
  ]},
};

tray('envelope: components round-trip, before the root, exactly as written', () => {
  const {instance} = open(STATE_SPEC);
  const dump = instance.dump();
  assert.deepEqual(Object.keys(dump), ['$schema', 'components', 'root'],
    'authoring order: what the form binds to is declared first');
  assert.deepEqual(dump.components, STATE_SPEC.components);
  assert.equal(instance.nodes().size, 3, 'the tray node is a node like any other');
  instance.dispose();
});

tray('envelope: parseSpec refuses a tray that is not a list of nodes', () => {
  const ctx = new SpecContext();
  assert.throws(() => renderSpec({$schema: 'dg-ui/1', root: {tag: 'div'}, components: {}}, ctx),
    /"components" must be an array/);
  assert.throws(() => renderSpec({$schema: 'dg-ui/1', root: {tag: 'div'}, components: [{}]}, ctx),
    /must be a node with a "tag"/);
});

tray('u2-state: the default step binds a text input two ways, end to end', () => {
  const {instance} = open(STATE_SPEC);
  const state = instance.node('draft');
  assert.equal(state instanceof StateSource, true);
  assert.equal(state.value.peek(), 'seed');

  const editor = editors(instance)[0];
  assert.equal(editor.value, 'seed', '$.draft resolves through bindStep(\'\')');
  state.value.value = 'from the source';
  assert.equal(editor.value, 'from the source');
  type(editor, 'typed');
  assert.equal(state.value.peek(), 'typed', 'and the input writes back');

  state.getFunctions().find((f) => f.name === 'clear').apply({});
  assert.equal(state.value.peek(), 'seed', 'clear returns to the initial value');
  assert.equal(editor.value, 'seed');
  instance.dispose();
});

tray('u2-state: an initial value that does not fit the type warns and starts empty', () => {
  const warnings = captureWarnings(() => {
    const source = new StateSource({type: 'int', initial: 'nope'});
    assert.equal(source.value.peek(), 0);
    source.dispose();
  });
  assert.equal(warnings.length, 1, warnings.join('; '));
  assert.match(warnings[0], /expects int/);
});

tray('a component that cannot be built maps the placeholder: counted, selectable, contained', () => {
  const {instance} = open({
    $schema: 'dg-ui/1',
    components: [{tag: 'u2-nope', name: 'broken'}, {tag: 'u2-button', name: 'visual'}],
    root: {tag: 'div', name: 'root'},
  });
  const [broken, visual] = instance.spec.components;
  assert.equal(brokenCount(instance), 2, 'both are placeholders, not exceptions');
  assert.equal(instance.nodes().has(broken), true, 'a designer must be able to select what it fixes');
  assert.match(instance.nodes().get(broken).textContent, /is not a registered component/);
  assert.match(instance.nodes().get(visual).textContent, /belongs in the form/);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0,
    'and nothing of the tray reaches the canvas');
  instance.dispose();
});

tray('a component warns about the children it drops and about having no name', () => {
  const warnings = captureWarnings(() => {
    const {instance} = open({
      $schema: 'dg-ui/1',
      components: [{tag: 'u2-fake-source', children: [{tag: 'div'}]}],
      root: {tag: 'div'},
    });
    instance.dispose();
  });
  assert.deepEqual(warnings, [
    'u2 spec: u2-fake-source: takes no children — 1 dropped',
    'u2 spec: u2-fake-source: a component without a name cannot be bound to',
  ]);
});

tray('construction is two-phase: a sub-bind resolves against an input declared after the source', () => {
  const {instance} = open({
    $schema: 'dg-ui/1',
    components: [{tag: 'u2-fake-source', name: 'src', bind: {'params.from': '$.later'}}],
    root: {tag: 'u2-form', name: 'form', children: [
      {tag: 'u2-text-input', name: 'later', props: {label: 'Later', value: 'ready'}},
    ]},
  });
  assert.equal(built.length, 1);
  assert.equal(built[0].started, 'ready', 'start() runs once the whole form is up');
  instance.dispose();
});

tray('add-component and remove-component: undo restores the exact dump and the same node', () => {
  const {instance, editor} = open({$schema: 'dg-ui/1', root: {tag: 'div', name: 'root'}});
  const before = instance.dump();
  const node = {tag: 'u2-state', name: 'draft'};

  editor.apply({op: 'add-component', index: 0, node});
  assert.equal(instance.node('draft') instanceof StateSource, true, 'it is built and named at once');
  assert.deepEqual(instance.dump().components, [node]);

  editor.undo();
  assert.deepEqual(instance.dump(), before, 'the emptied tray goes with its last entry');
  assert.equal(instance.nodes().has(node), false);
  assert.equal(instance.node('draft'), undefined);

  editor.redo();
  assert.equal(instance.spec.components[0] === node, true, 'the SAME node object comes back');

  editor.apply({op: 'remove-component', node});
  assert.deepEqual(instance.dump(), before);
  editor.undo();
  assert.equal(instance.spec.components[0] === node, true);
  assert.equal(instance.node('draft') instanceof StateSource, true);
  instance.dispose();
});

tray('canApply: a source belongs on the tray, a control in the form, an unknown tag in either', () => {
  const {instance, editor} = open(STATE_SPEC);
  const form = instance.spec.root;
  assert.match(editor.canApply({op: 'add', parent: form, index: 0, node: {tag: 'u2-state'}}),
    /belongs on the components tray/);
  assert.match(editor.canApply({op: 'add-component', index: 0, node: {tag: 'u2-button'}}),
    /belongs in the form/);
  assert.equal(editor.canApply({op: 'add-component', index: 0, node: {tag: 'u2-unknown'}}), null);
  assert.equal(editor.canApply({op: 'add', parent: form, index: 0, node: {tag: 'u2-unknown'}}), null);
  assert.match(editor.canApply({op: 'add-component', index: 9, node: {tag: 'u2-unknown'}}),
    /out of range for the tray/);
  assert.match(editor.canApply({op: 'remove', node: instance.spec.components[0]}),
    /the tray has its own add and remove/);
  instance.dispose();
});

tray('one namespace: a name is unique across both trees, and rename refuses a collision', () => {
  const {instance, editor} = open(STATE_SPEC);
  assert.equal(editor.uniqueName('draft'), 'draft1');
  assert.match(editor.canApply({op: 'rename', node: instance.spec.root, name: 'draft'}),
    /already taken/);
  instance.dispose();
});

tray('a patch on a component renders what references it again, and nothing else', () => {
  const {instance, editor} = open({
    $schema: 'dg-ui/1',
    components: [{tag: 'u2-fake-source', name: 'src', props: {text: 'one'}}],
    root: {tag: 'div', name: 'root', children: [
      {tag: 'u2-text-input', name: 'bound', props: {label: 'Bound'}, bind: {value: '$.src'}},
      {tag: 'u2-text-input', name: 'free', props: {label: 'Free'}},
    ]},
  });
  const [source] = instance.spec.components;
  const [boundNode, freeNode] = instance.spec.root.children;
  const boundBefore = instance.nodes().get(boundNode);
  const freeBefore = instance.nodes().get(freeNode);

  editor.apply({op: 'set-prop', node: source, name: 'text', value: 'two'});
  assert.equal(built.length, 2, 'the component itself is built again');
  assert.equal(instance.nodes().get(boundNode) === boundBefore, false,
    'the consumer captured the old signal and has to come back');
  assert.equal(instance.nodes().get(freeNode) === freeBefore, true,
    'a node that references nothing keeps its identity');
  assert.equal(editors(instance)[0].value, 'two', 'and it is wired to the live component');

  editor.undo();
  assert.equal(editors(instance)[0].value, 'one', 'undo rides the same path');
  instance.dispose();
});

tray('renaming a component rewrites its references and renders them again', () => {
  const {instance, editor} = open({
    $schema: 'dg-ui/1',
    components: [{tag: 'u2-fake-source', name: 'src', props: {text: 'one'}}],
    root: {tag: 'u2-text-input', name: 'bound', props: {label: 'Bound'}, bind: {value: '$.src'}},
  });
  const [source] = instance.spec.components;
  const before = instance.nodes().get(instance.spec.root);

  editor.apply({op: 'rename', node: source, name: 'orders'});
  assert.equal(instance.spec.root.bind.value, '$.orders');
  assert.equal(instance.nodes().get(instance.spec.root) === before, false);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0,
    'the consumer resolves against the new name');

  editor.undo();
  assert.equal(instance.spec.root.bind.value, '$.src');
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  instance.dispose();
});

tray('setDesignTime re-mounts every component, re-renders its consumers and keeps the history', () => {
  const {instance, editor} = open({
    $schema: 'dg-ui/1',
    components: [{tag: 'u2-fake-source', name: 'src', props: {text: 'one'}}],
    root: {tag: 'u2-text-input', name: 'bound', props: {label: 'Bound'}, bind: {value: '$.src'}},
  });
  editor.apply({op: 'set-prop', node: instance.spec.components[0], name: 'text', value: 'two'});
  const mounts = built.length;
  const consumer = instance.nodes().get(instance.spec.root);
  assert.equal(instance.designTime, false);
  assert.equal(built[mounts - 1].designTime, false);

  instance.setDesignTime(true);
  assert.equal(instance.designTime, true);
  assert.equal(built.length, mounts + 1, 'the source is built again for the new mode');
  assert.equal(built[mounts].designTime, true);
  assert.equal(instance.nodes().get(instance.spec.root) === consumer, false,
    'and its consumer with it — it held the old component\'s signal');

  instance.setDesignTime(true);
  assert.equal(built.length, mounts + 1, 'the same mode twice is a no-op');

  assert.equal(editor.canUndo.value, true, 'the mode is view state: the history survives it');
  editor.undo();
  assert.equal(instance.dump().components[0].props.text, 'one');
  instance.dispose();
});

/* A source whose param binds to ANOTHER source is a dependent that is itself a component: the
   toggle has already re-mounted it, so re-rendering it as a dependent would build it a second time
   and start it twice — and `start()` is not idempotent, so every effect it registers, and every
   request that follows, would double for the rest of the session. */
tray('a source bound to another source is started exactly once per toggle', () => {
  const {instance} = open({
    $schema: 'dg-ui/1',
    components: [
      {tag: 'u2-fake-source', name: 'a', props: {text: 'one'}},
      {tag: 'u2-fake-source', name: 'b', bind: {'params.from': '$.a.value'}},
    ],
    root: {tag: 'u2-text-input', name: 'bound', props: {label: 'Bound'}, bind: {value: '$.b'}},
  });
  assert.deepEqual(built.map((source) => source.starts), [1, 1], 'construction starts each once');

  instance.setDesignTime(true);
  assert.equal(built.length, 4, 'each component is re-mounted exactly once');
  assert.deepEqual(built.slice(2).map((source) => source.starts), [1, 1]);
  assert.equal(built[3].bound, built[2].value, 'and the second still resolves against the first');
  instance.dispose();
});

/* The panel builds its inputs inside its own scope, and an input committing on the way makes the
   patch land with that scope ambient. A tray component owned by it dies with the next
   `disposePanel()` — with the form still on screen and bound to it. */
tray('a component mounted while another scope is ambient belongs to the instance', () => {
  const {instance, editor} = open(STATE_SPEC);
  const panel = new Scope();
  Scope.runWith(panel, () =>
    editor.apply({op: 'set-prop', node: instance.spec.components[0], name: 'initial', value: 'two'}));
  panel.dispose();

  const state = instance.node('draft');
  assert.equal(state.scope.isDisposed, false, 'the panel that happened to be building took it down');
  state.value.value = 'from the source';
  assert.equal(editors(instance)[0].value, 'from the source', 'and the form still follows it');
  instance.dispose();
});

/* WO-16 — the crux of the two phases, on the way back: a re-mounted source must resolve its bound
   params against the form the SAME change leaves behind. A dependent's re-render promotes to the
   container that adopts it, which rebuilds every input under it — a source started before that
   would hold the signal of an input nothing shows any more. */
tray('a re-mounted source starts against the form the dependent re-render leaves behind', () => {
  const spec = {$schema: 'dg-ui/1',
    components: [{tag: 'u2-fake-source', name: 'src', bind: {'params.from': '$.field.value'}}],
    root: {tag: 'u2-form', name: 'form', children: [
      {tag: 'u2-text-input', name: 'field', props: {label: 'Field', value: 'seed'}},
      {tag: 'u2-text-input', name: 'echo', props: {label: 'Echo'}, bind: {value: '$.src.value'}},
    ]}};
  const {instance, editor} = open(spec);
  const field = () => instance.nodes().get(instance.spec.root.children[0]).value;
  assert.equal(built[0].started, 'seed');
  assert.equal(built[0].bound, field(), 'construction resolves against the form it just built');

  instance.setDesignTime(true);
  assert.equal(built.length, 2, 'the toggle re-mounts the source');
  assert.equal(built[1].bound, field(), 'and the toggle re-resolves against the rebuilt input');

  editor.apply({op: 'set-prop', node: instance.spec.components[0], name: 'text', value: 'two'});
  assert.equal(built.length, 3);
  assert.equal(built[2].bound, field(), 'as does a prop edit that re-renders the same dependents');
  instance.dispose();
});

tray('the structure tree carries the form first and the tray after it', () => {
  const {instance} = open(STATE_SPEC);
  const {roots, ids} = specTree(instance);
  assert.deepEqual(roots.map((r) => r.label), ['form', 'draft']);
  assert.equal(ids.get(instance.spec.components[0]), '1');
  assert.equal(roots[1].children, undefined);
  instance.dispose();
});
