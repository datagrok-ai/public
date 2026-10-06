/* The re-render binding tier (UB-1/UB-2): a bind on a declared non-bindable prop participates as
   the literal at build, and a source change re-creates the node at its render root — coalesced
   into one microtask, dependents included. Forward refs re-render once at the flush; childProps
   and tray config ride the same machinery. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Signal, batch} from '../src/core/signals.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {TextInput} from '../src/components/inputs/text-input.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {SpecEditor} from '../src/spec/editor.js';
import {registerAll} from '../src/spec/registrations.js';

/** Every test runs against a clean document and must leave the live-scope count where it was. */
function spec(name, body) {
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

function registry() {
  const reg = new Registry();
  registerAll(reg);
  return reg;
}

/** A fake input whose creations are counted — the rebuild meter every coalescing assert reads.
 * Capped so a runaway re-render loop fails the test instead of blowing the stack. */
function registerCounted(reg, counts) {
  reg.register({
    tag: 'u2-fake-input',
    description: 'Counted fake input for the re-render tests',
    props: [
      {name: 'label', type: 'string'},
      {name: 'value', type: 'string', bindable: true, twoWay: true},
    ],
    create: (props) => {
      counts.push({...props});
      if (counts.length > 20)
        throw new Error('runaway rebuild');
      const bound = props.value instanceof Signal;
      return new TextInput({
        label: typeof props.label === 'string' ? props.label : undefined,
        value: bound ? undefined : props.value,
        bind: bound ? props.value : undefined,
      });
    },
    example: {tag: 'u2-fake-input'},
  });
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

spec('a bind on a declared non-bindable prop builds with the resolved value', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {s: 'Save all'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [{tag: 'u2-button', name: 'go', bind: {text: '$.s'}}]},
  }, ctx, reg);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  assert.equal(instance.node('go').root.textContent, 'Save all');
  instance.dispose();
});

spec('a bind on an undeclared prop stays an error, with the new message', () => {
  const reg = registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [{tag: 'u2-button', bind: {ghost: '$.s'}}]},
  }, new SpecContext({data: {s: 'x'}}), reg);
  const errors = instance.root.querySelectorAll('.u2-spec-error');
  assert.equal(errors.length, 1);
  assert.match(errors[0].textContent, /u2-button: has no prop "ghost" to bind/);
  instance.dispose();
});

spec('a source change re-creates the node once per microtask, element identity changed', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const ctx = new SpecContext({data: {s: 'One'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [{tag: 'u2-fake-input', name: 'a', bind: {label: '$.s'}}]},
  }, ctx, reg);
  document.body.append(instance.root);
  assert.equal(counts.length, 1);
  assert.equal(counts[0].label, 'One');
  const before = instance.node('a').root;

  ctx.data.s.value = 'Two';
  ctx.data.s.value = 'Three';
  assert.equal(counts.length, 1, 'nothing rebuilds synchronously — the effect only enqueues');
  await flush();
  assert.equal(counts.length, 2, 'two writes, one microtask, one rebuild');
  assert.equal(counts[1].label, 'Three');
  assert.notEqual(instance.node('a').root, before, 'the node was re-created');
  assert.equal(instance.root.querySelector('.u2-input-label').textContent, 'Three');
  instance.dispose();
});

spec('three batched writes to a source shared by two nodes under one render root: one rebuild', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const ctx = new SpecContext({data: {s: 'One'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-form', name: 'form', children: [
      {tag: 'u2-fake-input', name: 'a', bind: {label: '$.s'}},
      {tag: 'u2-fake-input', name: 'b', bind: {label: '$.s'}},
    ]},
  }, ctx, reg);
  document.body.append(instance.root);
  assert.equal(counts.length, 2, 'two children built once each');
  const before = instance.node('form').root;

  batch(() => {
    ctx.data.s.value = 'Two';
    ctx.data.s.value = 'Three';
    ctx.data.s.value = 'Four';
  });
  assert.equal(counts.length, 2, 'nothing rebuilds synchronously — the effect only enqueues');
  await flush();
  assert.equal(counts.length, 4, 'one rebuild of the shared render root — each child re-created once');
  assert.equal(counts[2].label, 'Four');
  assert.equal(counts[3].label, 'Four');
  assert.notEqual(instance.node('form').root, before, 'the root itself was re-created (promotion)');
  await flush();
  await flush();
  assert.equal(counts.length, 4, 'and the queue is quiet afterwards');
  instance.dispose();
});

spec('a removed node ignores later source writes', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const ctx = new SpecContext({data: {s: 'One'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [{tag: 'u2-fake-input', name: 'a', bind: {label: '$.s'}}]},
  }, ctx, reg);
  const editor = new SpecEditor(instance);
  editor.apply({op: 'remove', node: instance.spec.root.children?.[0] ?? instance.spec.root});
  const built = counts.length;
  ctx.data.s.value = 'Two';
  await flush();
  assert.equal(counts.length, built, 'the disposed effect and the queue guard both stay silent');
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  instance.dispose();
});

spec('a named re-render-bound node drags its dependents with it', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const ctx = new SpecContext({data: {s: 'One'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'u2-fake-input', name: 'a', props: {value: 'seed'}, bind: {label: '$.s'}},
      {tag: 'u2-fake-input', name: 'b', bind: {value: '$.a'}},
    ]},
  }, ctx, reg);
  assert.equal(counts.length, 2);
  ctx.data.s.value = 'Two';
  await flush();
  assert.equal(counts.length, 4, 'the dependent re-rendered with the node it references');
  const b = instance.node('b');
  assert.equal(b.root.querySelector('input').value, 'seed',
    'the dependent follows the rebuilt node\'s signal, not the corpse');
  instance.dispose();
});

spec('re-renders never write signals, so the queue drains without looping', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const ctx = new SpecContext({data: {s: 'One'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [{tag: 'u2-fake-input', name: 'a', bind: {label: '$.s'}}]},
  }, ctx, reg);
  ctx.data.s.value = 'Two';
  await flush();
  await flush();
  await flush();
  assert.equal(counts.length, 2, 'one change, one rebuild, then quiet');
  instance.dispose();
});

spec('forward ref: a re-render bind to a later node resolves through one extra render', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'u2-fake-input', name: 'early', bind: {label: '$.later'}},
      {tag: 'u2-fake-input', name: 'later', props: {value: 'target'}},
    ]},
  }, new SpecContext(), reg);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0, 'built, never broken');
  assert.equal(instance.node('early').root.querySelector('.u2-input-label').textContent, 'target');
  const built = counts.length;
  await flush();
  await flush();
  assert.equal(counts.length, built, 'the second render resolved directly — no further re-renders');
  instance.dispose();
});

spec('forward ref under an adopting render root terminates with the value resolved', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-form', name: 'form', children: [
      {tag: 'u2-fake-input', name: 'early', bind: {label: '$.later'}},
      {tag: 'u2-fake-input', name: 'later', props: {value: 'target'}},
    ]},
  }, new SpecContext(), reg);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0, 'built, never broken');
  assert.equal(counts.length, 4,
    'one resolving rebuild of the promoted root, then quiet — the re-created target re-defers no loop');
  assert.equal(instance.node('early').root.querySelector('.u2-input-label').textContent, 'target',
    'the value resolved');
  await flush();
  await flush();
  assert.equal(counts.length, 4, 'no leaked deferrals — the queue is quiet');

  instance.node('later').value.value = 'typed';
  await flush();
  await flush();
  await flush();
  assert.equal(counts.length, 8,
    'a later source change converges too: rebuild, re-bake against the re-created target, wire');
  instance.dispose();
});

spec('forward ref under an adopting root to an identity-unstable value still terminates', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-form', name: 'form', children: [
      {tag: 'u2-fake-input', name: 'early', bind: {label: '$.later'}},
      // a list input's value is a FRESH array on every build — Object.is can never call it converged
      {tag: 'u2-list-input', name: 'later'},
    ]},
  }, new SpecContext(), reg);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0, 'built, never broken');
  assert.equal(counts.length, 3, 'the round cap stops the re-bake and wires the follow');
  await flush();
  await flush();
  assert.equal(counts.length, 3, 'quiet afterwards');
  instance.dispose();
});

spec('forward ref: a target that never builds warns once and the node stays rendered', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  let instance;
  const warnings = captureWarnings(() => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'div', children: [
        {tag: 'u2-fake-input', name: 'early', props: {label: 'Literal'}, bind: {label: '$.ghost'}},
        // ghost is declared but dropped: the parent takes no children
        {tag: 'u2-fake-input', children: [{tag: 'u2-fake-input', name: 'ghost'}]},
      ]},
    }, new SpecContext(), reg);
  });
  assert.equal(warnings.filter((w) => w.includes('declared but not built')).length, 1);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  assert.equal(instance.node('early').root.querySelector('.u2-input-label').textContent, 'Literal');
  await flush();
  instance.dispose();
});

spec('childProps: a bound pane title renders resolved, follows the source, unbinds to the literal', async () => {
  const reg = registry();
  const ctx = new SpecContext({data: {t: 'Hello'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-tabs', name: 'tabs', children: [
      {tag: 'u2-panel', name: 'pane', props: {title: 'Literal'}, bind: {title: '$.t'}},
    ]},
  }, ctx, reg);
  document.body.append(instance.root);
  const label = () => instance.root.querySelector('.u2-tabs-label').textContent;
  assert.equal(label(), 'Hello', 'the resolved bind wins over the literal');

  ctx.data.t.value = 'World';
  await flush();
  assert.equal(label(), 'World', 'the tabs parent re-rendered with the new label');

  const editor = new SpecEditor(instance);
  editor.apply({op: 'set-bind', node: instance.spec.root.children[0], name: 'title', path: undefined});
  assert.equal(label(), 'Literal', 'unbound again, the literal is back');
  instance.dispose();
});

spec('a childProp re-render keeps the active tab while the pane exists; an explicit activeTab wins',
  async () => {
    const reg = registry();
    const ctx = new SpecContext({data: {t: 'Hello'}});
    const instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'u2-tabs', name: 'tabs', children: [
        {tag: 'u2-panel', name: 'first', props: {title: 'Literal'}, bind: {title: '$.t'}},
        {tag: 'u2-panel', name: 'second', props: {title: 'Second'}},
      ]},
    }, ctx, reg);
    document.body.append(instance.root);
    instance.node('tabs').activeTab.value = 'tab-1';
    ctx.data.t.value = 'World';
    await flush();
    const carried = instance.node('tabs');
    assert.equal(carried.root.querySelector('.u2-tabs-label').textContent, 'World',
      'the strip re-rendered onto the new title');
    assert.equal(carried.activeTab.peek(), 'tab-1', 'with the user\'s tab still active');
    instance.dispose();

    const pinned = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'u2-tabs', name: 'tabs', props: {activeTab: 'tab-0'}, children: [
        {tag: 'u2-panel', name: 'first', props: {title: 'Literal'}, bind: {title: '$.t'}},
        {tag: 'u2-panel', name: 'second', props: {title: 'Second'}},
      ]},
    }, ctx, reg);
    document.body.append(pinned.root);
    pinned.node('tabs').activeTab.value = 'tab-1';
    ctx.data.t.value = 'Again';
    await flush();
    assert.equal(pinned.node('tabs').activeTab.peek(), 'tab-0',
      'a declared activeTab stays the document\'s to decide');
    pinned.dispose();
  });

spec('tray config: a bound source prop re-mounts the source and its dependents follow', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const ctx = new SpecContext({data: {seed: 'one'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    components: [{tag: 'u2-state', name: 'b', props: {type: 'string'}, bind: {initial: '$.seed'}}],
    root: {tag: 'div', children: [{tag: 'u2-fake-input', name: 'in', bind: {value: '$.b'}}]},
  }, ctx, reg);
  assert.equal(instance.node('in').root.querySelector('input').value, 'one');

  ctx.data.seed.value = 'two';
  await flush();
  assert.equal(instance.node('in').root.querySelector('input').value, 'two',
    'the re-mounted source carries the new initial and the input re-rendered onto it');
  instance.dispose();
});

spec('a tray source and a visual dependent queued in one microtask: the source re-mounts first', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const ctx = new SpecContext({data: {seed: 'one', s2: 'A'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    components: [{tag: 'u2-state', name: 'st', props: {type: 'string'}, bind: {initial: '$.seed'}}],
    root: {tag: 'div', children: [
      {tag: 'u2-fake-input', name: 'in', bind: {value: '$.st', label: '$.s2'}},
    ]},
  }, ctx, reg);
  // the dependent enqueues first: its own label bind fires before the tray's config bind
  ctx.data.s2.value = 'B';
  ctx.data.seed.value = 'two';
  await flush();
  assert.equal(instance.node('in').root.querySelector('input').value, 'two');
  instance.resolveBinding('$.st').signal.value = 'live';
  assert.equal(instance.node('in').root.querySelector('input').value, 'live',
    'the dependent ends bound to the live re-mounted source, not the corpse');
  instance.dispose();
});

spec('tray: a junk bind key breaks that chip alone, with the standard placeholder', () => {
  const reg = registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    components: [{tag: 'u2-state', name: 's', bind: {bogus: '$.x'}}],
    root: {tag: 'div'},
  }, new SpecContext(), reg);
  const built = instance.nodes().get(instance.spec.components[0]);
  assert.equal(built.classList.contains('u2-spec-error'), true);
  assert.match(built.textContent, /has no prop "bogus" to bind/);
  instance.dispose();
});

spec('tray: a dotted bind whose head takes no sub-binds breaks the chip', () => {
  const reg = registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    components: [{tag: 'u2-state', name: 's', bind: {'initial.x': '$.y'}}],
    root: {tag: 'div'},
  }, new SpecContext(), reg);
  const built = instance.nodes().get(instance.spec.components[0]);
  assert.equal(built.classList.contains('u2-spec-error'), true);
  assert.match(built.textContent, /prop "initial" does not take sub-binds/);
  instance.dispose();
});

/** The persona spec: a heading outside a form, live-bound into it — the shape every
 * patch-promotion test below rebuilds. */
function personaSpec(outside) {
  return {
    $schema: 'dg-ui/1',
    root: {tag: 'u2-panel', children: [
      outside,
      {tag: 'u2-form', name: 'form', children: [
        {tag: 'u2-text-input', name: 'first', props: {label: 'First', value: 'seed'}},
        {tag: 'u2-text-input', name: 'second', props: {label: 'Second'}},
      ]},
    ]},
  };
}

spec('a patch promoted to an adopting root drags outside HTML live binds with it', async () => {
  const reg = registry();
  const instance = renderSpec(
    personaSpec({tag: 'h3', name: 'title', bind: {text: '$.first.value'}}),
    new SpecContext(), reg);
  document.body.append(instance.root);
  const title = () => instance.root.querySelector('h3').textContent;
  assert.equal(title(), 'seed');

  const editor = new SpecEditor(instance);
  const second = instance.spec.root.children[1].children[1];
  editor.apply({op: 'set-bind', node: second, name: 'label', path: '$.first.value'});
  instance.node('first').value.value = 'typed';
  await flush();
  assert.equal(title(), 'typed', 'the heading follows the re-created input, not the corpse');
  instance.dispose();
});

spec('a patch promoted to an adopting root drags outside re-render binds with it', async () => {
  const reg = registry();
  const instance = renderSpec(
    personaSpec({tag: 'u2-button', name: 'go', bind: {text: '$.first.value'}}),
    new SpecContext(), reg);
  document.body.append(instance.root);

  const editor = new SpecEditor(instance);
  const second = instance.spec.root.children[1].children[1];
  editor.apply({op: 'set-bind', node: second, name: 'label', path: '$.first.value'});
  instance.node('first').value.value = 'typed';
  await flush();
  assert.equal(instance.node('go').root.textContent, 'typed',
    'the button re-rendered against the re-created input');
  instance.dispose();
});

spec('unbind after the promoted rebuild: the literal is back at once, and typing stops arriving', async () => {
  const reg = registry();
  const heading = {tag: 'h3', name: 'title', props: {text: 'Static'}, bind: {text: '$.first.value'}};
  const instance = renderSpec(personaSpec(heading), new SpecContext(), reg);
  document.body.append(instance.root);
  const title = () => instance.root.querySelector('h3').textContent;

  const editor = new SpecEditor(instance);
  const [titleNode, form] = instance.spec.root.children;
  editor.apply({op: 'set-bind', node: form.children[1], name: 'label', path: '$.first.value'});
  instance.node('first').value.value = 'typed';
  await flush();
  editor.apply({op: 'set-bind', node: titleNode, name: 'text', path: undefined});
  await flush();
  assert.equal(title(), 'Static', 'no stale bound value survives the unbind');
  instance.node('first').value.value = 'more';
  await flush();
  assert.equal(title(), 'Static');

  editor.apply({op: 'set-bind', node: titleNode, name: 'text', path: '$.first.value'});
  assert.equal(title(), 'more', 'a rebind resolves the current value immediately');
  instance.node('first').value.value = 'live';
  await flush();
  assert.equal(title(), 'live');
  instance.dispose();
});

spec('the dependent expansion chases a chain across render roots to a fixed point', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-panel', children: [
      {tag: 'u2-form', name: 'formC', children: [
        {tag: 'u2-fake-input', name: 'c', props: {value: 'c0'}},
        {tag: 'u2-fake-input', name: 'sib'},
      ]},
      {tag: 'u2-form', name: 'formB', children: [
        {tag: 'u2-fake-input', name: 'b', props: {value: 'b0'}, bind: {label: '$.c.value'}},
      ]},
      {tag: 'u2-fake-input', name: 'a', bind: {label: '$.b.value'}},
    ]},
  }, new SpecContext(), reg);
  document.body.append(instance.root);
  assert.equal(counts.length, 4);

  const editor = new SpecEditor(instance);
  const sib = instance.spec.root.children[0].children[1];
  editor.apply({op: 'set-prop', node: sib, name: 'label', value: 'X'});
  assert.equal(counts.length, 8, 'c and sib with their root, then b, then a — each exactly once');

  instance.node('c').value.value = 'c1';
  await flush();
  assert.equal(instance.node('b').root.querySelector('.u2-input-label').textContent, 'c1',
    'b watches the re-created c');
  instance.node('b').value.value = 'b1';
  await flush();
  assert.equal(instance.node('a').root.querySelector('.u2-input-label').textContent, 'b1',
    'a watches the re-created b');
  instance.dispose();
});
