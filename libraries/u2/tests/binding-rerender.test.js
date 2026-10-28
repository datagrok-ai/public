/* Binds never rebuild at run time (R7, P3, P4): a bind on a declared non-bindable prop builds with
   the resolved value and then writes every change through the component's own setter — the node,
   its siblings and its render root stay. A structural prop reads its bind once. Rebuilding is left
   to the designer's patches, which promote to the render root and drag their dependents along. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Signal, batch} from '../src/core/signals.js';
import {Scope} from '../src/core/scope.js';
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

/** A fake input whose creations are counted — the rebuild meter every assert reads. `label` is
 * declared non-bindable, so a bind on it goes through TextInput's `label` setter. Capped so a
 * runaway rebuild fails the test instead of blowing the stack. */
function registerCounted(reg, counts) {
  reg.register({
    tag: 'u2-fake-input',
    description: 'Counted fake input for the live-bind tests',
    props: [
      {name: 'label', type: 'string'},
      {name: 'value', type: 'string', bindable: true, twoWay: true},
      {name: 'mode', type: 'string', structural: true},
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

async function captureWarnings(body) {
  const warnings = [];
  const original = console.warn;
  console.warn = (message) => warnings.push(message);
  try {
    await body();
  } finally {
    console.warn = original;
  }
  return warnings;
}

const labelOf = (instance, name) => instance.node(name).root.querySelector('.u2-input-label').textContent;

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

spec('a source change writes through the setter at once: same node, same element, no rebuild', async () => {
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
  const before = instance.node('a');

  ctx.data.s.value = 'Two';
  assert.equal(labelOf(instance, 'a'), 'Two', 'synchronous — no microtask, no queue');
  ctx.data.s.value = 'Three';
  await flush();
  assert.equal(counts.length, 1, 'never re-created');
  assert.equal(instance.node('a'), before, 'the same component');
  assert.equal(labelOf(instance, 'a'), 'Three');
  instance.dispose();
});

spec('batched writes to a source shared by two nodes under an adopting root: no rebuild, the last value lands',
  async () => {
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
    const form = instance.node('form').root;
    batch(() => {
      ctx.data.s.value = 'Two';
      ctx.data.s.value = 'Three';
      ctx.data.s.value = 'Four';
    });
    await flush();
    assert.equal(counts.length, 2, 'each child built once, ever');
    assert.equal(instance.node('form').root, form, 'no promotion to the render root');
    assert.equal(labelOf(instance, 'a'), 'Four');
    assert.equal(labelOf(instance, 'b'), 'Four');
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
  assert.equal(counts.length, built, 'the follower died with the node');
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  instance.dispose();
});

spec('a node followed by others is never re-created by its own binds, so its dependents stay bound', async () => {
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
  ctx.data.s.value = 'Two';
  await flush();
  assert.equal(counts.length, 2, 'nothing re-created');
  instance.node('a').value = 'typed';
  assert.equal(instance.node('b').root.querySelector('input').value, 'typed', 'b still follows a');
  instance.dispose();
});

spec('forward ref: a setter bind to a later node applies at the flush, built once', async () => {
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
  assert.equal(labelOf(instance, 'early'), 'target');
  assert.equal(counts.length, 2, 'no resolving render');
  instance.node('later').value = 'typed';
  assert.equal(labelOf(instance, 'early'), 'typed', 'and follows from then on');
  instance.dispose();
});

spec('forward ref under an adopting render root: resolved at the flush, no promotion, no loop', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-form', name: 'form', children: [
      {tag: 'u2-fake-input', name: 'early', bind: {label: '$.later'}},
      {tag: 'u2-list-input', name: 'list'},
      {tag: 'u2-fake-input', name: 'later', props: {value: 'target'}},
    ]},
  }, new SpecContext(), reg);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0, 'built, never broken');
  assert.equal(counts.length, 2);
  assert.equal(labelOf(instance, 'early'), 'target', 'the value resolved');
  instance.node('later').value = 'typed';
  await flush();
  await flush();
  assert.equal(counts.length, 2, 'a later source change re-creates nothing');
  assert.equal(labelOf(instance, 'early'), 'typed');
  instance.dispose();
});

spec('forward ref: a target that never builds warns once and the node stays rendered', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  let instance;
  const warnings = await captureWarnings(() => {
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
  assert.equal(labelOf(instance, 'early'), 'Literal');
  await flush();
  instance.dispose();
});

spec('a structural prop reads its bind once and says so; a later change leaves the node alone', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const ctx = new SpecContext({data: {m: 'compact'}});
  let instance;
  const warnings = await captureWarnings(async () => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'div', children: [{tag: 'u2-fake-input', name: 'a', bind: {mode: '$.m'}}]},
    }, ctx, reg);
    ctx.data.m.value = 'wide';
    await flush();
  });
  assert.equal(counts.length, 1);
  assert.equal(counts[0].mode, 'compact', 'built with the value at render');
  assert.deepEqual(warnings, ['u2 spec: u2-fake-input.mode is structural — a bind is read once']);
  instance.dispose();
});

spec('childProps: a bound pane title renders resolved, follows the source in place, unbinds to the literal',
  async () => {
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
    const strip = instance.node('tabs');
    const header = strip.root.querySelector('.u2-tabs-tab');

    ctx.data.t.value = 'World';
    assert.equal(label(), 'World', 'the tab strip renamed the tab');
    assert.equal(instance.node('tabs'), strip, 'the strip was not rebuilt');
    assert.equal(strip.root.querySelector('.u2-tabs-tab'), header, 'nor the header');

    const editor = new SpecEditor(instance);
    editor.apply({op: 'set-bind', node: instance.spec.root.children[0], name: 'title', path: undefined});
    assert.equal(label(), 'Literal', 'unbound again (a designer patch), the literal is back');
    instance.dispose();
  });

spec('a bound tab title keeps the user\'s tab; a designer patch carries it unless activeTab is declared',
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
    assert.equal(instance.node('tabs').activeTab.peek(), 'tab-1', 'nothing rebuilt, nothing to lose');

    const editor = new SpecEditor(instance);
    editor.apply({op: 'set-prop', node: instance.spec.root.children[1], name: 'title', value: 'Renamed'});
    const carried = instance.node('tabs');
    assert.equal(carried.root.querySelectorAll('.u2-tabs-label')[1].textContent, 'Renamed',
      'the patch rebuilt the strip');
    assert.equal(carried.activeTab.peek(), 'tab-1', 'with the user\'s tab still active');
    instance.dispose();

    const pinned = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'u2-tabs', name: 'tabs', props: {activeTab: 'tab-0'}, children: [
        {tag: 'u2-panel', name: 'first', props: {title: 'Literal'}},
        {tag: 'u2-panel', name: 'second', props: {title: 'Second'}},
      ]},
    }, ctx, reg);
    document.body.append(pinned.root);
    pinned.node('tabs').activeTab.value = 'tab-1';
    new SpecEditor(pinned).apply({op: 'set-prop', node: pinned.spec.root.children[1], name: 'title', value: 'X'});
    assert.equal(pinned.node('tabs').activeTab.peek(), 'tab-0',
      'a declared activeTab stays the document\'s to decide');
    pinned.dispose();
  });

spec('tray config: a bound structural source prop is read once — the source is never re-mounted', async () => {
  const reg = registry();
  const counts = [];
  registerCounted(reg, counts);
  const ctx = new SpecContext({data: {seed: 'one'}});
  let instance;
  const warnings = await captureWarnings(async () => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      components: [{tag: 'u2-state', name: 'b', props: {type: 'string'}, bind: {initial: '$.seed'}}],
      root: {tag: 'div', children: [{tag: 'u2-fake-input', name: 'in', bind: {value: '$.b'}}]},
    }, ctx, reg);
    assert.equal(instance.node('in').root.querySelector('input').value, 'one');
    const source = instance.nodes().get(instance.spec.components[0]);
    ctx.data.seed.value = 'two';
    await flush();
    assert.equal(instance.nodes().get(instance.spec.components[0]), source, 'not re-mounted');
  });
  assert.equal(instance.node('in').root.querySelector('input').value, 'one');
  assert.equal(counts.length, 1);
  assert.equal(warnings.filter((w) => w.includes('u2-state.initial is structural')).length, 1);
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
  instance.node('first').value = 'typed';
  await flush();
  assert.equal(title(), 'typed', 'the heading follows the re-created input, not the corpse');
  instance.dispose();
});

spec('a patch promoted to an adopting root drags outside setter binds with it', async () => {
  const reg = registry();
  const instance = renderSpec(
    personaSpec({tag: 'u2-button', name: 'go', bind: {text: '$.first.value'}}),
    new SpecContext(), reg);
  document.body.append(instance.root);

  const editor = new SpecEditor(instance);
  const second = instance.spec.root.children[1].children[1];
  editor.apply({op: 'set-bind', node: second, name: 'label', path: '$.first.value'});
  instance.node('first').value = 'typed';
  await flush();
  assert.equal(instance.node('go').root.textContent, 'typed',
    'the button follows the re-created input');
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
  instance.node('first').value = 'typed';
  await flush();
  editor.apply({op: 'set-bind', node: titleNode, name: 'text', path: undefined});
  await flush();
  assert.equal(title(), 'Static', 'no stale bound value survives the unbind');
  instance.node('first').value = 'more';
  await flush();
  assert.equal(title(), 'Static');

  editor.apply({op: 'set-bind', node: titleNode, name: 'text', path: '$.first.value'});
  assert.equal(title(), 'more', 'a rebind resolves the current value immediately');
  instance.node('first').value = 'live';
  await flush();
  assert.equal(title(), 'live');
  instance.dispose();
});

spec('the dependent expansion of a designer patch chases a chain across render roots to a fixed point',
  async () => {
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

    instance.node('c').value = 'c1';
    await flush();
    assert.equal(labelOf(instance, 'b'), 'c1', 'b follows the re-created c');
    instance.node('b').value = 'b1';
    await flush();
    assert.equal(labelOf(instance, 'a'), 'b1', 'a follows the re-created b');
    assert.equal(counts.length, 8, 'and following re-created nothing');
    instance.dispose();
  });
