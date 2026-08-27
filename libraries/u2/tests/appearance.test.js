/* The shared Appearance prop group: injected at registration for every visual tag, applied to the
   root style once after create (literals) or through an effect (bound signals), filtered back out
   of the per-component manifest entries, and routed through the HTML-node path. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, htmlProps, renderSpec} from '../src/spec/spec.js';
import {SpecEditor} from '../src/spec/editor.js';
import {registerAll} from '../src/spec/registrations.js';
import {APPEARANCE_CATEGORY, APPEARANCE_PROPS} from '../src/spec/appearance.js';

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

spec('injection: every visual meta answers the shared group, a non-visual one none', () => {
  const reg = registry();
  for (const meta of reg.metas()) {
    const injected = meta.props.filter((p) => p.category === APPEARANCE_CATEGORY);
    if (meta.visual === false) {
      assert.equal(injected.length, 0, `${meta.tag} is non-visual and got appearance props`);
      continue;
    }
    assert.ok(injected.length >= 12, `${meta.tag} misses part of the appearance group`);
    for (const p of injected) {
      assert.equal(p.bindable, true, `${meta.tag}.${p.name} must be bindable`);
      assert.equal(p.twoWay, undefined, `${meta.tag}.${p.name} must never be twoWay`);
    }
  }
  assert.equal(reg.get('u2-state').props.some((p) => p.category === APPEARANCE_CATEGORY), false);
});

spec('collision: a component-own prop wins, the shared one is skipped and never applied', () => {
  const reg = new Registry();
  const own = {name: 'width', type: 'string', description: 'Own width'};
  reg.register({
    tag: 'u2-fake-box',
    description: 'Fake box for the collision test',
    props: [own],
    create: () => new Control(),
    example: {tag: 'u2-fake-box'},
  });
  const meta = reg.get('u2-fake-box');
  const injected = meta.props.filter((p) => p.category === APPEARANCE_CATEGORY).map((p) => p.name);
  assert.equal(injected.length, 12);
  assert.equal(injected.includes('width'), false);
  assert.equal(meta.props.find((p) => p.name === 'width'), own, 'the own meta object survives');

  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-box', name: 'box', props: {width: '200px', color: 'red'}},
  }, new SpecContext(), reg);
  const box = instance.node('box');
  assert.equal(box.root.style.color, 'red');
  assert.equal(box.root.style.width, undefined, 'the applier never touches the colliding prop');
  instance.dispose();
});

spec('appearance: false keeps the group off the tag entirely', () => {
  const reg = new Registry();
  reg.register({
    tag: 'u2-fake-viewer',
    appearance: false,
    description: 'Fake viewer for the opt-out test',
    props: [{name: 'table', type: 'string'}],
    create: () => new Control(),
    example: {tag: 'u2-fake-viewer'},
  });
  const meta = reg.get('u2-fake-viewer');
  assert.deepEqual(meta.props.map((p) => p.name), ['table']);
});

spec('literal application: color, opacity and font land on the root style, warning-free', () => {
  const reg = registry();
  let instance;
  const warnings = captureWarnings(() => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'u2-text-input', name: 'name', props: {
        label: 'Name', color: 'var(--dg-accent)', opacity: 0.5, font: 'bold 14px Roboto'}},
    }, new SpecContext(), reg);
  });
  assert.deepEqual(warnings, []);
  const el = instance.node('name').root;
  assert.equal(el.style.color, 'var(--dg-accent)');
  assert.equal(el.style.opacity, '0.5');
  assert.equal(el.style.font, 'bold 14px Roboto');
  instance.dispose();
});

spec('bound appearance: the style follows the context signal, and the effect dies with the node', async () => {
  const reg = registry();
  const ctx = new SpecContext({data: {c: 'red'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-text-input', name: 'name', props: {label: 'Name'}, bind: {color: '$.c'}},
  }, ctx, reg);
  const el = instance.node('name').root;
  assert.equal(el.style.color, 'red');
  ctx.data.c.value = 'blue';
  await flush();
  assert.equal(el.style.color, 'blue');
  instance.dispose();
  ctx.data.c.value = 'green';
  await flush();
  assert.equal(el.style.color, 'blue', 'a disposed node stops following');
});

spec('forward reference: color bound to a later input restyles when the input is edited', async () => {
  const reg = registry();
  let instance;
  const warnings = captureWarnings(() => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'div', children: [
        {tag: 'u2-text-input', name: 'styled', props: {label: 'Styled'}, bind: {color: '$.in1.value'}},
        {tag: 'u2-text-input', name: 'in1', props: {label: 'Color', value: 'red'}},
      ]},
    }, new SpecContext(), reg);
  });
  assert.deepEqual(warnings, []);
  const el = instance.node('styled').root;
  assert.equal(el.style.color, 'red', 'seeded from the later node once the pass is over');
  type(instance.node('in1').root.querySelector('input'), 'blue');
  await flush();
  assert.equal(el.style.color, 'blue', 'editing the input restyles');
  instance.dispose();
});

spec('editor round-trip: set + clear leaves the document byte-identical, undo and redo restyle', () => {
  const reg = registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-text-input', name: 'name', props: {label: 'Name'}},
  }, new SpecContext(), reg);
  const before = JSON.stringify(instance.dump());
  const editor = new SpecEditor(instance);
  const node = instance.spec.root;

  editor.apply({op: 'set-prop', node, name: 'color', value: 'red'});
  assert.equal(instance.node('name').root.style.color, 'red');
  const styled = JSON.stringify(instance.dump());
  assert.equal(JSON.parse(styled).root.props.color, 'red');

  // applyAll: a separate undo entry — apply would coalesce the clear into the set
  editor.applyAll([{op: 'set-prop', node, name: 'color', value: undefined}]);
  assert.equal(JSON.stringify(instance.dump()), before, 'the cleared prop leaves no trace');
  assert.equal(instance.node('name').root.style.color, undefined, 'the rebuilt node carries no style');

  editor.undo();
  assert.equal(JSON.stringify(instance.dump()), styled);
  assert.equal(instance.node('name').root.style.color, 'red', 'undo of the clear restyles');
  editor.undo();
  assert.equal(JSON.stringify(instance.dump()), before);
  editor.redo();
  assert.equal(JSON.stringify(instance.dump()), styled);
  assert.equal(instance.node('name').root.style.color, 'red');
  editor.redo();
  assert.equal(JSON.stringify(instance.dump()), before);
  assert.equal(instance.node('name').root.style.color, undefined);
  instance.dispose();
});

spec('canApply: a choice outside the list and a bind on an HTML node are refused', () => {
  const reg = registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [{tag: 'u2-text-input', name: 'name', props: {label: 'Name'}}]},
  }, new SpecContext(), reg);
  const editor = new SpecEditor(instance);
  const input = instance.spec.root.children[0];
  assert.match(editor.canApply({op: 'set-prop', node: input, name: 'cursor', value: 'nope'}),
    /must be one of/);
  assert.equal(editor.canApply({op: 'set-bind', node: instance.spec.root, name: 'color', path: '$.x'}),
    'prop "color" is not bindable');
  instance.dispose();
});

spec('manifest: the shared block once at the top level, filtered out of every component', () => {
  const manifest = JSON.parse(JSON.stringify(registry().manifest()));
  assert.deepEqual(manifest.appearance.map((p) => p.name), APPEARANCE_PROPS.map((p) => p.name));
  assert.equal(manifest.appearance.length, 13);
  for (const meta of manifest.components) {
    assert.equal(meta.props.some((p) => p.category === APPEARANCE_CATEGORY), false,
      `${meta.tag} duplicates the shared block`);
  }

  const reg = new Registry();
  const own = {name: 'width', type: 'string', description: 'Own width'};
  reg.register({tag: 'u2-fake-box', description: 'Fake', props: [own],
    create: () => new Control(), example: {tag: 'u2-fake-box'}});
  reg.register({tag: 'u2-fake-viewer', appearance: false, description: 'Fake', props: [],
    create: () => new Control(), example: {tag: 'u2-fake-viewer'}});
  const metas = new Map(JSON.parse(JSON.stringify(reg.manifest())).components.map((c) => [c.tag, c]));
  assert.deepEqual(metas.get('u2-fake-box').props, [own], 'the colliding own prop survives the filter');
  assert.equal(metas.get('u2-fake-viewer').appearance, false, 'the opt-out is flagged');
});

spec('html node: appearance props style the element with no warning, opacity as a double', () => {
  const reg = registry();
  let instance;
  const warnings = captureWarnings(() => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'div', name: 'card', props: {
        backgroundColor: 'var(--dg-bg-secondary)', padding: '8px 12px', borderRadius: '4px', opacity: 0.4}},
    }, new SpecContext(), reg);
  });
  assert.deepEqual(warnings, []);
  const el = instance.node('card');
  assert.equal(el.style.backgroundColor, 'var(--dg-bg-secondary)');
  assert.equal(el.style.padding, '8px 12px');
  assert.equal(el.style.borderRadius, '4px');
  assert.equal(el.style.opacity, '0.4');
  instance.dispose();

  const props = htmlProps('div');
  const appearance = props.filter((p) => p.category === APPEARANCE_CATEGORY);
  assert.equal(appearance.length, 13);
  for (const p of appearance)
    assert.equal(p.bindable, undefined, `${p.name} must not advertise bindable on an HTML node`);
});
