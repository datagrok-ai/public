/* Where a drag may land (WO-9): the pure half of drag-and-drop. Legality comes from registry child
   metadata, geometry from the rendered elements — and the shim lays nothing out, so every test
   feeds the rects itself (`el.rect` is what getBoundingClientRect answers). */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {accepts, resolveDrop} from '../src/dg/designer/dnd.js';

function dnd(name, body) {
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
  reg.register({
    tag: 'u2-fake-box', description: 'Fake container', props: [], acceptsChildren: true,
    create: () => new Control(), example: {tag: 'u2-fake-box'},
  });
  reg.register({
    tag: 'u2-fake-form', description: 'Fake adopter', props: [],
    create: () => new Control(),
    adopt: (parent, child) => parent.root.append(child instanceof Control ? child.root : child),
    example: {tag: 'u2-fake-form'},
  });
  reg.register({
    tag: 'u2-fake-tabs', description: 'Fake constructor-children container', props: [],
    create: () => new Control(),
    createWithChildren: (_props, children) => {
      const control = new Control();
      for (const child of children)
        control.root.append(child instanceof Control ? child.root : child);
      return control;
    },
    example: {tag: 'u2-fake-tabs'},
  });
  reg.register({
    tag: 'u2-fake-leaf', description: 'Fake leaf', props: [],
    create: () => new Control(), example: {tag: 'u2-fake-leaf'},
  });
  return reg;
}

function render(root, reg) {
  const warn = console.warn;
  console.warn = () => {};
  try {
    return renderSpec({$schema: 'dg-ui/1', root}, new SpecContext(), reg);
  } finally {
    console.warn = warn;
  }
}

/** The shim has no layout: this is where a node "is" for the resolver. */
function place(instance, node, x, y, width, height) {
  const built = instance.nodes().get(node);
  const el = built instanceof Control ? built.root : built;
  el.rect = new DOMRect(x, y, width, height);
  return el;
}

function stack(instance, parent, direction) {
  const built = instance.nodes().get(parent);
  (built instanceof Control ? built.root : built).style.flexDirection = direction;
}

dnd('accepts(): a container, an adopter and a constructor-children parent take children; a leaf does not', () => {
  const reg = registry();
  assert.equal(accepts(reg, {tag: 'u2-fake-box'}), true);
  assert.equal(accepts(reg, {tag: 'u2-fake-form'}), true, 'adopt implies children');
  assert.equal(accepts(reg, {tag: 'u2-fake-tabs'}), true, 'createWithChildren implies children');
  assert.equal(accepts(reg, {tag: 'div'}), true, 'a supported HTML tag is a container');
  assert.equal(accepts(reg, {tag: 'u2-fake-leaf'}), false);
  assert.equal(accepts(reg, {tag: 'u2-nope'}), false, 'an unknown tag takes nothing');
});

dnd('resolveDrop(): a container takes the drop into itself, at the end of what it holds', () => {
  const reg = registry();
  const empty = {tag: 'u2-fake-box', name: 'empty'};
  const one = {tag: 'u2-fake-leaf', name: 'one'};
  const box = {tag: 'u2-fake-box', name: 'box', children: [one]};
  const root = {tag: 'u2-fake-box', name: 'root', children: [empty, box]};
  const instance = render(root, reg);
  place(instance, empty, 0, 0, 200, 40);
  place(instance, box, 0, 40, 200, 40);

  const into = resolveDrop(instance, empty, 100, 20);
  assert.deepEqual(into, {parent: empty, index: 0, kind: 'into', rect: {x: 0, y: 0, width: 200, height: 40}});
  assert.equal(resolveDrop(instance, box, 100, 60).index, 1, 'after what the container already holds');
  instance.dispose();
});

dnd('resolveDrop(): a leaf splits its parent at the midpoint, on the parent\'s own flow axis', () => {
  const reg = registry();
  const first = {tag: 'u2-fake-leaf', name: 'first'};
  const second = {tag: 'u2-fake-leaf', name: 'second'};
  const column = {tag: 'u2-fake-box', name: 'column', children: [first, second]};
  const instance = render(column, reg);
  place(instance, column, 0, 0, 200, 80);
  place(instance, first, 0, 0, 200, 40);
  place(instance, second, 0, 40, 200, 40);

  stack(instance, column, 'column');
  const before = resolveDrop(instance, second, 100, 45);
  assert.equal(before.kind, 'line');
  assert.equal(before.index, 1, 'above the midpoint drops before the leaf');
  assert.deepEqual(before.rect, {x: 0, y: 40, width: 200, height: 2});
  const after = resolveDrop(instance, second, 100, 75);
  assert.equal(after.index, 2, 'below the midpoint drops after it');
  assert.deepEqual(after.rect, {x: 0, y: 80, width: 200, height: 2});

  stack(instance, column, 'row');
  assert.equal(resolveDrop(instance, first, 50, 20).index, 0, 'a row splits on x, not y');
  const right = resolveDrop(instance, first, 150, 20);
  assert.equal(right.index, 1);
  assert.deepEqual(right.rect, {x: 200, y: 0, width: 2, height: 40});
  instance.dispose();
});

dnd('resolveDrop(): nothing to drop into — a leaf root, and a leaf under a parent that takes none', () => {
  const reg = registry();
  const leaf = {tag: 'u2-fake-leaf', name: 'leaf'};
  const instance = render(leaf, reg);
  place(instance, leaf, 0, 0, 200, 40);
  assert.equal(resolveDrop(instance, leaf, 100, 20), null, 'the root has no parent to drop beside');
  instance.dispose();

  const child = {tag: 'u2-fake-leaf', name: 'child'};
  const dropped = render({tag: 'u2-fake-leaf', name: 'outer', children: [child]}, reg);
  assert.equal(dropped.nodes().has(child), false, 'a non-accepting parent never rendered it');
  dropped.dispose();
});

dnd('resolveDrop(): a move refuses itself, its own subtree and every position it is already in', () => {
  const reg = registry();
  const inner = {tag: 'u2-fake-leaf', name: 'inner'};
  const group = {tag: 'u2-fake-box', name: 'group', children: [inner]};
  const first = {tag: 'u2-fake-leaf', name: 'first'};
  const last = {tag: 'u2-fake-leaf', name: 'last'};
  const root = {tag: 'u2-fake-box', name: 'root', children: [first, group, last]};
  const instance = render(root, reg);
  place(instance, root, 0, 0, 200, 120);
  place(instance, first, 0, 0, 200, 40);
  place(instance, group, 0, 40, 200, 40);
  place(instance, inner, 0, 40, 200, 40);
  place(instance, last, 0, 80, 200, 40);
  stack(instance, root, 'column');

  assert.equal(resolveDrop(instance, group, 100, 60, group), null, 'into itself');
  assert.equal(resolveDrop(instance, inner, 100, 50, group), null, 'into its own subtree');
  assert.equal(resolveDrop(instance, first, 100, 5, first), null, 'before itself is where it is');
  assert.equal(resolveDrop(instance, first, 100, 35, first), null, 'and so is right after itself');
  assert.equal(resolveDrop(instance, last, 100, 115, first).index, 2,
    'the index a move needs counts the node as already spliced out');
  assert.equal(resolveDrop(instance, group, 100, 60, first).index, 1, 'into a sibling container');
  assert.equal(resolveDrop(instance, root, 100, 10, last), null, 'into the parent it already ends');
  assert.equal(resolveDrop(instance, root, 100, 10, first).index, 2, 'the tail of the same parent');
  instance.dispose();
});

dnd('resolveDrop(): a null hit — canvas padding — appends into the root when it takes children', () => {
  const reg = registry();
  const child = {tag: 'u2-fake-leaf', name: 'only'};
  const root = {tag: 'u2-fake-box', name: 'root', children: [child]};
  const instance = render(root, reg);
  place(instance, root, 0, 0, 400, 300);

  const target = resolveDrop(instance, null, 500, 500);
  assert.equal(target.parent, root, 'no rendered node under the pointer still means the root');
  assert.equal(target.index, 1, 'appended at the end');
  assert.equal(target.kind, 'into');
  assert.equal(resolveDrop(instance, null, 500, 500, child), null,
    'moving the last child to the end is nothing to do');
  instance.dispose();

  const leaf = render({tag: 'u2-fake-leaf', name: 'solo'}, reg);
  assert.equal(resolveDrop(leaf, null, 0, 0), null, 'a root that takes no children is no fallback');
  leaf.dispose();
});
