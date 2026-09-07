/* The control inspector (src/dg/shell/control-inspector.ts): a u2 control made
   `grok.shell.o` renders as a two-way PropertyGrid over its registry props, read-only where a
   property has no setter, and each render replaces the previous panel's scope. `DG.ObjectHandler`
   comes from tests/dg-stub.mjs. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {TextInput} from '../src/components/inputs/text-input.js';

register('./dg-stub.mjs', import.meta.url);
const {controlPropDescriptors, controlProperties, registerControlInspector, disposePanel,
  controlAt, noControl} = await import('../src/dg/shell/control-inspector.js');
const DG = await import('datagrok-api/dg');

registerControlInspector();
registerControlInspector();

function inspector(name, body) {
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

inspector('a TextInput answers editable descriptors; non-editable types are excluded', () => {
  const input = new TextInput({label: 'Name', value: 'Aspirin'});
  const names = controlPropDescriptors(input).map((p) => p.name);
  assert.ok(names.includes('label') && names.includes('value'), `got ${names.join(', ')}`);

  const bare = new Control();
  bare.title = 'Fake';
  bare.componentMeta = {tag: 'u2-fake', props: [
    {name: 'title', type: 'string'},
    {name: 'sizes', type: 'object'},
  ]};
  assert.deepEqual(controlPropDescriptors(bare).map((p) => p.name), ['title'],
    'the object-typed prop is not a grid row');
  input.dispose();
  bare.dispose();
});

inspector('a declared prop whose value cannot be read is not a row at all', () => {
  // the registry describes every TextInput option, but this one was built by hand and given none
  // of them: `password` and `autoResize` read nothing, so a row for either could only be blank
  const plain = new TextInput({label: 'Name'});
  const names = controlPropDescriptors(plain).map((p) => p.name);
  assert.ok(!names.includes('password') && !names.includes('autoResize'),
    `unreadable props still listed: ${names.join(', ')}`);

  const given = new TextInput({label: 'Name', password: true});
  assert.ok(controlPropDescriptors(given).map((p) => p.name).includes('password'),
    'a prop the constructor was given reads, so it stays a row');
  plain.dispose();
  given.dispose();
});

inspector('a prop without a setter maps to readonly', () => {
  const input = new TextInput({search: true});
  const by = Object.fromEntries(controlPropDescriptors(input).map((p) => [p.name, p]));
  assert.equal(by.search.readonly, true, 'an options-only prop has no setter');
  assert.equal(by.label.readonly, false, 'label writes through the accessor');
  assert.equal(by.value.readonly, false, 'value writes through the signal');
  input.dispose();
});

inspector('a grid edit writes through to the control', async () => {
  const input = new TextInput({label: 'Name', value: 'Aspirin'});
  document.body.append(input.root);
  const grid = controlProperties(input);
  document.body.append(grid.root);
  const field = (name) => [...grid.root.querySelectorAll('.u2-propgrid-row')]
    .find((r) => r.querySelector('.u2-propgrid-name')?.textContent === name)
    ?.querySelector('input');

  const value = field('value');
  assert.equal(value.value, 'Aspirin', 'the grid reads the live value');
  value.value = 'Ibuprofen';
  fire(value, 'input');
  await flush();
  assert.equal(input.value.peek(), 'Ibuprofen', 'the value edit reached the signal');

  const label = field('label');
  label.value = 'Renamed';
  fire(label, 'input');
  await flush();
  assert.equal(input.label, 'Renamed', 'the label edit reached the rendered label');
  grid.dispose();
  input.dispose();
});

inspector('renderProperties twice disposes the first panel scope', () => {
  assert.equal(DG.ObjectHandler.registered.filter((h) => h.type === 'u2-control').length, 1,
    'a second registerControlInspector() call stacks no duplicate');
  const input = new TextInput({label: 'Name', value: 'Aspirin'});
  const handler = DG.ObjectHandler.forEntity(input);
  assert.equal(handler?.type, 'u2-control', 'the handler claims a u2 control');

  const first = handler.renderProperties(input);
  const live = Scope.liveCount;
  const second = handler.renderProperties(input);
  assert.equal(Scope.liveCount, live, 'the second render replaced the first panel scope');
  assert.notEqual(first, second);
  assert.ok(second.querySelector('[data-u2="property-grid"]'), 'the panel holds the grid');
  assert.equal(second.querySelector('h3')?.textContent, 'u2-text-input',
    'the header names the control tag');
  disposePanel();
  input.dispose();
});

inspector('controlAt resolves only inside the given subtree', () => {
  const outer = new Control();
  const host = document.createElement('div');
  outer.root.append(host);
  const inner = new TextInput({label: 'Name'});
  const plain = document.createElement('div');
  host.append(inner.root, plain);

  assert.equal(controlAt(inner.root.querySelector('input') ?? inner.root, host), inner,
    'the nearest control inside the subtree');
  assert.equal(controlAt(plain, host), undefined,
    'an element under no control does not fall through to the ancestor control');
  assert.equal(controlAt(host, host), undefined, 'the boundary itself is never the answer');
  assert.equal(controlAt(outer.root, host), undefined, 'nothing outside the subtree');
  inner.dispose();
  outer.dispose();
});

inspector('a control the registry does not describe says so instead of rendering an empty grid', () => {
  const bare = new Control();
  bare.root.dataset.u2 = 'list';
  const handler = DG.ObjectHandler.forEntity(bare);
  const panel = handler.renderProperties(bare);
  assert.match(panel.textContent, /No inspectable properties/);
  assert.equal(panel.querySelector('[data-u2="property-grid"]'), null);

  const empty = DG.ObjectHandler.forEntity(noControl);
  assert.equal(empty?.type, 'u2-no-control', 'the empty state has a handler of its own');
  assert.match(empty.renderProperties(noControl).textContent, /not a u2 component/);
  disposePanel();
  bare.dispose();
});
