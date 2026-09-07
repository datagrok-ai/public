/* The live property surface: meta-declared props backed by accessors are writable through
   getProperties() (the context-panel/automation door), signal-backed ones stay writable as
   before, and construction-only options stay read-only. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope, registry, registerAll, TextInput, ChoiceInput, TabStrip, RangeSlider}
  from '../src/index.js';

registerAll();

function smoke(name, body) {
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

function stamped(control, tag) {
  control.componentMeta = registry.get(tag);
  document.body.append(control.root);
  return control;
}

function prop(control, name) {
  return control.getProperties().find((p) => p.name === name);
}

smoke('accessor-backed props are writable; construction-only options stay read-only', async () => {
  const input = stamped(new TextInput({label: 'Name', placeholder: 'type here'}), 'u2-text-input');

  const label = prop(input, 'label');
  assert.equal(label.get(input), 'Name');
  assert.ok(label.set, 'the label accessor makes it writable');
  label.set(input, 'Compound');
  assert.equal(input.root.querySelector('.u2-input-label').textContent, 'Compound');

  const placeholder = prop(input, 'placeholder');
  placeholder.set(input, 'e.g. aspirin');
  assert.equal(input.root.querySelector('input').placeholder, 'e.g. aspirin');

  const enabled = prop(input, 'enabled');
  assert.equal(enabled.get(input), true);
  enabled.set(input, false);
  assert.ok(input.root.classList.contains('u2-input-disabled'));
  assert.ok(input.root.querySelector('input').disabled);

  assert.equal(prop(input, 'inline').set, undefined, 'a plain option has no live setter');
  prop(input, 'value').set(input, 'aspirin');
  assert.equal(input.value.value, 'aspirin', 'signal-backed props write as before');
  input.dispose();
});

smoke('a label set on an input created without one materializes the element', async () => {
  const input = stamped(new TextInput({}), 'u2-text-input');
  assert.equal(input.root.querySelector('.u2-input-label'), null);
  prop(input, 'label').set(input, 'Late label');
  assert.equal(input.root.querySelector('.u2-input-label').textContent, 'Late label');
  assert.equal(input.label, 'Late label');

  const inline = stamped(new TextInput({inline: true}), 'u2-text-input');
  prop(inline, 'label').set(inline, 'Never');
  assert.equal(inline.root.querySelector('.u2-input-label'), null, 'inline inputs take no label');
  input.dispose();
  inline.dispose();
});

smoke('enabled arrives as a construction option too', async () => {
  const input = new TextInput({label: 'Off', enabled: false});
  document.body.append(input.root);
  assert.ok(input.root.querySelector('input').disabled);
  input.dispose();
});

smoke('choice items rewrite the select and prune a vanished value as a system write', async () => {
  const input = stamped(new ChoiceInput({label: 'Series', items: ['a', 'b'], value: 'b'}),
    'u2-choice-input');
  prop(input, 'items').set(input, ['a', 'c']);
  const values = [...input.root.querySelectorAll('option')].map((o) => o.value);
  assert.deepEqual(values, ['', 'a', 'c']);
  assert.equal(input.value.value, null, 'the vanished value was pruned');
  assert.deepEqual(prop(input, 'items').get(input), ['a', 'c']);
  input.dispose();
});

smoke('activeTab is a declared, two-way property of the tab strip', async () => {
  const tabs = stamped(new TabStrip({tabs: [
    {id: 'data', label: 'Data', content: document.createElement('div')},
    {id: 'style', label: 'Style', content: document.createElement('div')},
  ]}), 'u2-tabs');
  await flush();

  const active = prop(tabs, 'activeTab');
  assert.ok(active, 'declared in the registry meta');
  assert.equal(active.get(tabs), 'data');
  active.set(tabs, 'style');
  await flush();
  const panels = tabs.root.querySelectorAll('.u2-tabs-panel');
  assert.equal(panels[0].style.display, 'none');
  assert.equal(panels[1].style.display, '');
  assert.equal(tabs.activeTab.value, 'style', 'the write went through the signal itself');
  tabs.dispose();
});

smoke('range-slider bounds are live: a min write re-renders the handles', async () => {
  const slider = stamped(new RangeSlider({min: 0, max: 100, lo: 50, hi: 100}), 'u2-range-slider');
  await flush();
  assert.equal(slider.root.querySelector('.u2-range-handle-lo').style.left, '50%');

  prop(slider, 'min').set(slider, 50);
  await flush();
  assert.equal(slider.root.querySelector('.u2-range-handle-lo').style.left, '0%',
    'the same lo value is now the left edge');
  assert.equal(prop(slider, 'min').get(slider), 50);
  assert.equal(prop(slider, 'vertical').set, undefined, 'orientation stays construction-only');
  slider.dispose();
});
