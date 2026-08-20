/* The real registrations: every v1 component registers with metadata that survives JSON, and every
   registered example renders through the spec renderer without a placeholder and disposes clean. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
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

function components(reg) {
  return reg.manifest().components;
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

spec('registerAll: fills a fresh registry and is idempotent', () => {
  const reg = registry();
  const tags = components(reg).map((c) => c.tag);
  assert.equal(tags.length >= 12, true, `expected at least 12 components, got ${tags.length}`);
  assert.equal(new Set(tags).size, tags.length, 'no tag registered twice');
  for (const tag of tags) {
    assert.match(tag, /^u2-/);
    const meta = reg.get(tag);
    assert.equal((meta.visual === false ? meta.createComponent : meta.create) instanceof Function,
      true, `${tag} carries no factory for what it is`);
  }
  registerAll(reg);
  assert.equal(components(reg).length, tags.length, 'a second call registers nothing new');
});

spec('manifest: survives JSON and describes every component', () => {
  const manifest = registry().manifest();
  const restored = JSON.parse(JSON.stringify(manifest));
  assert.deepEqual(restored, manifest);
  assert.equal(restored.schema, 'dg-ui/1');
  for (const meta of restored.components) {
    assert.equal('create' in meta, false, `${meta.tag} leaks its factory into the manifest`);
    assert.equal(typeof meta.description === 'string' && meta.description.length > 0, true,
      `${meta.tag} has no description`);
    assert.equal(meta.example.tag, meta.tag, `${meta.tag} has no example of itself`);
    assert.equal(Array.isArray(meta.props), true);
    for (const prop of [...meta.props, ...meta.childProps ?? []])
      assert.match(prop.type, /^(string|int|double|bool|string_list|object)$/);
  }
});

spec('every registered example renders without a placeholder and disposes clean', () => {
  const reg = registry();
  const live = Scope.liveCount;
  const instances = [];
  const warnings = captureWarnings(() => {
    for (const meta of components(reg)) {
      // a non-visual example belongs on the tray: rendered as the root it is a placeholder by design
      const spec = meta.visual === false ?
        {$schema: 'dg-ui/1', root: {tag: 'div'}, components: [meta.example]} :
        {$schema: 'dg-ui/1', root: meta.example};
      const instance = renderSpec(spec, new SpecContext(), reg);
      const errors = instance.root.querySelectorAll('.u2-spec-error');
      assert.equal(errors.length, 0, `${meta.tag}: ${errors.map((e) => e.textContent).join('; ')}`);
      instances.push(instance);
    }
  });
  assert.deepEqual(warnings, [], 'every example is fully described by its metadata');

  for (const instance of instances)
    instance.dispose();
  assert.equal(Scope.liveCount, live, 'every example released its scopes');
});

spec('u2-text-input: the bound value drives the context signal and back', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {name: 'Aspirin'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-text-input', props: {label: 'Name'}, bind: {value: '$.name'}},
  }, ctx, reg);
  const editor = instance.root.querySelector('input');

  assert.equal(editor.value, 'Aspirin');
  type(editor, 'Ibuprofen');
  assert.equal(ctx.data.name.value, 'Ibuprofen', 'typing reaches the context signal');
  ctx.data.name.value = 'Paracetamol';
  assert.equal(editor.value, 'Paracetamol', 'and the signal reaches the editor');
  instance.dispose();
});

spec('u2-form: a misspelled child is a placeholder in its own row, its sibling still renders', () => {
  const reg = registry();
  let instance;
  const warnings = captureWarnings(() => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'u2-form', children: [
        {tag: 'u2-txt-input', props: {label: 'Typo'}},
        {tag: 'u2-text-input', props: {label: 'Name', value: 'Aspirin'}},
      ]},
    }, new SpecContext(), reg);
  });

  const rows = instance.root.querySelector('.u2-form-rows');
  assert.equal(rows.children.length, 2, 'the placeholder is a row, not a trailer on the form root');
  assert.equal(rows.children[0].classList.contains('u2-spec-error'), true, 'and it holds its child index');
  assert.match(rows.children[0].textContent, /^u2-txt-input: /);
  assert.equal(instance.root.querySelector('input').value, 'Aspirin');
  assert.equal(warnings.length, 0, 'a child rendered in place warrants no warning');
  instance.dispose();
});

spec('u2-form: DOM order follows spec order across inputs and plain children alike', () => {
  const reg = registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-form', children: [
      {tag: 'u2-text-input', name: 'a', props: {label: 'A'}},
      {tag: 'u2-button', name: 'b', props: {text: 'Go'}},
      {tag: 'u2-text-input', name: 'c', props: {label: 'C'}},
    ]},
  }, new SpecContext(), reg);

  const rows = instance.root.querySelector('.u2-form-rows');
  assert.deepEqual([...rows.children].map((el) => el.dataset.u2Name), ['a', 'b', 'c']);
  instance.dispose();
});

spec('u2-form: a child broken by a patch renders its placeholder in place, and undo heals it there', () => {
  const reg = registry();
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'u2-form', name: 'form', children: [
      {tag: 'u2-text-input', name: 'a', props: {label: 'A'}},
      {tag: 'u2-text-input', name: 'mid', props: {label: 'Mid'}},
      {tag: 'u2-text-input', name: 'z', props: {label: 'Z'}},
    ]},
  };
  const instance = renderSpec(source, new SpecContext(), reg);
  const editor = new SpecEditor(instance);
  const rows = () => instance.root.querySelector('.u2-form-rows');

  editor.apply({op: 'set-bind', node: source.root.children[1], name: 'value', path: '$.nowhere'});
  assert.equal(rows().children[1].classList.contains('u2-spec-error'), true,
    'the broken child keeps its row between its siblings');
  assert.deepEqual([...rows().children].map((el) => el.dataset.u2Name), ['a', 'mid', 'z']);

  editor.undo();
  assert.equal(rows().children[1].classList.contains('u2-spec-error'), false);
  assert.deepEqual([...rows().children].map((el) => el.dataset.u2Name), ['a', 'mid', 'z']);
  instance.dispose();
});

spec('u2-splitter: spec children become the panels, and the sash still resizes', () => {
  const reg = registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-splitter', props: {direction: 'horizontal', sizes: [0.3, 0.7]}, children: [
      {tag: 'u2-panel', children: [{tag: 'u2-text-input', props: {label: 'Name'}}]},
      {tag: 'u2-panel', children: [{tag: 'u2-text-area', props: {label: 'Notes'}}]},
    ]},
  }, new SpecContext(), reg);
  document.body.append(instance.root);

  const panels = instance.root.querySelectorAll('.u2-splitter-panel');
  const sash = instance.root.querySelector('.u2-sash');
  assert.equal(panels.length, 2);
  assert.equal(panels[0].querySelector('input') !== null, true, 'the first child rendered into its panel');
  assert.equal(panels[1].querySelector('textarea') !== null, true);
  assert.equal(sash.getAttribute('aria-valuenow'), '30', 'the sizes prop reached the constructor');

  instance.root.querySelector('.u2-splitter').clientWidth = 600;
  fire(sash, 'keydown', {key: 'ArrowRight'});
  assert.equal(sash.getAttribute('aria-valuenow'), '32');
  instance.dispose();
});

spec('u2-form: input children join the form and share its label column', async () => {
  const reg = registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-form', children: [
      {tag: 'u2-text-input', props: {label: 'Name', value: 'Aspirin'}},
      {tag: 'u2-number-input', props: {label: 'Amount', value: 10}},
    ]},
  }, new SpecContext(), reg);
  document.body.append(instance.root);

  const rows = instance.root.querySelectorAll('.u2-form-rows .u2-input-root');
  assert.equal(rows.length, 2, 'both inputs went through Form.add');

  const labels = instance.root.querySelectorAll('.u2-input-label');
  labels[0].offsetWidth = 40;
  labels[1].offsetWidth = 64;
  await flush();
  const form = instance.root.querySelector('.u2-form');
  assert.equal(form.style.getPropertyValue('--u2-form-label-width'), '64px', 'the labels are aligned');

  type(instance.root.querySelectorAll('input')[1], 'nonsense');
  assert.equal(instance.root.querySelectorAll('.u2-invalid').length, 1, 'validity reaches the row');
  instance.dispose();
});

spec('u2-form: adopt takes inputs into the form and aggregates their validity', () => {
  const reg = registry();
  const meta = reg.get('u2-form');
  const form = meta.create({});
  const name = new TextInput({label: 'Name'});
  name.addValidator((value) => value ? null : 'Name is required');
  const stray = document.createElement('div');

  const warnings = captureWarnings(() => {
    meta.adopt(form, name, 0);
    meta.adopt(form, stray, 1);
  });
  assert.deepEqual([...form.inputs], [name]);
  assert.equal(form.validity.value, 'Name is required', 'an invalid child invalidates the form');
  name.value.value = 'Aspirin';
  assert.equal(form.validity.value, null);
  assert.deepEqual([...form.root.querySelector('.u2-form-rows').children], [name.root, stray],
    'a non-input is a row in adopt order');
  assert.equal(warnings.length, 0);
  form.dispose();
  name.dispose();
});

spec('u2-accordion: child nodes title the panes, and the headers still toggle', () => {
  const reg = registry();
  let instance;
  const warnings = captureWarnings(() => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'u2-accordion', children: [
        {tag: 'u2-panel', props: {title: 'General'}, children: [{tag: 'u2-text-input', props: {label: 'Name'}}]},
        {tag: 'u2-panel', children: [{tag: 'u2-bool-input', props: {label: 'Active'}}]},
      ]},
    }, new SpecContext(), reg);
  });
  document.body.append(instance.root);

  assert.deepEqual(warnings, [], 'a child title is declared metadata, not an unknown prop');
  const titles = instance.root.querySelectorAll('.u2-accordion-title').map((el) => el.textContent);
  assert.deepEqual(titles, ['General', 'Pane 2']);

  const contents = instance.root.querySelectorAll('.u2-accordion-content');
  assert.equal(contents[0].style.display, 'none');
  fire(instance.root.querySelector('.u2-accordion-header'), 'click');
  assert.equal(contents[0].style.display, '');
  assert.equal(contents[0].querySelector('input') !== null, true, 'the pane holds the rendered child');
  instance.dispose();
});

spec('u2-tabs: child nodes label the tabs, and activation switches panels', () => {
  const reg = registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-tabs', children: [
      {tag: 'u2-panel', props: {title: 'Data'}, children: [{tag: 'u2-text-input', props: {label: 'Table'}}]},
      {tag: 'u2-panel', children: [{tag: 'u2-bool-input', props: {label: 'Legend'}}]},
    ]},
  }, new SpecContext(), reg);
  document.body.append(instance.root);

  const labels = instance.root.querySelectorAll('.u2-tabs-label').map((el) => el.textContent);
  assert.deepEqual(labels, ['Data', 'Tab 2']);
  const panels = instance.root.querySelectorAll('.u2-tabs-panel');
  assert.equal(panels[0].style.display, '');
  assert.equal(panels[1].style.display, 'none');

  fire(instance.root.querySelectorAll('.u2-tabs-tab')[1], 'pointerdown');
  assert.equal(panels[0].style.display, 'none');
  assert.equal(panels[1].style.display, '');
  assert.equal(panels[1].querySelector('input') !== null, true, 'the panel holds the rendered child');
  instance.dispose();
});

spec('u2-tabs and u2-accordion: the child icon prop renders, the tabs props pick orientation and variant', () => {
  const reg = registry();
  let instance;
  const warnings = captureWarnings(() => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'u2-div-v', children: [
        {tag: 'u2-tabs', props: {orientation: 'vertical', variant: 'document'}, children: [
          {tag: 'u2-panel', props: {title: 'Data', icon: 'table'}},
          {tag: 'u2-panel', props: {title: 'Style'}},
        ]},
        {tag: 'u2-accordion', children: [
          {tag: 'u2-panel', props: {title: 'General', icon: 'cog'}},
          {tag: 'u2-panel', props: {title: 'Notes'}},
        ]},
      ]},
    }, new SpecContext(), reg);
  });
  document.body.append(instance.root);

  assert.deepEqual(warnings, [], 'icon is declared child metadata on both hosts');
  const tabs = instance.root.querySelector('.u2-tabs');
  assert.equal(tabs.classList.contains('u2-tabs-vertical'), true);
  assert.equal(tabs.classList.contains('u2-tabs-document'), true);
  const tabIcons = tabs.querySelectorAll('.u2-tabs-tab').map((t) => t.querySelector('.u2-tabs-icon .fa-table'));
  assert.equal(tabIcons[0] !== null, true, 'the first tab carries its icon');
  assert.equal(tabIcons[1], null);
  const paneIcons = instance.root.querySelectorAll('.u2-accordion-header')
    .map((h) => h.querySelector('.u2-accordion-icon .fa-cog'));
  assert.equal(paneIcons[0] !== null, true, 'the first pane carries its icon');
  assert.equal(paneIcons[1], null);
  instance.dispose();
});

spec('u2-property-grid: its JSON payload renders without a single warning', () => {
  const reg = registry();
  let instance;
  const warnings = captureWarnings(() => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'u2-property-grid', props: {
        properties: [{name: 'Title', type: 'string'}, {name: 'Rows', type: 'int', min: 0}],
        values: {Title: 'Demo', Rows: 20},
      }},
    }, new SpecContext(), reg);
  });

  assert.deepEqual(warnings, []);
  assert.equal(instance.root.querySelectorAll('.u2-propgrid-row').length, 2);
  assert.equal(instance.root.querySelector('input').value, 'Demo');
  instance.dispose();
});

spec('manifest: the child hooks stay out of it, the metadata they need stays in', () => {
  const reg = registry();
  assert.equal(typeof reg.get('u2-splitter').createWithChildren, 'function');
  assert.equal(typeof reg.get('u2-form').adopt, 'function');

  const manifest = JSON.parse(JSON.stringify(reg.manifest()));
  const metas = new Map(manifest.components.map((c) => [c.tag, c]));
  for (const tag of ['u2-splitter', 'u2-accordion', 'u2-tabs', 'u2-form']) {
    const meta = metas.get(tag);
    assert.equal('createWithChildren' in meta, false, `${tag} leaks a hook into the manifest`);
    assert.equal('adopt' in meta, false, `${tag} leaks a hook into the manifest`);
    assert.equal(meta.acceptsChildren, true, `${tag} does not advertise its children`);
  }
  assert.deepEqual(metas.get('u2-tabs').childProps.map((p) => p.name), ['title', 'icon']);
  assert.deepEqual(metas.get('u2-accordion').childProps.map((p) => p.name), ['title', 'icon']);
  const tabProps = Object.fromEntries(metas.get('u2-tabs').props.map((p) => [p.name, p.choices]));
  assert.deepEqual(tabProps, {orientation: ['horizontal', 'vertical'], variant: ['platform', 'document']});
  assert.deepEqual(metas.get('u2-tabs').defaults, {orientation: 'horizontal', variant: 'platform'});
  assert.deepEqual(metas.get('u2-property-grid').props.map((p) => p.type), ['object', 'object']);
});

spec('dump: a splitter of a form and an accordion round-trips', () => {
  const reg = registry();
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'u2-splitter', props: {direction: 'vertical'}, children: [
      {tag: 'u2-form', children: [{tag: 'u2-text-input', props: {label: 'Name', value: 'Aspirin'}}]},
      {tag: 'u2-accordion', children: [
        {tag: 'u2-panel', props: {title: 'Notes'}, children: [{tag: 'u2-text-area', props: {label: 'Notes'}}]},
      ]},
    ]},
  };
  const instance = renderSpec(source, new SpecContext(), reg);
  assert.deepEqual(instance.dump(), source, 'an untouched tree dumps back to what it was');

  type(instance.root.querySelector('input'), 'Ibuprofen');
  const dumped = instance.dump();
  assert.equal(dumped.root.children[0].children[0].props.value, 'Ibuprofen');

  const restored = renderSpec(dumped, new SpecContext(), reg);
  assert.equal(restored.root.querySelector('input').value, 'Ibuprofen');
  assert.deepEqual(restored.dump(), dumped);
  instance.dispose();
  restored.dispose();
});
