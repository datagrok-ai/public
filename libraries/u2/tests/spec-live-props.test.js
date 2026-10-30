/* Every registered prop is live (R7, WO-C15): a bind follows its source in place — bindable props as
   signals, every other prop through the component's own setter — and only props marked structural
   read a bind once. Covers the registry walk that keeps registrations honest, in-place updates
   inside an adopting form (P3), the promotion loop that can no longer happen (P4), malformed named
   nodes (P27) and the spec side of the wizard (P93). */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {signal} from '../src/core/signals.js';
import {Scope} from '../src/core/scope.js';
import {Component} from '../src/core/component.js';
import {Registry} from '../src/spec/registry.js';
import {APPEARANCE_PROPS} from '../src/spec/appearance.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {SpecEditor} from '../src/spec/editor.js';
import {registerAll} from '../src/spec/registrations.js';

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

/** Counts every registry build, so "nothing was rebuilt" is a number. */
function countCreates(reg) {
  const counter = {n: 0};
  for (const meta of reg.metas()) {
    for (const hook of ['create', 'createWithChildren', 'createComponent']) {
      const fn = meta[hook];
      if (typeof fn === 'function')
        meta[hook] = (...args) => (counter.n++, fn(...args));
    }
  }
  return counter;
}

async function quietly(body) {
  const warnings = [];
  const original = console.warn;
  console.warn = (message) => warnings.push(String(message));
  try {
    await body();
  } finally {
    console.warn = original;
  }
  return warnings;
}

spec('registry walk: every declared prop of every visual tag is bindable, settable or structural', async () => {
  const reg = registry();
  const missing = [];
  await quietly(() => {
    for (const meta of reg.metas()) {
      if (meta.visual === false)
        continue;
      const instance = renderSpec({$schema: 'dg-ui/1', root: meta.example}, new SpecContext(), reg);
      const built = instance.nodes().get(instance.spec.root);
      assert.ok(Component.is(built), `${meta.tag}: the example builds`);
      const props = built.getProperties();
      for (const p of meta.props) {
        if (APPEARANCE_PROPS.includes(p) || p.bindable || p.structural)
          continue;
        if (props.find((x) => x.name === p.name)?.set == null)
          missing.push(`${meta.tag}.${p.name}`);
      }
      instance.dispose();
    }
  });
  assert.deepEqual(missing, [], 'add a setter (an accessor or a same-named signal) or mark the prop structural');
});

spec('every visual tag declares visible, outside the appearance group; the tray does not', () => {
  const reg = registry();
  for (const meta of reg.metas()) {
    const visible = meta.props.find((p) => p.name === 'visible');
    if (meta.visual === false)
      assert.equal(visible, undefined, meta.tag);
    else {
      assert.equal(visible?.type, 'bool', meta.tag);
      assert.equal(APPEARANCE_PROPS.includes(visible), false, meta.tag);
    }
  }
  const manifest = reg.manifest();
  assert.match(manifest.binding, /follows it live, except props marked structural, which read a bind once/);
  const button = manifest.components.find((c) => c.tag === 'u2-button');
  assert.ok(button.props.some((p) => p.name === 'visible'), 'the manifest lists it');
  const wizard = manifest.components.find((c) => c.tag === 'u2-wizard');
  assert.equal('setChildProp' in wizard, false, 'closures stay out of the manifest');
});

spec('P3: button text, a tab title, a section title and a number max update in place inside a form',
  async () => {
    const reg = registry();
    const creates = countCreates(reg);
    const ctx = new SpecContext({data: {s: 'Save', t: 'General', h: 'Study', m: 10}});
    const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-tabs', name: 'tabs', children: [
      {tag: 'u2-panel', name: 'pane', bind: {title: '$.t'}, children: [
        {tag: 'u2-form', name: 'form', children: [
          {tag: 'u2-text-input', name: 'first', props: {label: 'First'}},
          {tag: 'u2-number-input', name: 'amount', props: {label: 'Amount', value: 5, min: 0}, bind: {max: '$.m'}},
          {tag: 'u2-button', name: 'go', bind: {text: '$.s'}},
          {tag: 'u2-section', name: 'section', bind: {title: '$.h'},
            children: [{tag: 'u2-text-input', name: 'inner', props: {label: 'Inner'}}]},
        ]},
      ]},
    ]}}, ctx, reg);
    document.body.append(instance.root);
    const form = instance.node('form');
    const children = [...form.root.querySelectorAll('*')];
    const editor = instance.node('first').root.querySelector('input');
    editor.focus();
    const built = creates.n;

    instance.node('amount').value = 15;
    assert.match(instance.node('amount').validity.peek(), /at most 10/);
    ctx.data.s.value = 'Save all';
    ctx.data.t.value = 'Details';
    ctx.data.h.value = 'Trial';
    ctx.data.m.value = 20;
    await flush();

    assert.equal(creates.n, built, 'nothing was rebuilt');
    assert.equal(instance.node('form'), form);
    assert.deepEqual([...form.root.querySelectorAll('*')], children, 'every element under the form is the same');
    assert.equal(document.activeElement, editor, 'the focused input keeps focus');
    assert.equal(instance.node('go').root.textContent, 'Save all');
    assert.equal(instance.root.querySelector('.u2-tabs-label').textContent, 'Details');
    assert.equal(instance.node('section').root.textContent.includes('Trial'), true);
    assert.equal(instance.node('amount').max, 20);
    assert.equal(instance.node('amount').validity.peek(), null, '15 is within the new max');
    ctx.data.m.value = 12;
    assert.match(instance.node('amount').validity.peek(), /at most 12/, 'the range check follows the bound max');
    instance.dispose();
  });

spec('P4: form + filter input + stat card + tray source settle in two flushes with no rebuild', async () => {
  const reg = registry();
  let queries = 0;
  // the stand-in for a domain source: `query` re-runs, `state` reports
  reg.register({
    tag: 'u2-fake-source',
    visual: false,
    description: 'Fake source for the P4 reproduction',
    props: [{name: 'query', type: 'string', bindable: true}],
    createComponent: (_props, env) => {
      const source = new Component();
      source.state = signal('idle');
      source.bindStep = (name) => name === 'state' ? source.state : null;
      source.bindProps = () => [{name: 'state', type: 'string'}];
      source.start = () => {
        const query = env.resolve(env.subBinds.query);
        source.scope.effect(() => {
          queries++;
          source.state.value = `loaded ${query?.value ?? ''}`;
        });
      };
      return source;
    },
    example: {tag: 'u2-fake-source', name: 'issues'},
  });
  const creates = countCreates(reg);
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    components: [{tag: 'u2-fake-source', name: 'issues', bind: {query: '$.q.value'}}],
    root: {tag: 'u2-form', name: 'form', children: [
      {tag: 'u2-text-input', name: 'q', props: {label: 'Search', search: true}},
      {tag: 'u2-stat-card', name: 'card', props: {value: '0'}, bind: {label: '$.issues.state'}},
    ]},
  }, new SpecContext(), reg);
  document.body.append(instance.root);
  const built = creates.n;
  const label = () => instance.node('card').root.querySelector('[data-u2-part="label"]').textContent;
  assert.equal(label(), 'loaded ');

  instance.node('q').value = 'open';
  await flush();
  await flush();
  assert.equal(label(), 'loaded open', 'the card follows the source');
  assert.equal(creates.n, built, 'no rebuild, no re-mount');
  assert.equal(queries, 2, 'one query per change — no loop');
  instance.dispose();
});

spec('P27: a malformed named node becomes a placeholder with a diagnostic, never a constructor throw', () => {
  const reg = registry();
  let instance;
  assert.doesNotThrow(() => {
    instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-panel', name: 'root', children: [
      {tag: 'u2-text-input', name: 'badBind', bind: {label: 42}},
      {tag: 'u2-panel', name: 'badChildren', children: {tag: 'u2-text-input'}},
      {tag: 'u2-panel', name: 'nullChild', children: [null]},
      {tag: 'u2-text-input', name: 'fine', props: {label: 'Fine'}},
    ]}}, new SpecContext(), reg);
  });
  const errors = [...instance.root.querySelectorAll('.u2-spec-error')].map((e) => e.textContent);
  assert.deepEqual(errors, [
    'u2-text-input: bind "label" must be a "$." path string',
    'u2-panel: "children" must be an array of nodes',
    'u2-panel: children[0] is not a node with a "tag"',
  ]);
  assert.ok(instance.node('fine'), 'siblings render');
  assert.deepEqual(instance.dump().root.children.map((c) => c?.children ?? null),
    [null, {tag: 'u2-text-input'}, [null], null], 'dump keeps the mistakes as written');
  instance.dispose();
});

spec('P93: the wizard declares its events and currentStep; a spec reacts to step and finish', async () => {
  const reg = registry();
  const meta = reg.get('u2-wizard');
  assert.deepEqual(meta.events, ['finish', 'cancel', 'step']);
  assert.equal(meta.props.find((p) => p.name === 'currentStep').bindable, true);
  const fired = [];
  const ctx = new SpecContext({commands: {finished: () => fired.push('finish'), moved: () => fired.push('step')}});
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-panel', children: [
    {tag: 'h3', name: 'where', bind: {text: '$.wiz.currentStep'}},
    {tag: 'u2-wizard', name: 'wiz', on: {finish: 'cmd:finished', step: 'cmd:moved'}, children: [
      {tag: 'u2-panel', props: {title: 'Source', id: 'source'}},
      {tag: 'u2-panel', props: {title: 'Options', id: 'options'}},
    ]},
  ]}}, ctx, reg);
  document.body.append(instance.root);
  const wiz = instance.node('wiz');
  assert.equal(instance.root.querySelector('h3').textContent, 'source');
  wiz.next();
  assert.equal(instance.root.querySelector('h3').textContent, 'options', 'currentStep is a bind source');
  wiz.next();
  assert.deepEqual(fired, ['step', 'finish']);
  instance.dispose();
});

spec('P93: a bound step title follows in place; a designer patch keeps the step and the visited set',
  async () => {
    const reg = registry();
    const ctx = new SpecContext({data: {t: 'Source'}});
    const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-wizard', name: 'wiz', children: [
      {tag: 'u2-panel', bind: {title: '$.t'}, children: [{tag: 'u2-text-input', name: 'file', props: {label: 'File'}}]},
      {tag: 'u2-panel', props: {title: 'Options'}},
      {tag: 'u2-panel', props: {title: 'Review'}},
    ]}}, ctx, reg);
    document.body.append(instance.root);
    const wiz = instance.node('wiz');
    instance.node('file').value = 'typed.csv';
    wiz.next();
    wiz.next();
    wiz.back();
    ctx.data.t.value = 'Upload';
    assert.equal(instance.node('wiz'), wiz, 'never rebuilt by the bind');
    assert.equal(wiz.root.querySelector('.u2-wizard-title').textContent, 'Upload');
    assert.equal(wiz.currentStep.peek(), 'step-1', 'the user is still on step 2');
    assert.equal(instance.node('file').peek(), 'typed.csv');

    new SpecEditor(instance).apply({op: 'set-prop', node: instance.spec.root.children[2], name: 'title', value: 'Check'});
    const rebuilt = instance.node('wiz');
    assert.notEqual(rebuilt, wiz, 'the patch rebuilt the wizard');
    assert.equal(rebuilt.currentStep.peek(), 'step-1', 'the carry kept the open step');
    rebuilt.goTo('step-2');
    assert.equal(rebuilt.currentStep.peek(), 'step-2', 'and the visited set: a visited step is reachable');
    instance.dispose();
  });

spec('visible, readOnly and required bound on inputs follow their sources in place', async () => {
  const reg = registry();
  const ctx = new SpecContext({data: {show: true, locked: false, needed: false}});
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-form', name: 'form', children: [
    {tag: 'u2-text-input', name: 'name', props: {label: 'Name'},
      bind: {visible: '$.show', readOnly: '$.locked', required: '$.needed'}},
    {tag: 'u2-button', name: 'hidden', props: {text: 'Never', visible: false}},
  ]}}, ctx, reg);
  const input = instance.node('name');
  assert.equal(instance.node('hidden').root.hidden, true, 'a literal visible: false hides at build');
  ctx.data.show.value = false;
  ctx.data.locked.value = true;
  ctx.data.needed.value = true;
  assert.equal(input.visible, false);
  assert.equal(input.readOnly, true);
  assert.equal(input.required, true);
  assert.equal(input.root.querySelector('input').getAttribute('aria-required'), 'true');
  assert.equal(instance.node('name'), input);
  instance.dispose();
});

spec('the setters behind the live props: radio items, table columns, grid properties, cls, icon size, alt',
  async () => {
    const reg = registry();
    const creates = countCreates(reg);
    const ctx = new SpecContext({data: {items: ['a', 'b'], cols: ['name'], props: [{name: 'Title', type: 'string'}],
      cls: 'one', size: 'small', alt: 'Logo', heading: 'Group', empty: 'Nothing'}});
    const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-div-v', name: 'box', bind: {cls: '$.cls'}, children: [
      {tag: 'u2-radio-input', name: 'radio', props: {label: 'Stage'}, bind: {items: '$.items'}},
      {tag: 'u2-data-table', name: 'table', props: {items: [{name: 'Aspirin', cas: '50-78-2'}]}, bind: {columns: '$.cols'}},
      {tag: 'u2-property-grid', name: 'grid', bind: {properties: '$.props'}},
      {tag: 'u2-icon', name: 'star', props: {name: 'star'}, bind: {size: '$.size'}},
      {tag: 'u2-image-input', name: 'logo', props: {label: 'Logo'}, bind: {alt: '$.alt'}},
      {tag: 'u2-access-grid', name: 'access', props: {label: 'Access'}, bind: {principalLabel: '$.heading'}},
      {tag: 'u2-chips-input', name: 'chips', props: {label: 'Chips', items: []}, bind: {emptyText: '$.empty'}},
    ]}}, ctx, reg);
    document.body.append(instance.root);
    const built = creates.n;
    ctx.data.items.value = ['a', 'b', 'c'];
    ctx.data.cols.value = ['name', 'cas:CAS'];
    ctx.data.props.value = [{name: 'Title', type: 'string'}, {name: 'Rows', type: 'int'}];
    ctx.data.cls.value = 'two';
    ctx.data.size.value = 'large';
    ctx.data.alt.value = 'Brand';
    ctx.data.heading.value = 'Who';
    ctx.data.empty.value = 'None yet';
    await flush();
    assert.equal(creates.n, built, 'nothing rebuilt');
    assert.equal(instance.node('radio').root.querySelectorAll('input[type=radio]').length, 3);
    const table = instance.node('table');
    assert.deepEqual([...table.root.querySelectorAll('.u2-data-table-head')].map((h) => h.textContent), ['name', 'CAS']);
    assert.equal(instance.node('grid').root.querySelectorAll('.u2-propgrid-row').length, 2);
    const box = instance.node('box').root;
    assert.equal(box.classList.contains('two') && !box.classList.contains('one'), true);
    assert.equal(instance.node('star').root.classList.contains('u2-icon-large'), true);
    assert.equal(instance.node('logo').root.querySelector('img').alt, 'Brand');
    assert.equal(instance.node('access').root.querySelector('th').textContent, 'Who');
    assert.equal(instance.node('chips').root.querySelector('.u2-chips-empty').textContent, 'None yet');
    instance.dispose();
  });

spec('forward ref on a bindable prop the component holds no signal for: the proxy follows the target', () => {
  const reg = registry();
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-panel', children: [
    {tag: 'u2-stat-card', name: 'card', props: {value: '1'}, bind: {label: '$.later.value'}},
    {tag: 'u2-text-input', name: 'later', props: {label: 'Later', value: 'first'}},
  ]}}, new SpecContext(), reg);
  const label = () => instance.node('card').root.querySelector('[data-u2-part="label"]').textContent;
  assert.equal(label(), 'first');
  instance.node('later').value = 'second';
  assert.equal(label(), 'second');
  instance.dispose();
});
