/* dg-ui/1 spec tests: render a spec against a context, edit it from both ends, contain broken
   nodes, and round-trip through dump(). Every test registers its own fakes in a private registry,
   so nothing here depends on which components the library has registered. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Signal, computed, signal} from '../src/core/signals.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {TextInput} from '../src/components/inputs/text-input.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, checkProp, parseSpec, renderSpec} from '../src/spec/spec.js';
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

const FAKE_INPUT = {
  tag: 'u2-fake-input',
  description: 'Fake text input for the spec tests',
  props: [
    {name: 'label', type: 'string'},
    {name: 'value', type: 'string', bindable: true, twoWay: true},
  ],
  events: ['input'],
  example: {tag: 'u2-fake-input', props: {label: 'Name'}},
};

/** Takes the bound signal as its own value signal — what a real registration does. */
function registerFake(reg, seen) {
  reg.register({
    ...FAKE_INPUT,
    create: (props) => {
      if (seen)
        seen.push(props);
      const bound = props.value instanceof Signal;
      return new TextInput({
        label: props.label,
        value: bound ? undefined : props.value,
        bind: bound ? props.value : undefined,
      });
    },
  });
}

/** Keeps a signal of its own, so the two-way bridge has something to bridge. */
function registerEcho(reg) {
  reg.register({...FAKE_INPUT, tag: 'u2-fake-echo', create: () => new TextInput({})});
}

function editors(instance) {
  return instance.root.querySelectorAll('input');
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

spec('registry: rejects duplicates and manifests the metadata minus create', () => {
  const reg = new Registry();
  registerFake(reg);
  assert.throws(() => registerFake(reg), /already registered/);

  const manifest = reg.manifest();
  assert.equal(manifest.schema, 'dg-ui/1');
  assert.equal(manifest.components.length, 1);
  const meta = manifest.components[0];
  assert.equal('create' in meta, false);
  assert.deepEqual(meta, FAKE_INPUT);
  assert.equal(reg.get('u2-fake-input').create instanceof Function, true);
  assert.equal(reg.get('u2-nope'), undefined);
});

spec('renderSpec: three levels, a two-way binding and a command', () => {
  const reg = new Registry();
  registerFake(reg);
  let saves = 0;
  const ctx = new SpecContext({data: {name: 'world'}, commands: {save: () => saves++}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {
      tag: 'div',
      props: {cls: 'u2-spec-form'},
      children: [
        {tag: 'h1', props: {text: 'Details'}},
        {tag: 'div', children: [
          {tag: 'u2-fake-input', props: {label: 'Name'}, bind: {value: '$.name'}},
          {tag: 'span', props: {text: 'Save', cls: 'save-btn'}, on: {click: 'cmd:save'}},
        ]},
      ],
    },
  }, ctx, reg);
  document.body.append(instance.root);

  assert.equal(instance.root.querySelector('h1').textContent, 'Details');
  assert.equal(instance.root.querySelector('.u2-spec-form') !== null, true);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);

  const editor = editors(instance)[0];
  assert.equal(editor.value, 'world');
  ctx.data.name.value = 'u2';
  assert.equal(editor.value, 'u2', 'the context signal drives the input');

  type(editor, 'typed');
  assert.equal(ctx.data.name.value, 'typed', 'and the input drives the context signal');

  const save = instance.root.querySelector('.save-btn');
  fire(save, 'click');
  assert.equal(saves, 1);

  instance.dispose();
  fire(save, 'click');
  assert.equal(saves, 1, 'listeners die with the instance');
});

spec('renderSpec: bridges a two-way prop when the component keeps its own signal', () => {
  const reg = new Registry();
  registerEcho(reg);
  const ctx = new SpecContext({data: {name: 'world'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-echo', props: {label: 'Name'}, bind: {value: '$.name'}},
  }, ctx, reg);

  const editor = editors(instance)[0];
  assert.equal(editor.value, 'world', 'the context value wins over whatever create seeded');
  ctx.data.name.value = 'u2';
  assert.equal(editor.value, 'u2');

  type(editor, 'typed');
  assert.equal(ctx.data.name.value, 'typed');
  assert.equal(editor.value, 'typed', 'the echo back into the editor is suppressed');
  instance.dispose();
});

spec('renderSpec: broken nodes become placeholders, siblings still render', () => {
  const reg = new Registry();
  const seen = [];
  registerFake(reg, seen);
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'u2-missing', props: {label: 'Nope'}},
      {tag: 'blink', props: {text: 'retro'}},
      {tag: 'u2-fake-input', props: {label: 42}},
      {tag: 'u2-fake-input', props: {label: 'Bound'}, bind: {ghost: '$.name'}},
      {tag: 'u2-fake-input', props: {label: 'Fine', bogus: 'kept'}},
    ]},
  }, new SpecContext(), reg);

  const errors = instance.root.querySelectorAll('.u2-spec-error').map((el) => el.textContent);
  assert.equal(errors.length, 4);
  assert.match(errors[0], /^u2-missing: /);
  assert.match(errors[1], /^blink: /);
  assert.match(errors[2], /prop "label" expects string/);
  assert.match(errors[3], /has no prop "ghost" to bind/);
  assert.equal(editors(instance).length, 1, 'the healthy sibling rendered');
  assert.equal(seen[seen.length - 1].bogus, 'kept', 'unknown props are warned about, not dropped');
  instance.dispose();
});

spec('adopt: children are rendered first, then handed over one by one with their index', () => {
  const reg = new Registry();
  registerFake(reg);
  const adopted = [];
  reg.register({
    tag: 'u2-fake-box',
    description: 'Fake adopting container for the spec tests',
    props: [],
    acceptsChildren: true,
    create: () => new Control(),
    adopt: (parent, child, index) => {
      adopted.push([child instanceof Control ? 'component' : child.tagName, index]);
      parent.root.append(child instanceof Control ? child.root : child);
    },
    example: {tag: 'u2-fake-box'},
  });
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-box', children: [
      {tag: 'u2-fake-input', props: {label: 'Name', value: 'Aspirin'}},
      {tag: 'span', props: {text: 'static', cls: 'static-note'}},
    ]},
  }, new SpecContext(), reg);

  assert.deepEqual(adopted, [['component', 0], ['SPAN', 1]]);
  assert.equal(editors(instance)[0].value, 'Aspirin', 'the adopted input is live in the tree');
  assert.equal(instance.root.querySelector('.static-note').textContent, 'static');
  instance.dispose();
});

spec('createWithChildren: rendered children and their spec nodes reach the constructor', () => {
  const reg = new Registry();
  registerFake(reg);
  const seen = [];
  reg.register({
    tag: 'u2-fake-strip',
    description: 'Fake container that takes its children at construction',
    props: [{name: 'gap', type: 'int'}],
    childProps: [{name: 'title', type: 'string'}],
    acceptsChildren: true,
    create: () => new Control(),
    createWithChildren: (props, children, nodes) => {
      seen.push({gap: props.gap, titles: nodes.map((n) => n.props?.title),
        ready: children.map((c) => c instanceof Control && c.root.querySelector('input') !== null)});
      const component = new Control();
      for (const child of children)
        component.root.append(child instanceof Control ? child.root : child);
      return component;
    },
    example: {tag: 'u2-fake-strip'},
  });

  const warnings = captureWarnings(() => {
    const instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'u2-fake-strip', props: {gap: 4}, children: [
        {tag: 'u2-fake-input', props: {title: 'One', label: 'Name'}},
        {tag: 'u2-fake-input', props: {label: 'Alias'}},
      ]},
    }, new SpecContext(), reg);
    assert.deepEqual(seen, [{gap: 4, titles: ['One', undefined], ready: [true, true]}]);
    assert.equal(editors(instance).length, 2);
    instance.dispose();
  });
  assert.deepEqual(warnings, [], 'a prop the parent declares is not an unknown prop');

  const bad = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-strip', children: [{tag: 'u2-fake-input', props: {title: 42}}]},
  }, new SpecContext(), reg);
  assert.match(bad.root.querySelector('.u2-spec-error').textContent, /prop "title" expects string/);
  bad.dispose();
});

spec('json props: any JSON payload passes through, functions and cycles do not', () => {
  const reg = new Registry();
  const seen = [];
  reg.register({
    tag: 'u2-fake-json',
    description: 'Fake holder of a JSON payload',
    props: [{name: 'payload', type: 'object'}],
    create: (props) => {
      seen.push(props.payload);
      return new Control();
    },
    example: {tag: 'u2-fake-json', props: {payload: {rows: 3}}},
  });

  const payload = {rows: [1, 'two', true, null], nested: {a: 1}};
  const good = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-json', props: {payload}},
  }, new SpecContext(), reg);
  assert.equal(good.root.querySelectorAll('.u2-spec-error').length, 0);
  assert.deepEqual(seen, [payload], 'the payload reaches create untouched');

  const cycle = {};
  cycle.self = cycle;
  for (const value of [{run: () => {}}, cycle, undefined]) {
    const bad = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'u2-fake-json', props: {payload: value}},
    }, new SpecContext(), reg);
    assert.match(bad.root.querySelector('.u2-spec-error').textContent, /prop "payload" expects object/);
    bad.dispose();
  }
  assert.equal(seen.length, 1, 'a rejected payload never reaches create');
  good.dispose();
});

spec('renderSpec: a missing command warns once instead of throwing', () => {
  const reg = new Registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'span', props: {text: 'Save'}, on: {click: 'cmd:save'}},
  }, new SpecContext(), reg);
  const save = instance.root.querySelector('span');

  const warnings = captureWarnings(() => {
    fire(save, 'click');
    fire(save, 'click');
  });
  assert.equal(warnings.length, 1);
  assert.match(warnings[0], /no command "save"/);

  const bad = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'span', on: {click: 'alert(1)'}},
  }, new SpecContext(), reg);
  assert.match(bad.root.querySelector('.u2-spec-error').textContent, /must be 'cmd:' followed by a command name/);
  instance.dispose();
  bad.dispose();
});

spec('dump: the document as it is — a Run-mode edit never folds in, a bound prop stays a binding', () => {
  const reg = new Registry();
  registerFake(reg);
  const ctx = new SpecContext({data: {alias: 'ASA'}, commands: {touch: () => {}}});
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'div', props: {cls: 'u2-spec-form'}, children: [
      {tag: 'h1', props: {text: 'Compound'}},
      {tag: 'u2-fake-input', props: {label: 'Name', value: 'Aspirin'}},
      {tag: 'u2-fake-input', props: {label: 'Alias'}, bind: {value: '$.alias'}, on: {input: 'cmd:touch'}},
    ]},
  };
  const instance = renderSpec(source, ctx, reg);
  assert.deepEqual(instance.dump(), source, 'an untouched tree dumps back to what it was');

  type(editors(instance)[0], 'Ibuprofen');
  ctx.data.alias.value = 'IBU';
  const dumped = instance.dump();
  assert.equal(dumped.root.children[1].props.value, 'Aspirin', 'the document is the authority, not the live value');
  assert.deepEqual(dumped, source);
  assert.notEqual(dumped.root.children[1].props, source.root.children[1].props, 'a copy, not the document object');
  assert.deepEqual(dumped.root.children[2].bind, {value: '$.alias'}, 'a bound prop stays a binding');
  assert.equal('value' in dumped.root.children[2].props, false);

  const restored = renderSpec(dumped, ctx, reg);
  assert.deepEqual(editors(restored).map((e) => e.value), ['Aspirin', 'IBU'], 'restored from the document');
  assert.deepEqual(restored.dump(), dumped, 'and the restored tree dumps the same again');
  instance.dispose();
  restored.dispose();
});

spec('names: round-trip through dump, and a duplicate warns on the node that repeats it', () => {
  const reg = new Registry();
  registerFake(reg);
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'u2-fake-input', name: 'nameInput', props: {label: 'Name', value: 'Aspirin'}},
      {tag: 'span', name: 'note', props: {text: 'static'}},
    ]},
  };
  const instance = renderSpec(source, new SpecContext(), reg);
  assert.deepEqual(instance.dump(), source, 'names survive the round-trip');
  instance.dispose();

  const warnings = captureWarnings(() => {
    const duplicate = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'div', children: [
        {tag: 'u2-fake-input', name: 'field', props: {label: 'One'}},
        {tag: 'u2-fake-input', name: 'field', props: {label: 'Two'}},
      ]},
    }, new SpecContext(), reg);
    assert.equal(duplicate.root.querySelectorAll('.u2-spec-error').length, 0, 'a duplicate name still renders');
    assert.equal(editors(duplicate).length, 2);
    duplicate.dispose();
  });
  assert.equal(warnings.length, 1);
  assert.match(warnings[0], /duplicate name "field"/);
});

spec('data-u2-name: named nodes carry the automation id, unnamed ones carry nothing', () => {
  const reg = new Registry();
  registerFake(reg);
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'div', name: 'form', children: [
      {tag: 'u2-fake-input', name: 'nameInput', props: {label: 'Name', value: 'Aspirin'}},
      {tag: 'span', props: {text: 'unnamed'}},
      {tag: 'u2-missing', name: 'brokenNode'},
    ]},
  };
  const instance = renderSpec(source, new SpecContext(), reg);
  document.body.append(instance.root);

  assert.equal(instance.root.querySelector('[data-u2-name="form"]'), instance.node('form'));
  assert.equal(instance.root.querySelector('[data-u2-name="nameInput"]'), instance.node('nameInput').root);
  const broken = instance.root.querySelector('[data-u2-name="brokenNode"]');
  assert.equal(broken.classList.contains('u2-spec-error'), true, 'a placeholder is addressable too');
  assert.equal(instance.root.querySelectorAll('[data-u2-name]').length, 3, 'the unnamed span carries none');
  assert.deepEqual(instance.dump(), source, 'stamping leaves the spec alone');
  instance.dispose();
});

spec('nodes, node and nodeAt: every rendered node is reachable, and elements hit-test to theirs', () => {
  const reg = new Registry();
  registerFake(reg);
  const input = {tag: 'u2-fake-input', name: 'nameInput', props: {label: 'Name'}};
  const inner = {tag: 'span', name: 'note', props: {text: 'static', cls: 'note'}};
  const middle = {tag: 'div', props: {cls: 'middle'}, children: [input, inner]};
  const root = {tag: 'div', children: [middle]};
  const instance = renderSpec({$schema: 'dg-ui/1', root}, new SpecContext(), reg);
  document.body.append(instance.root);

  const nodes = instance.nodes();
  assert.equal(nodes.size, 4, 'components and HTML tags alike are tracked');
  assert.equal(nodes.get(input).root.querySelector('input') !== null, true);
  assert.equal(nodes.get(inner).tagName, 'SPAN');
  assert.equal(instance.node('nameInput'), nodes.get(input));
  assert.equal(instance.node('note'), nodes.get(inner));
  assert.equal(instance.node('nope'), undefined);

  assert.equal(instance.nodeAt(editors(instance)[0]), input, 'an element inside a component finds it');
  assert.equal(instance.nodeAt(instance.root.querySelector('.note')), inner);
  assert.equal(instance.nodeAt(instance.root.querySelector('.middle')), middle);
  assert.equal(instance.nodeAt(document.createElement('div')), null, 'nothing outside the spec matches');
  instance.dispose();
});

spec('rerender: the node comes back new in the same place, and the maps follow it', () => {
  const reg = new Registry();
  registerFake(reg);
  const input = {tag: 'u2-fake-input', name: 'nameInput', props: {label: 'Name', value: 'Aspirin'}};
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'span', props: {text: 'before'}},
      input,
      {tag: 'span', props: {text: 'after'}},
    ]},
  }, new SpecContext(), reg);
  document.body.append(instance.root);

  const live = Scope.liveCount;
  const old = instance.nodes().get(input);
  input.props.label = 'Alias';
  instance.rerender(input);

  const built = instance.nodes().get(input);
  assert.notEqual(built, old, 'a fresh component');
  assert.equal(old.scope.isDisposed, true, 'and the one it replaced is disposed');
  assert.equal(Scope.liveCount, live, 'a re-render accumulates nothing');

  const container = instance.root.querySelector('div');
  assert.deepEqual([...container.children].map((el) => el.textContent === 'before' ? 'before' :
    el === built.root ? 'input' : 'after'), ['before', 'input', 'after'], 'same position');
  assert.equal(instance.nodeAt(built.root), input, 'the element map points at the new element');
  assert.equal(instance.node('nameInput'), built, 'and the name at the new component');
  assert.equal(built.root.dataset.u2Name, 'nameInput', 'the automation id is restamped');
  assert.equal(instance.dump().root.children[1].props.label, 'Alias');
  assert.throws(() => instance.rerender({tag: 'div'}), /not rendered/);
  instance.dispose();
});

spec('rerender: a duplicate name does not take the survivor out of node(name) with it', () => {
  const reg = new Registry();
  registerFake(reg);
  const second = {tag: 'u2-fake-input', name: 'field', props: {label: 'Two'}};
  let instance;
  captureWarnings(() => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'div', children: [
        {tag: 'u2-fake-input', name: 'field', props: {label: 'One'}},
        second,
      ]},
    }, new SpecContext(), reg);
  });
  const survivor = instance.node('field');
  captureWarnings(() => instance.rerender(second));
  assert.equal(instance.node('field'), survivor, 'the first node still owns the name');
  assert.notEqual(instance.nodes().get(second), survivor);
  instance.dispose();
});

spec('rerender: a warning fires once per instance, however often the node is rebuilt', () => {
  const reg = new Registry();
  registerFake(reg);
  const node = {tag: 'u2-fake-input', props: {label: 'Name', bogus: 'kept'}};
  let instance;
  const first = captureWarnings(() => {
    instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'div', children: [node]}},
      new SpecContext(), reg);
  });
  assert.equal(first.length, 1);
  assert.match(first[0], /unknown prop "bogus"/);

  const again = captureWarnings(() => {
    instance.rerender(node);
    instance.rerender(node);
  });
  assert.deepEqual(again, [], 'the same message is not re-spammed');
  instance.dispose();
});

spec('checkProp: a literal value for a prop with choices must be one of them; a bind is exempt', () => {
  const direction = {name: 'direction', type: 'string', choices: ['horizontal', 'vertical']};
  assert.equal(checkProp(direction, 'vertical'), null);
  assert.equal(checkProp(direction, 'diagonal'),
    'prop "direction" must be one of "horizontal", "vertical"');
  assert.match(checkProp(direction, 42), /expects string/, 'the type rule still comes first');
  assert.equal(checkProp({name: 'families', type: 'string_list', choices: ['a', 'b']}, ['a', 'c']),
    null, 'the rule is for scalars — on a string_list, choices are the picker\'s, not the array\'s');

  const reg = new Registry();
  reg.register({
    tag: 'u2-fake-choice',
    description: 'Fake choices holder for the spec tests',
    props: [{...direction, bindable: true}],
    create: () => new Control(),
    example: {tag: 'u2-fake-choice'},
  });
  // a bound value is the context's, not the spec's — the renderer never routes it through checkProp
  const bound = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-choice', bind: {direction: '$.dir'}},
  }, new SpecContext({data: {dir: 'diagonal'}}), reg);
  assert.equal(bound.root.querySelectorAll('.u2-spec-error').length, 0);
  bound.dispose();

  const literal = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-choice', props: {direction: 'diagonal'}},
  }, new SpecContext(), reg);
  assert.match(literal.root.querySelector('.u2-spec-error').textContent, /must be one of/);
  literal.dispose();
});

/** The one test here over the real registrations: the designer metadata they declare. */
spec('designer metadata: actions read the node state, and the manifest keeps defaults but not closures', () => {
  const reg = new Registry();
  registerAll(reg);

  const tabs = reg.get('u2-tabs').designerActions;
  assert.deepEqual(tabs.map((a) => a.name), ['Add tab']);
  assert.deepEqual(tabs[0].produce({tag: 'u2-tabs', children: [{tag: 'u2-panel'}, {tag: 'u2-panel'}]}),
    {op: 'add-child', node: {tag: 'u2-panel', props: {title: 'Tab 3'}}},
    'the title counts the children the node has now');
  assert.deepEqual(tabs[0].produce({tag: 'u2-tabs'}),
    {op: 'add-child', node: {tag: 'u2-panel', props: {title: 'Tab 1'}}});
  assert.deepEqual(tabs[0].produce({tag: 'u2-tabs', children: [{tag: 'u2-panel', props: {title: 'Tab 1'}}]}),
    {op: 'add-child', node: {tag: 'u2-panel', props: {title: 'Tab 2'}}},
    'the seeded pane counts: the first Add tab never repeats it');

  assert.deepEqual(reg.get('u2-accordion').designerActions[0]
    .produce({tag: 'u2-accordion', children: [{tag: 'u2-panel'}]}),
  {op: 'add-child', node: {tag: 'u2-panel', props: {title: 'Pane 2'}}});
  assert.deepEqual(reg.get('u2-splitter').designerActions[0].produce({tag: 'u2-splitter'}),
    {op: 'add-child', node: {tag: 'u2-panel'}});

  const crumbs = reg.get('u2-breadcrumbs').designerActions[0];
  assert.deepEqual(crumbs.produce({tag: 'u2-breadcrumbs', props: {items: ['Item 1', 'Item 2']}}),
    {op: 'set-prop', name: 'items', value: ['Item 1', 'Item 2', 'Item 3']},
    'the appended item counts the items the node has now');
  assert.deepEqual(crumbs.produce({tag: 'u2-breadcrumbs'}),
    {op: 'set-prop', name: 'items', value: ['Item 1']});

  const metas = new Map(reg.manifest().components.map((c) => [c.tag, c]));
  assert.deepEqual(metas.get('u2-tabs').defaultChildren, [{tag: 'u2-panel', props: {title: 'Tab 1'}}]);
  assert.deepEqual(metas.get('u2-accordion').defaultChildren,
    [{tag: 'u2-panel', props: {title: 'Pane 1'}}]);
  assert.deepEqual(metas.get('u2-splitter').defaultChildren, [{tag: 'u2-panel'}, {tag: 'u2-panel'}],
    'a one-panel splitter has no sash to drag — it seeds the pair');
  assert.equal('defaultChildren' in metas.get('u2-form'), false,
    'a container that is not a multi-host declares none');
  assert.deepEqual(metas.get('u2-breadcrumbs').defaults, {items: ['Item 1', 'Item 2', 'Item 3']});
  assert.deepEqual(metas.get('u2-choice-input').defaults, {items: ['Item 1', 'Item 2', 'Item 3']});
  assert.deepEqual(metas.get('u2-button').defaults, {text: 'Button'});
  assert.deepEqual(metas.get('u2-splitter').props.find((p) => p.name === 'direction').choices,
    ['horizontal', 'vertical']);
  assert.deepEqual(metas.get('u2-number-input').props.find((p) => p.name === 'mode').choices,
    ['int', 'float']);
  for (const [tag, meta] of metas)
    assert.equal('designerActions' in meta, false, `${tag} leaks designer closures into the manifest`);
});

/* WO-11 — the binding substrate: full path resolution, named nodes as bind sources, writability,
   cycles, and the cmd: tiers with fire-time argument resolution. */

spec('resolveBinding: a named node is a bind source — explicit prop, default step, precedence', () => {
  const reg = new Registry();
  registerEcho(reg);
  // the context deliberately carries a same-named key: the node's declaration is closer and wins
  const ctx = new SpecContext({data: {doseInput: 'shadowed'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'u2-fake-echo', name: 'doseInput', props: {label: 'Dose'}},
      {tag: 'u2-fake-echo', name: 'viaValue', bind: {value: '$.doseInput.value'}},
      {tag: 'u2-fake-echo', name: 'viaDefault', bind: {value: '$.doseInput'}},
    ]},
  }, ctx, reg);
  document.body.append(instance.root);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);

  const source = instance.node('doseInput');
  const resolved = instance.resolveBinding('$.doseInput.value');
  assert.equal(resolved.signal, source.value, 'the node shadows the same-named context key');
  assert.equal(resolved.writable, true);
  assert.equal(instance.resolveBinding('$.doseInput').signal, source.value,
    'the default step is the value signal');
  assert.deepEqual(source.bindProps().map((p) => [p.name, p.writable]), [['value', true]],
    'only signal-backed meta props enumerate');

  const [dose, viaValue, viaDefault] = editors(instance);
  type(dose, 'typed');
  assert.equal(viaValue.value, 'typed', 'two-way bridge over the explicit step');
  assert.equal(viaDefault.value, 'typed', 'and over the default step');
  type(viaValue, 'back');
  assert.equal(dose.value, 'back', 'edits flow back into the named node');
  instance.dispose();
});

spec('resolveBinding: a BindSource in the context is passed through and walked', () => {
  const reg = new Registry();
  const col = signal(42);
  const row = {bindStep: (n) => n === 'weight' ? col : null,
    bindProps: () => [{name: 'weight', type: 'int', writable: true}]};
  const orders = {
    bindStep: (n) => n === 'currentRow' ? row : n === '' ? col : null,
    bindProps: () => [{name: 'currentRow', walkable: true}],
  };
  const ctx = new SpecContext({data: {orders, plain: 'x'}});
  assert.equal(ctx.data.orders, orders, 'passed through un-wrapped');

  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'div'}}, ctx, reg);
  const deep = instance.resolveBinding('$.orders.currentRow.weight');
  assert.equal(deep.signal, col);
  assert.equal(deep.writable, true);
  assert.equal(instance.resolveBinding('$.orders').signal, col,
    'exhausted at the source: the default step');
  assert.throws(() => instance.resolveBinding('$.orders.nope'), /no binding step "nope"/);
  assert.throws(() => instance.resolveBinding('$.orders.currentRow'),
    /needs one more step — one of: weight/);
  assert.throws(() => instance.resolveBinding('$.orders.currentRow.weight.deeper'),
    /the path continues past a signal/);
  assert.throws(() => instance.resolveBinding('$.plain.x'), /the path continues past a signal/);
  assert.throws(() => instance.resolveBinding('$.missing'), /nothing bound at "\$\.missing"/);
  assert.throws(() => instance.resolveBinding('$.orders['), /is not a valid bind path/);
  instance.dispose();
});

spec('resolveBinding: a default step that answers its own source is bounded, not a hang', () => {
  const itself = {bindStep: () => itself, bindProps: () => [{name: 'round', walkable: true}]};
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'div'}},
    new SpecContext({data: {itself}}), new Registry());
  assert.throws(() => instance.resolveBinding('$.itself'), /needs one more step — one of: round/);
  instance.dispose();
});

spec('resolveBinding: a plain HTML node is not a bind source, before or after the node binding it', () => {
  const reg = new Registry();
  registerEcho(reg);
  let instance;
  const warnings = captureWarnings(() => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'div', children: [
        {tag: 'u2-fake-echo', name: 'toLater', bind: {value: '$.later'}},
        {tag: 'span', name: 'static', props: {text: 'x'}},
        {tag: 'u2-fake-echo', name: 'toHtml', bind: {value: '$.static'}},
        {tag: 'span', name: 'later', props: {text: 'y'}},
      ]},
    }, new SpecContext(), reg);
  });
  const errors = instance.root.querySelectorAll('.u2-spec-error').map((el) => el.textContent);
  assert.equal(errors.length, 1);
  assert.match(errors[0], /"static" is not a bind source/, 'a plain HTML node has no signals');
  assert.deepEqual(warnings, ['u2 spec: u2-fake-echo: "later" is not a bind source'],
    'declared later: the node builds, the link is refused at the flush');
  assert.equal(editors(instance).length, 1);
  instance.dispose();
});

/* R-a (VP-26) — forward references: a bind to a visual node declared later links after the render
   pass; document-order-valid specs take exactly the construction-time path they always did. */

const FORWARD = {
  $schema: 'dg-ui/1',
  root: {tag: 'div', children: [
    {tag: 'u2-fake-input', name: 'early', props: {label: 'Early'}, bind: {value: '$.late.value'}},
    {tag: 'u2-fake-input', name: 'late', props: {label: 'Late', value: 'Aspirin'}},
  ]},
};

spec('forward reference: a two-way bind to a later sibling links after the render pass, both ways', () => {
  const reg = new Registry();
  registerFake(reg);
  let instance;
  const warnings = captureWarnings(() => {
    instance = renderSpec(FORWARD, new SpecContext(), reg);
  });
  assert.deepEqual(warnings, []);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  const [early, late] = editors(instance);
  assert.equal(early.value, 'Aspirin', 'seeded from the later node once the pass is over');
  type(late, 'Ibuprofen');
  assert.equal(early.value, 'Ibuprofen', 'the later node drives the earlier one');
  type(early, 'Back');
  assert.equal(late.value, 'Back', 'twoWay: the edit reaches the later node');
  assert.deepEqual(instance.dump(), FORWARD, 'the document is untouched');
  instance.dispose();
});

spec('forward reference: a one-way bind follows the later node and never writes back', () => {
  const reg = new Registry();
  registerFake(reg);
  reg.register({...FAKE_INPUT, tag: 'u2-fake-view',
    props: [{name: 'label', type: 'string'}, {name: 'value', type: 'string', bindable: true}],
    create: (props) => new TextInput({bind: props.value})});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'u2-fake-view', name: 'view', bind: {value: '$.late'}},
      {tag: 'u2-fake-input', name: 'late', props: {value: 'A'}},
    ]},
  }, new SpecContext(), reg);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  const [view, late] = editors(instance);
  assert.equal(view.value, 'A');
  type(late, 'B');
  assert.equal(view.value, 'B', 'the forward direction is live');
  type(view, 'C');
  assert.equal(late.value, 'B', 'nothing flows back');
  type(late, 'D');
  assert.equal(view.value, 'D');
  instance.dispose();
});

spec('forward reference: a target that fails to build leaves the referencing node built, with one warning', () => {
  const reg = new Registry();
  registerFake(reg);
  let instance;
  const warnings = captureWarnings(() => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'div', children: [
        {tag: 'u2-fake-input', name: 'early', bind: {value: '$.broken.value'}},
        {tag: 'u2-fake-input', name: 'broken', props: {label: 42}},
      ]},
    }, new SpecContext(), reg);
  });
  assert.equal(warnings.length, 1);
  assert.match(warnings[0], /not built/);
  assert.equal(instance.node('early') instanceof TextInput, true, 'built, unlinked');
  const errors = instance.root.querySelectorAll('.u2-spec-error').map((el) => el.textContent);
  assert.equal(errors.length, 1, 'the target\'s own placeholder names the real problem');
  assert.match(errors[0], /prop "label" expects string/);
  instance.dispose();
});

spec('forward reference: a loop through a later node is still a cycle — both members placeholders', () => {
  const reg = new Registry();
  registerEcho(reg);
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'u2-fake-echo', name: 'a', bind: {value: '$.c'}},
      {tag: 'u2-fake-echo', name: 'b'},
      {tag: 'u2-fake-echo', name: 'c', bind: {value: '$.a'}},
    ]},
  }, new SpecContext(), reg);
  const errors = instance.root.querySelectorAll('.u2-spec-error').map((el) => el.textContent);
  assert.equal(errors.length, 2);
  assert.match(errors[0], /binding cycle: a → c → a/);
  assert.match(errors[1], /binding cycle: a → c → a/);
  assert.equal(editors(instance).length, 1, 'the node off the loop renders');
  instance.dispose();
});

spec('forward reference: rerender of the referencing node alone links it again', () => {
  const reg = new Registry();
  registerFake(reg);
  const source = JSON.parse(JSON.stringify(FORWARD));
  const instance = renderSpec(source, new SpecContext(), reg);
  const live = Scope.liveCount;
  const old = instance.node('early');
  instance.rerender(source.root.children[0]);
  assert.notEqual(instance.node('early'), old);
  assert.equal(Scope.liveCount, live, 'a re-render accumulates nothing');
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  const [early, late] = editors(instance);
  assert.equal(early.value, 'Aspirin');
  type(late, 'Ibuprofen');
  assert.equal(early.value, 'Ibuprofen', 'the later node still reaches the re-rendered one');
  type(early, 'Back');
  assert.equal(late.value, 'Back');
  instance.dispose();
});

spec('forward reference: a node re-rendered again before the flush links only its live component', () => {
  const reg = new Registry();
  const seen = [];
  registerFake(reg, seen);
  const source = JSON.parse(JSON.stringify(FORWARD));
  const instance = renderSpec(source, new SpecContext(), reg);
  instance._batching = true;
  instance.rerender(source.root);
  const corpse = seen.at(-2).value;
  instance.rerender(source.root.children[0]);
  instance._batching = false;
  instance.rerenderAll([]);
  const [early, late] = editors(instance);
  assert.equal(early.value, 'Aspirin', 'the live component is linked');
  assert.equal(corpse.peek(), undefined, 'the corpse is not');
  type(late, 'Ibuprofen');
  assert.equal(early.value, 'Ibuprofen');
  assert.equal(corpse.peek(), undefined);
  instance.dispose();
});

spec('twoWay to a read-only leaf: the one-way bridge still applies, with a single warning', () => {
  const reg = new Registry();
  registerEcho(reg);
  const base = signal('A');
  const ro = computed(() => base.value.toLowerCase());
  let instance;
  const warnings = captureWarnings(() => {
    instance = renderSpec({
      $schema: 'dg-ui/1',
      root: {tag: 'div', children: [
        {tag: 'u2-fake-echo', name: 'one', bind: {value: '$.ro'}},
        {tag: 'u2-fake-echo', name: 'two', bind: {value: '$.ro'}},
      ]},
    }, new SpecContext({data: {ro}}), reg);
  });
  assert.deepEqual(warnings, ['u2 spec: bind "$.ro" is read-only — edits will not flow back']);
  assert.equal(instance.resolveBinding('$.ro').writable, false);

  const [one, two] = editors(instance);
  assert.equal(one.value, 'a');
  base.value = 'B';
  assert.equal(one.value, 'b', 'the forward direction is live');
  type(two, 'edited');
  assert.equal(ro.value, 'b', 'nothing flows back into the computed');
  assert.equal(one.value, 'b', 'and the sibling keeps the source value');
  instance.dispose();
});

spec('a binding cycle marks its members broken, naming the loop', () => {
  const reg = new Registry();
  registerEcho(reg);
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'u2-fake-echo', name: 'a', bind: {value: '$.b'}},
      {tag: 'u2-fake-echo', name: 'b', bind: {value: '$.a'}},
      {tag: 'u2-fake-echo', name: 'c'},
    ]},
  }, new SpecContext(), reg);
  const errors = instance.root.querySelectorAll('.u2-spec-error').map((el) => el.textContent);
  assert.equal(errors.length, 2);
  assert.match(errors[0], /binding cycle: a → b → a/);
  assert.match(errors[1], /binding cycle: a → b → a/);
  assert.equal(editors(instance).length, 1, 'the node off the loop renders');
  instance.dispose();
});

spec('cmd tiers: bare command, platform fallback, colon, component function, reserved #', () => {
  const reg = new Registry();
  const calls = [];
  class Source extends Control {
    constructor() {
      super();
      this.registerFunction({name: 'refresh', inputs: [],
        apply: (params) => calls.push(['refresh', params])});
    }
  }
  reg.register({tag: 'u2-fake-source', description: 'Fake function owner', props: [],
    create: () => new Source(), example: {tag: 'u2-fake-source'}});

  let saves = 0;
  const platform = [];
  const amount = signal(5);
  const ctx = new SpecContext({
    data: {amount},
    commands: {save: () => saves++},
    callFunction: async (name, args) => platform.push([name, args]),
  });
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'u2-fake-source', name: 'orders'},
      {tag: 'span', props: {cls: 'bare'}, on: {click: 'cmd:save'}},
      {tag: 'span', props: {cls: 'fallback'}, on: {click: 'cmd:Sin'}},
      {tag: 'span', props: {cls: 'platform'}, on: {click: {cmd: 'cmd:MyPkg:SaveOrder',
        args: {amount: '$.amount', tag: '$$.literal', note: 'plain'}}}},
      {tag: 'span', props: {cls: 'component'}, on: {click: {cmd: 'cmd:orders.refresh', args: {n: 1}}}},
      {tag: 'span', props: {cls: 'reserved'}, on: {click: 'cmd:#local'}},
      {tag: 'span', props: {cls: 'nowhere'}, on: {click: {cmd: 'cmd:orders.nope'}}},
    ]},
  }, ctx, reg);
  document.body.append(instance.root);
  const click = (cls) => fire(instance.root.querySelector('.' + cls), 'click');

  click('bare');
  assert.equal(saves, 1, 'a context command wins the bare tier');
  assert.equal(platform.length, 0);

  click('fallback');
  assert.deepEqual(platform, [['Sin', {}]], 'a bare name falls back to the platform tier');

  click('platform');
  assert.deepEqual(platform[1], ['MyPkg:SaveOrder', {amount: 5, tag: '$.literal', note: 'plain'}],
    'args: the path peeked, $$. unescaped once, the literal as-is');
  amount.value = 7;
  click('platform');
  assert.equal(platform[2][1].amount, 7, 'paths resolve at fire time');

  click('component');
  assert.deepEqual(calls, [['refresh', {n: 1}]], 'the dotted tier applies the component function');

  const warnings = captureWarnings(() => {
    click('reserved');
    click('nowhere');
  });
  assert.equal(warnings.length, 2);
  assert.match(warnings[0], /"cmd:#local" is reserved — code-behind commands arrive with dg-form/);
  assert.match(warnings[1], /no function "nope" on "orders"/);
  instance.dispose();
});

spec('dump: structured event entries deep-copy, args never shared with the document', () => {
  const reg = new Registry();
  const entry = {cmd: 'cmd:Pkg:Save', args: {ids: ['a'], note: '$.note'}};
  const source = {$schema: 'dg-ui/1',
    root: {tag: 'span', props: {text: 'Go'}, on: {click: entry}}};
  const instance = renderSpec(source, new SpecContext({data: {note: 'x'}}), reg);
  const dumped = instance.dump();
  assert.deepEqual(dumped.root.on.click, entry);
  assert.notEqual(dumped.root.on.click, entry, 'a copy, not the document object');
  dumped.root.on.click.args.ids.push('b');
  assert.deepEqual(entry.args.ids, ['a'], 'and deep — mutations never reach back');

  const bad = renderSpec({$schema: 'dg-ui/1', root: {tag: 'span', on: {click: {cmd: 'save'}}}},
    new SpecContext(), reg);
  assert.match(bad.root.querySelector('.u2-spec-error').textContent,
    /must be 'cmd:' followed by a command name/);
  instance.dispose();
  bad.dispose();
});

spec('SpecContext.bindProps: data keys typed from current values, sources walkable', () => {
  const base = signal(2);
  const ctx = new SpecContext({data: {
    name: 'x', count: 3, ratio: 1.5, on: true, tags: ['a'], blob: {a: 1},
    ro: computed(() => base.value), src: {bindStep: () => null, bindProps: () => []},
  }});
  assert.deepEqual(ctx.bindProps(), [
    {name: 'name', type: 'string', writable: true},
    {name: 'count', type: 'int', writable: true},
    {name: 'ratio', type: 'double', writable: true},
    {name: 'on', type: 'bool', writable: true},
    {name: 'tags', type: 'string_list', writable: true},
    {name: 'blob', type: 'object', writable: true},
    {name: 'ro', type: 'int', writable: false},
    {name: 'src', walkable: true},
  ]);
});

spec('parseSpec: takes JSON or objects and rejects a wrong envelope', () => {
  const source = {$schema: 'dg-ui/1', root: {tag: 'div'}};
  assert.deepEqual(parseSpec(JSON.stringify(source)), source);
  assert.throws(() => parseSpec({$schema: 'dg-ui/2', root: {tag: 'div'}}), /\$schema must be "dg-ui\/1"/);
  assert.throws(() => parseSpec({root: {tag: 'div'}}), /\$schema must be "dg-ui\/1"/);
  assert.throws(() => parseSpec({$schema: 'dg-ui/1'}), /"root" must be a node/);
  assert.throws(() => parseSpec({$schema: 'dg-ui/1', root: {}}), /"root" must be a node/);
  assert.throws(() => parseSpec('not json'), SyntaxError);
});

/** A property-tier control (VP-1): `level` lives on a plain object, not on a signal member. */
class TierGauge extends Control {
  constructor(target) {
    super();
    this.target = target;
  }

  get propertyTier() {
    return true;
  }

  getProperties() {
    return [{name: 'level', type: 'int', get: () => this.target.level, set: (_t, v) => this.target.level = v}];
  }

  bump() {
    this.fireEvent('bumped', {by: 1});
  }
}

function registerGauge(reg, target) {
  reg.register({
    tag: 'u2-gauge',
    description: 'Property-tier control',
    props: [{name: 'level', type: 'int', bindable: true, twoWay: true}],
    events: ['bumped'],
    create: () => new TierGauge(target),
    example: {tag: 'u2-gauge'},
  });
}

spec('renderSpec: a property-tier prop binds two-way through link', () => {
  const reg = new Registry();
  const target = {level: 1};
  registerGauge(reg, target);
  const ctx = new SpecContext({data: {lvl: 3}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-gauge', name: 'gauge', bind: {level: '$.lvl'}},
  }, ctx, reg);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  assert.equal(target.level, 3, 'the context value reached the object');
  ctx.data.lvl.value = 5;
  assert.equal(target.level, 5);
  instance.node('gauge').bindStep('level').value = 8;
  assert.equal(ctx.data.lvl.value, 8, 'and the step writes back');
  assert.equal(instance.resolveBinding('$.gauge.level').writable, true);
  instance.dispose();
});

spec('forward reference: a property-tier step follows a later input after the flush, both ways', () => {
  const reg = new Registry();
  const target = {level: 1};
  registerGauge(reg, target);
  registerFake(reg);
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'u2-gauge', name: 'gauge', bind: {level: '$.late.value'}},
      {tag: 'u2-fake-input', name: 'late', props: {value: '7'}},
    ]},
  }, new SpecContext(), reg);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  assert.equal(target.level, '7', 'the later input reached the object');
  const [late] = editors(instance);
  type(late, '9');
  assert.equal(target.level, '9');
  instance.node('gauge').bindStep('level').value = '3';
  assert.equal(late.value, '3', 'and the step writes back');
  instance.dispose();
});

spec('on: a DOM event on a u2-button and a component event on a control both run their command', () => {
  const reg = new Registry();
  registerAll(reg);
  registerGauge(reg, {level: 1});
  const ran = [];
  const ctx = new SpecContext({commands: {x: () => ran.push('x'), y: () => ran.push('y')}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'u2-button', name: 'btn', props: {text: 'Go'}, on: {click: 'cmd:x'}},
      {tag: 'u2-gauge', name: 'gauge', on: {bumped: 'cmd:y'}},
    ]},
  }, ctx, reg);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  fire(instance.root.querySelector('[data-u2-name="btn"]'), 'click');
  assert.deepEqual(ran, ['x']);
  const gauge = instance.node('gauge');
  gauge.bump();
  assert.deepEqual(ran, ['x', 'y']);
  instance.dispose();
  gauge.bump();
  assert.deepEqual(ran, ['x', 'y'], 'listeners die with the instance');
});
