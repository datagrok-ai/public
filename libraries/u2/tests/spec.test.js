/* dg-ui/1 spec tests: render a spec against a context, edit it from both ends, contain broken
   nodes, and round-trip through dump(). Every test registers its own fakes in a private registry,
   so nothing here depends on which components the library has registered. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Signal} from '../src/core/signals.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {TextInput} from '../src/components/text-input.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, parseSpec, renderSpec} from '../src/spec/spec.js';

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
      {tag: 'u2-fake-input', props: {label: 'Bound'}, bind: {label: '$.name'}},
      {tag: 'u2-fake-input', props: {label: 'Fine', bogus: 'kept'}},
    ]},
  }, new SpecContext(), reg);

  const errors = instance.root.querySelectorAll('.u2-spec-error').map((el) => el.textContent);
  assert.equal(errors.length, 4);
  assert.match(errors[0], /^u2-missing: /);
  assert.match(errors[1], /^blink: /);
  assert.match(errors[2], /prop "label" expects string/);
  assert.match(errors[3], /not bindable/);
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
  assert.match(bad.root.querySelector('.u2-spec-error').textContent, /must be "cmd:<name>"/);
  instance.dispose();
  bad.dispose();
});

spec('dump: round-trips the spec with live values folded in', () => {
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
  assert.equal(dumped.root.children[1].props.value, 'Ibuprofen', 'the edited value is in the dump');
  assert.deepEqual(dumped.root.children[2].bind, {value: '$.alias'}, 'a bound prop stays a binding');
  assert.equal('value' in dumped.root.children[2].props, false);

  const restored = renderSpec(dumped, ctx, reg);
  assert.deepEqual(editors(restored).map((e) => e.value), ['Ibuprofen', 'IBU']);
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

spec('parseSpec: takes JSON or objects and rejects a wrong envelope', () => {
  const source = {$schema: 'dg-ui/1', root: {tag: 'div'}};
  assert.deepEqual(parseSpec(JSON.stringify(source)), source);
  assert.throws(() => parseSpec({$schema: 'dg-ui/2', root: {tag: 'div'}}), /\$schema must be "dg-ui\/1"/);
  assert.throws(() => parseSpec({root: {tag: 'div'}}), /\$schema must be "dg-ui\/1"/);
  assert.throws(() => parseSpec({$schema: 'dg-ui/1'}), /"root" must be a node/);
  assert.throws(() => parseSpec({$schema: 'dg-ui/1', root: {}}), /"root" must be a node/);
  assert.throws(() => parseSpec('not json'), SyntaxError);
});
