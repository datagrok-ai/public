/* The patch engine (WO-8): every gesture the designer will make is one of these patches, and this
   is where they are proven — what applies, what is refused, what the re-render touches, and that
   undo puts the document back exactly where it was. Platform-free, so it runs headless as is. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {TextInput} from '../src/components/text-input.js';
import {TabStrip} from '../src/components/tabs.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {SpecEditor, nameForTag, renameRefs, seedNode} from '../src/spec/editor.js';

/** Every test runs against a clean document and must leave the live-scope count where it was. */
function edit(name, body) {
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

/** Every component the fakes build, in construction order — how a test sees what re-rendered. */
const created = [];

function element(built) {
  return built instanceof Control ? built.root : built;
}

function register(reg) {
  reg.register({
    tag: 'u2-fake-input',
    description: 'Fake input',
    props: [
      {name: 'label', type: 'string'},
      {name: 'size', type: 'int'},
      {name: 'value', type: 'string', bindable: true, twoWay: true},
    ],
    events: ['input'],
    create: (props) => {
      created.push('u2-fake-input');
      return new TextInput({label: props.label, value: props.value});
    },
    example: {tag: 'u2-fake-input'},
  });
  reg.register({
    tag: 'u2-fake-box',
    description: 'Fake container',
    props: [],
    acceptsChildren: true,
    create: () => {
      created.push('u2-fake-box');
      return new Control();
    },
    example: {tag: 'u2-fake-box'},
  });
  reg.register({
    tag: 'u2-fake-src',
    description: 'Fake source — the object params prop takes sub-binds',
    props: [
      {name: 'func', type: 'string'},
      {name: 'params', type: 'object', subBindable: true},
    ],
    create: () => {
      created.push('u2-fake-src');
      return new Control();
    },
    example: {tag: 'u2-fake-src'},
  });
  reg.register({
    tag: 'u2-fake-leaf',
    description: 'Fake leaf',
    props: [],
    create: () => {
      created.push('u2-fake-leaf');
      return new Control();
    },
    example: {tag: 'u2-fake-leaf'},
  });
  reg.register({
    tag: 'u2-fake-form',
    description: 'Fake form — wires every child into itself',
    props: [],
    acceptsChildren: true,
    create: () => {
      created.push('u2-fake-form');
      return new Control();
    },
    adopt: (parent, child) => parent.root.append(element(child)),
    example: {tag: 'u2-fake-form'},
  });
  reg.register({
    tag: 'u2-fake-tabs',
    description: 'Fake tab strip — takes its children at construction',
    props: [],
    childProps: [{name: 'title', type: 'string'}],
    acceptsChildren: true,
    create: () => new Control(),
    createWithChildren: (_props, children, nodes) => {
      created.push('u2-fake-tabs');
      const strip = new Control();
      for (let i = 0; i < children.length; i++) {
        const label = document.createElement('span');
        label.className = 'tab-label';
        label.textContent = nodes[i].props?.title ?? '';
        strip.root.append(label);
        strip.root.append(element(children[i]));
      }
      return strip;
    },
    example: {tag: 'u2-fake-tabs'},
  });
  reg.register({
    tag: 'u2-lazy-tabs',
    description: 'The real TabStrip — a pane\'s content attaches only when its tab first activates',
    props: [],
    childProps: [{name: 'title', type: 'string'}],
    acceptsChildren: true,
    create: () => new TabStrip(),
    createWithChildren: (_props, children, nodes) => {
      const strip = new TabStrip();
      for (let i = 0; i < children.length; i++)
        strip.addTab({id: `tab-${i}`, label: nodes[i].props?.title ?? `Tab ${i + 1}`,
          content: element(children[i])});
      return strip;
    },
    example: {tag: 'u2-lazy-tabs'},
  });
}

function fresh() {
  return {
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-box', name: 'layout', children: [
      {tag: 'u2-fake-input', name: 'nameInput', props: {label: 'Name', value: 'Aspirin'}},
      {tag: 'u2-fake-box', name: 'details', children: [
        {tag: 'span', name: 'note', props: {text: 'static'}},
      ]},
      {tag: 'u2-fake-leaf', name: 'leaf'},
    ]},
  };
}

function open(source = fresh(), data, commands) {
  const reg = new Registry();
  register(reg);
  created.length = 0;
  const instance = renderSpec(source, new SpecContext({data, commands}), reg);
  document.body.append(instance.root);
  return {source, instance, editor: new SpecEditor(instance)};
}

/** The nodes the DOM shows under `node`, by their automation id. */
function domNames(instance, node) {
  return [...element(instance.nodes().get(node)).children]
    .map((el) => el.dataset.u2Name ?? el.tagName);
}

function find(source, name) {
  const walk = (node) => {
    if (node.name === name)
      return node;
    for (const child of node.children ?? []) {
      const found = walk(child);
      if (found)
        return found;
    }
    return null;
  };
  return walk(source.root);
}

function errors(instance) {
  return instance.root.querySelectorAll('.u2-spec-error').length;
}

edit('add: the node lands where the patch says, and undo/redo move the same object back and forth', () => {
  const {source, instance, editor} = open();
  const before = instance.dump();
  const node = {tag: 'u2-fake-input', name: 'aliasInput', props: {label: 'Alias'}};

  editor.apply({op: 'add', parent: source.root, index: 1, node});
  assert.deepEqual(source.root.children.map((c) => c.name), ['nameInput', 'aliasInput', 'details', 'leaf']);
  assert.deepEqual(domNames(instance, source.root), ['nameInput', 'aliasInput', 'details', 'leaf']);
  assert.equal(instance.node('aliasInput') !== undefined, true);
  assert.equal(editor.canUndo.value, true);

  assert.equal(editor.undo(), true);
  assert.deepEqual(instance.dump(), before, 'the document is exactly what it was');
  assert.deepEqual(domNames(instance, source.root), ['nameInput', 'details', 'leaf']);
  assert.equal(instance.nodes().has(node), false);

  assert.equal(editor.redo(), true);
  assert.equal(instance.nodes().has(node), true, 'the very same node object came back');
  assert.equal(editor.canRedo.value, false);
  instance.dispose();
});

edit('remove: the node and its subtree leave the maps, and undo restores the document', () => {
  const {source, instance, editor} = open();
  const details = find(source, 'details');
  const note = find(source, 'note');
  const before = instance.dump();

  editor.apply({op: 'remove', node: details});
  assert.deepEqual(domNames(instance, source.root), ['nameInput', 'leaf']);
  assert.equal(instance.nodes().has(details), false);
  assert.equal(instance.nodes().has(note), false, 'the subtree went with it');
  assert.equal(instance.node('note'), undefined);

  assert.equal(editor.undo(), true);
  assert.deepEqual(instance.dump(), before);
  assert.equal(instance.nodes().has(details), true);
  assert.equal(instance.node('note') !== undefined, true);
  instance.dispose();
});

edit('move: DOM order follows children order, into another parent and back', () => {
  const {source, instance, editor} = open();
  const input = find(source, 'nameInput');
  const details = find(source, 'details');
  const before = instance.dump();

  editor.apply({op: 'move', node: input, parent: details, index: 0});
  assert.deepEqual(domNames(instance, source.root), ['details', 'leaf']);
  assert.deepEqual(domNames(instance, details), ['nameInput', 'note']);

  editor.undo();
  assert.deepEqual(instance.dump(), before);
  assert.deepEqual(domNames(instance, source.root), ['nameInput', 'details', 'leaf']);

  editor.apply({op: 'move', node: input, parent: source.root, index: 2});
  assert.deepEqual(domNames(instance, source.root), ['details', 'leaf', 'nameInput']);
  editor.undo();
  assert.deepEqual(domNames(instance, source.root), ['nameInput', 'details', 'leaf']);
  instance.dispose();
});

edit('promotion: a parent that wires its children in is rebuilt whole; a plain one is not', () => {
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-box', name: 'layout', children: [
      {tag: 'u2-fake-tabs', name: 'tabs', children: [
        {tag: 'u2-fake-box', name: 'first', props: {title: 'One'},
          children: [{tag: 'u2-fake-input', name: 'inTab'}]},
      ]},
      {tag: 'u2-fake-form', name: 'form', children: [{tag: 'u2-fake-input', name: 'inForm'}]},
      {tag: 'u2-fake-input', name: 'plain', props: {label: 'Plain'}},
    ]},
  };
  const {instance, editor} = open(source);
  const tabs = instance.nodes().get(find(source, 'tabs'));
  const form = instance.nodes().get(find(source, 'form'));
  const layout = instance.nodes().get(source.root);

  editor.apply({op: 'set-prop', node: find(source, 'first'), name: 'title', value: 'Renamed'});
  assert.notEqual(instance.nodes().get(find(source, 'tabs')), tabs, 'the per-child prop rebuilt the strip');
  assert.equal(instance.root.querySelector('.tab-label').textContent, 'Renamed');

  editor.apply({op: 'set-prop', node: find(source, 'inForm'), name: 'label', value: 'Edited'});
  assert.notEqual(instance.nodes().get(find(source, 'form')), form, 'the form was rebuilt around its child');

  created.length = 0;
  editor.apply({op: 'set-prop', node: find(source, 'plain'), name: 'label', value: 'Edited'});
  assert.deepEqual(created, ['u2-fake-input'], 'a plain container re-renders the child alone');
  assert.equal(instance.nodes().get(source.root), layout, 'and keeps its own identity');
  instance.dispose();
});

edit('promotion: a move between two constructing parents rebuilds both, within one rebuilds it once', () => {
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-box', name: 'layout', children: [
      {tag: 'u2-fake-tabs', name: 'left', children: [
        {tag: 'u2-fake-box', name: 'a', props: {title: 'A'}},
        {tag: 'u2-fake-box', name: 'b', props: {title: 'B'}},
      ]},
      {tag: 'u2-fake-tabs', name: 'right', children: [
        {tag: 'u2-fake-box', name: 'c', props: {title: 'C'}},
      ]},
    ]},
  };
  const {instance, editor} = open(source);

  created.length = 0;
  editor.apply({op: 'move', node: find(source, 'a'), parent: find(source, 'right'), index: 0});
  assert.equal(created.filter((tag) => tag === 'u2-fake-tabs').length, 2, 'both strips came back');
  assert.deepEqual(domNames(instance, find(source, 'right')), ['SPAN', 'a', 'SPAN', 'c']);

  created.length = 0;
  editor.apply({op: 'move', node: find(source, 'a'), parent: find(source, 'right'), index: 1});
  assert.equal(created.filter((tag) => tag === 'u2-fake-tabs').length, 1, 'one parent, one rebuild');
  assert.deepEqual(domNames(instance, find(source, 'right')), ['SPAN', 'c', 'SPAN', 'a']);
  instance.dispose();
});

edit('canApply: the refusal matrix — structure, names, prop types, bind and event shape', () => {
  const {source, instance, editor} = open();
  const input = find(source, 'nameInput');
  const details = find(source, 'details');
  const leaf = find(source, 'leaf');
  const node = {tag: 'u2-fake-input'};

  assert.match(editor.canApply({op: 'add', parent: leaf, index: 0, node}), /takes no children/);
  assert.match(editor.canApply({op: 'add', parent: details, index: 5, node}), /out of range/);
  assert.equal(editor.canApply({op: 'add', parent: details, index: 1, node}), null);
  assert.match(editor.canApply({op: 'remove', node: source.root}), /root node cannot be/);
  assert.match(editor.canApply({op: 'move', node: source.root, parent: details, index: 0}),
    /root node cannot be/);
  assert.match(editor.canApply({op: 'move', node: details, parent: leaf, index: 0}),
    /takes no children/);
  assert.match(editor.canApply({op: 'move', node: details, parent: details, index: 0}),
    /cannot move into itself/);
  assert.match(editor.canApply({op: 'move', node: details, parent: find(source, 'note'), index: 0}),
    /cannot move into itself/);
  assert.match(editor.canApply({op: 'set-prop', node: input, name: 'size', value: 1.5}),
    /prop "size" expects int/);
  assert.match(editor.canApply({op: 'set-prop', node: input, name: 'nope', value: 'x'}),
    /has no prop "nope"/);
  assert.equal(editor.canApply({op: 'set-prop', node: input, name: 'label', value: undefined}), null);
  assert.match(editor.canApply({op: 'set-bind', node: input, name: 'label', path: '$.x'}),
    /not bindable/);
  assert.match(editor.canApply({op: 'set-bind', node: input, name: 'value', path: 'x'}),
    /must start with "\$\."/);
  const eventRefusal = editor.canApply({op: 'set-event', node: input, event: 'input', command: 'alert(1)'});
  assert.equal(eventRefusal,
    'an event must name a command as \'cmd:\' followed by the command name — got \'alert(1)\'',
    'complete and free of angle brackets — the platform balloon sanitizer eats <...> with everything after it');
  assert.match(editor.canApply({op: 'rename', node: input, name: '2bad'}), /is not a valid name/);
  assert.match(editor.canApply({op: 'rename', node: input, name: 'details'}), /already taken/);
  assert.equal(editor.canApply({op: 'rename', node: input, name: 'nameInput'}), null,
    'a node keeps its own name');
  assert.throws(() => editor.apply({op: 'remove', node: source.root}), /u2 editor: /);
  assert.equal(editor.canUndo.value, false, 'a refused patch is not history');

  const html = find(source, 'note');
  assert.equal(editor.canApply({op: 'set-prop', node: html, name: 'text', value: 'x'}), null);
  assert.match(editor.canApply({op: 'set-prop', node: html, name: 'src', value: 'x'}),
    /has no prop "src"/);
  instance.dispose();
});

edit('a patch that breaks the node is contained: it applies, renders the placeholder, and undoes', () => {
  const {source, instance, editor} = open(fresh(), {here: 'ok'});
  const input = find(source, 'nameInput');
  assert.equal(errors(instance), 0);

  editor.apply({op: 'set-bind', node: input, name: 'value', path: '$.nowhere'});
  assert.equal(errors(instance), 1);
  assert.equal(instance.nodes().has(input), true, 'and stays selectable');
  assert.equal(instance.dump().root.children[0].bind.value, '$.nowhere');

  editor.apply({op: 'set-bind', node: input, name: 'value', path: '$.here'});
  assert.equal(errors(instance), 0, 'a live binding brings the component back');
  assert.equal(instance.node('nameInput') instanceof TextInput, true);

  editor.undo();
  editor.undo();
  assert.equal(errors(instance), 0);
  assert.equal('bind' in instance.dump().root.children[0], false, 'an emptied container goes too');
  instance.dispose();
});

edit('coalescing: one field is one undo step, and an undone step never merges', () => {
  const {source, instance, editor} = open();
  const input = find(source, 'nameInput');
  const set = (value) => editor.apply({op: 'set-prop', node: input, name: 'label', value});

  set('A');
  set('B');
  set('C');
  assert.equal(input.props.label, 'C');
  assert.equal(editor.undo(), true);
  assert.equal(input.props.label, 'Name', 'three keystrokes, one step back to the original');
  assert.equal(editor.canUndo.value, false);

  editor.redo();
  set('D');
  editor.apply({op: 'set-prop', node: input, name: 'size', value: 4});
  set('E');
  editor.undo();
  assert.equal(input.props.label, 'D', 'a patch on another key broke the merge');
  editor.undo();
  assert.equal(input.props.size, undefined);
  editor.undo();
  assert.equal(input.props.label, 'C', 'and the redone step was its own entry');
  instance.dispose();
});

edit('dirty: position against the saved mark, not bare canUndo', () => {
  const {source, instance, editor} = open();
  const input = find(source, 'nameInput');
  const set = (value) => editor.apply({op: 'set-prop', node: input, name: 'label', value});

  assert.equal(editor.dirty.value, false, 'a fresh document is clean');
  set('A');
  assert.equal(editor.dirty.value, true);
  editor.markClean();
  assert.equal(editor.dirty.value, false, 'saved — clean while canUndo stays true');
  assert.equal(editor.canUndo.value, true);

  editor.undo();
  assert.equal(editor.dirty.value, true, 'undoing past the mark is a different document');
  editor.redo();
  assert.equal(editor.dirty.value, false, 'redoing back to the mark is clean again');

  // a keystroke that COALESCES into the marked entry still changes the document
  editor.apply({op: 'set-prop', node: input, name: 'size', value: 4});
  editor.markClean();
  assert.equal(editor.dirty.value, false);
  editor.apply({op: 'set-prop', node: input, name: 'size', value: 5});
  assert.equal(input.props.size, 5);
  assert.equal(editor.dirty.value, true, 'a merged keystroke dirties the marked entry');
  editor.undo();
  assert.equal(input.props.size, undefined, 'the merge is still one undo step');
  assert.equal(editor.dirty.value, true, 'and the marked mid-word state is simply gone');
  instance.dispose();
});

edit('rename: first-segment references move with the name, and nothing else does', () => {
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-box', name: 'layout', children: [
      {tag: 'u2-fake-input', name: 'old', props: {label: 'Old'}},
      {tag: 'u2-fake-input', name: 'consumer', bind: {value: '$.old'},
        on: {input: 'cmd:old.run', change: 'cmd:old'}},
      {tag: 'u2-fake-input', name: 'deep', bind: {value: '$.old.x'}, on: {input: 'cmd:other.old'}},
      {tag: 'u2-fake-input', name: 'bracket', bind: {value: '$.old[\'a b\']'}},
      {tag: 'u2-fake-input', name: 'sibling', bind: {value: '$.older'}},
    ]},
  };
  const {instance, editor} = open(source, {old: 'value', older: 'other'});
  const before = instance.dump();

  editor.apply({op: 'rename', node: find(source, 'old'), name: 'renamed'});
  const dumped = instance.dump();
  assert.deepEqual(dumped.root.children.map((c) => c.bind?.value),
    [undefined, '$.renamed', '$.renamed.x', '$.renamed[\'a b\']', '$.older']);
  assert.deepEqual(dumped.root.children[1].on, {input: 'cmd:renamed.run', change: 'cmd:old'},
    'a bare cmd: names a context command, not a node');
  assert.deepEqual(dumped.root.children[2].on, {input: 'cmd:other.old'});
  assert.equal(instance.node('renamed') !== undefined, true);
  assert.equal(instance.node('old'), undefined);

  editor.undo();
  assert.deepEqual(instance.dump(), before, 'the inverse rewrites every reference back');

  assert.match(editor.canApply({op: 'rename', node: find(source, 'old'), name: 'consumer'}),
    /already taken/);
  instance.dispose();
});

edit('renameRefs: the grammar decides what a first segment is, and it answers what it touched', () => {
  const child = {tag: 'div', bind: {a: '$.older'}};
  const root = {tag: 'div', bind: {a: '$.old', b: '$.old.x', c: '$.old[\'a b\']', d: '$.older',
    e: '$.x.old'}, on: {i: 'cmd:old.run', j: 'cmd:old', k: 'cmd:oldish.run'}, children: [child]};
  assert.deepEqual(renameRefs(root, 'old', 'next'), [root], 'the child references nothing that moved');
  assert.deepEqual(root.bind, {a: '$.next', b: '$.next.x', c: '$.next[\'a b\']', d: '$.older',
    e: '$.x.old'});
  assert.deepEqual(root.on, {i: 'cmd:next.run', j: 'cmd:old', k: 'cmd:oldish.run'});
  assert.deepEqual(renameRefs({tag: 'div'}, 'old', 'next'), []);
});

edit('rename: a node whose references were rewritten is rendered again, and so is undoing it', async () => {
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-box', name: 'layout', children: [
      {tag: 'u2-fake-input', name: 'old', props: {label: 'Old', value: 'A'}},
      {tag: 'u2-fake-input', name: 'consumer', bind: {value: '$.old'}},
    ]},
  };
  // `$.old` is a node-to-node bind — the named input is the bind source — so what the consumer
  // is wired to has to follow the document through the rename, in both directions
  const {instance, editor} = open(source);
  const consumer = find(source, 'consumer');
  await flush();
  const built = instance.nodes().get(consumer);
  assert.equal(built.value.peek(), 'A', 'bound through the named node\'s default step');
  assert.equal(errors(instance), 0);

  editor.apply({op: 'rename', node: find(source, 'old'), name: 'renamed'});
  await flush();
  const rewritten = instance.nodes().get(consumer);
  assert.notEqual(rewritten, built, 'the node whose bind path the rename rewrote came back');
  assert.equal(errors(instance), 0, 'and the moved path still resolves');
  const renamed = instance.nodes().get(find(source, 'renamed'));
  renamed.value.value = 'B';
  assert.equal(rewritten.value.peek(), 'B', 'wired to the renamed node, not to a corpse');

  editor.undo();
  await flush();
  const undone = instance.nodes().get(consumer);
  assert.notEqual(undone, rewritten, 'and the reverse rewrite renders it again too');
  assert.equal(errors(instance), 0);
  assert.equal(undone.value.peek(), 'A', 'seeded from the re-rendered source once more');
  instance.dispose();
});

edit('rename: a context data key or a command is not a free name — the rewrite would capture it', () => {
  const {source, instance, editor} = open(fresh(), {orders: 'x'}, {refresh: () => {}});
  const input = find(source, 'nameInput');

  assert.equal(editor.canApply({op: 'rename', node: input, name: 'orders'}),
    '"orders" is already a context data key or command');
  assert.equal(editor.canApply({op: 'rename', node: input, name: 'refresh'}),
    '"refresh" is already a context data key or command');
  assert.equal(editor.canApply({op: 'rename', node: input, name: 'toString'}), null,
    'what every object inherits is nobody\'s name');
  assert.throws(() => editor.apply({op: 'rename', node: input, name: 'orders'}), /u2 editor: /);
  assert.equal(input.name, 'nameInput');
  instance.dispose();
});

edit('uniqueName, nameForTag and duplicate: what the palette and Duplicate name things', () => {
  const {source, instance, editor} = open();
  assert.equal(nameForTag('u2-text-input'), 'textInput');
  assert.equal(nameForTag('u2-form'), 'form');
  assert.equal(nameForTag('div'), 'div');

  assert.equal(editor.uniqueName('textInput'), 'textInput1');
  assert.equal(editor.uniqueName('nameInput'), 'nameInput1');

  const unique = editor.uniqueNames();
  assert.deepEqual([unique('textInput'), unique('textInput'), unique('note')],
    ['textInput1', 'textInput2', 'note1'], 'a run reserves what it hands out');
  assert.equal(editor.uniqueName('textInput'), 'textInput1',
    'and reserves nothing in the document — only what a single seeding hands out');

  const details = find(source, 'details');
  const clone = editor.duplicate(details);
  assert.equal(clone.name, 'details1');
  assert.equal(clone.children[0].name, 'note1');
  assert.notEqual(clone.children[0], details.children[0], 'a deep copy, not a shared subtree');
  assert.equal(instance.nodes().has(clone), false, 'and detached until an add patch takes it');

  editor.apply({op: 'add', parent: source.root, index: 3, node: clone});
  assert.equal(editor.duplicate(details).name, 'details2');
  assert.deepEqual(domNames(instance, source.root), ['nameInput', 'details', 'leaf', 'details1']);
  instance.dispose();
});

edit('nodePath and nodeAt: the serializable address a patch from outside arrives with', () => {
  const {source, instance, editor} = open();
  assert.deepEqual(editor.nodePath(source.root), []);
  assert.deepEqual(editor.nodePath(find(source, 'details')), [1]);
  assert.deepEqual(editor.nodePath(find(source, 'note')), [1, 0]);
  assert.equal(editor.nodeAt([1, 0]), find(source, 'note'));
  assert.equal(editor.nodeAt([]), source.root);
  assert.equal(editor.nodeAt([9]), null);
  assert.equal(editor.nodePath({tag: 'div'}), null, 'a detached node has no address');
  instance.dispose();
});

edit('the undo stack is capped: entry 201 pushes the oldest one off', () => {
  const {source, instance, editor} = open();
  const input = find(source, 'nameInput');
  for (let i = 1; i <= 201; i++) {
    editor.apply(i % 2 === 1 ?
      {op: 'set-prop', node: input, name: 'label', value: `v${i}`} :
      {op: 'set-prop', node: input, name: 'size', value: i});
  }
  let undone = 0;
  while (editor.undo())
    undone++;
  assert.equal(undone, 200);
  assert.equal(input.props.label, 'v1', 'the first patch has no inverse left to undo it');
  assert.equal(editor.canUndo.value, false);
  instance.dispose();
});

edit('seedNode: defaults deep-cloned, the label humanized off the name, and defaults winning', () => {
  const crumbs = {
    tag: 'u2-fake-crumbs',
    description: 'Fake crumbs',
    props: [{name: 'items', type: 'string_list'}],
    defaults: {items: ['Item 1', 'Item 2']},
    create: () => new Control(),
    example: {tag: 'u2-fake-crumbs'},
  };
  const a = seedNode(crumbs, 'u2-fake-crumbs', 'crumbs1');
  const b = seedNode(crumbs, 'u2-fake-crumbs', 'crumbs2');
  assert.deepEqual(a, {tag: 'u2-fake-crumbs', name: 'crumbs1', props: {items: ['Item 1', 'Item 2']}});
  assert.notEqual(a.props.items, crumbs.defaults.items, 'cloned off the meta');
  assert.notEqual(a.props.items, b.props.items, 'and never shared between two inserts');

  const labeled = {...crumbs, props: [{name: 'label', type: 'string'}], defaults: undefined};
  assert.deepEqual(seedNode(labeled, 'u2-fake-input', 'textInput1'),
    {tag: 'u2-fake-input', name: 'textInput1', props: {label: 'Text input 1'}},
    'camelCase and the trailing counter split into sentence case');
  assert.deepEqual(seedNode(labeled, 'u2-fake-input', 'form2').props, {label: 'Form 2'});

  assert.equal(seedNode({...labeled, defaults: {label: 'Preset'}}, 'u2-fake-input', 'textInput1')
    .props.label, 'Preset', 'a default label wins over the humanization');
  assert.deepEqual(seedNode({...crumbs, defaults: undefined}, 'u2-fake-crumbs', 'crumbs1'),
    {tag: 'u2-fake-crumbs', name: 'crumbs1'}, 'no label prop, no defaults — no props at all');
  assert.deepEqual(seedNode(undefined, 'div', 'div1'), {tag: 'div', name: 'div1'},
    'an unregistered tag seeds bare');
});

edit('seedNode: a multi-host seeds its default children, cloned per insert and uniquely named', () => {
  const {source, instance, editor} = open();
  const host = {
    tag: 'u2-fake-host',
    description: 'Fake multi-host',
    props: [],
    childProps: [{name: 'title', type: 'string'}],
    acceptsChildren: true,
    defaultChildren: [{tag: 'u2-fake-box', props: {title: 'Pane 1'}}, {tag: 'u2-fake-box'}],
    create: () => new Control(),
    example: {tag: 'u2-fake-host'},
  };
  // the document already carries the first free name the seeding would reach for
  editor.apply({op: 'rename', node: find(source, 'details'), name: 'fakeBox1'});

  const seed = () => {
    const unique = editor.uniqueNames();
    return seedNode(host, 'u2-fake-host', unique(nameForTag('u2-fake-host')), unique);
  };
  const a = seed();
  assert.deepEqual(a, {tag: 'u2-fake-host', name: 'fakeHost1', children: [
    {tag: 'u2-fake-box', props: {title: 'Pane 1'}, name: 'fakeBox2'},
    {tag: 'u2-fake-box', name: 'fakeBox3'},
  ]}, 'siblings of one insert never collide, and neither does the document');

  const b = seed();
  assert.notEqual(a.children[0], host.defaultChildren[0], 'cloned off the meta');
  assert.notEqual(a.children[0].props, b.children[0].props, 'and never shared between two inserts');
  a.children[0].props.title = 'Renamed';
  assert.equal(host.defaultChildren[0].props.title, 'Pane 1', 'editing an insert never edits the meta');

  assert.equal('children' in seedNode({...host, defaultChildren: undefined}, 'u2-fake-host',
    'fakeHost9', editor.uniqueNames()), false, 'a registration that declares none seeds none');
  instance.dispose();
});

edit('a re-rendered node keeps its command wiring — alone, promoted through a form, and on undo', () => {
  const runs = [];
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-box', name: 'layout', children: [
      {tag: 'u2-fake-form', name: 'form', children: [
        {tag: 'u2-fake-input', name: 'inner', props: {label: 'Inner'}, on: {input: 'cmd:ping'}},
      ]},
      {tag: 'u2-fake-input', name: 'outer', props: {label: 'Outer'}, on: {input: 'cmd:ping'}},
    ]},
  };
  const {instance, editor} = open(source, undefined, {ping: () => runs.push('ping')});
  const poke = (name) => fire(element(instance.nodes().get(find(source, name))), 'input');

  poke('outer');
  assert.equal(runs.length, 1, 'wired as rendered');

  editor.apply({op: 'set-prop', node: find(source, 'outer'), name: 'label', value: 'Outer 2'});
  poke('outer');
  assert.equal(runs.length, 2, 'the element a patch rebuilt is wired');

  editor.apply({op: 'set-prop', node: find(source, 'inner'), name: 'label', value: 'Inner 2'});
  poke('inner');
  assert.equal(runs.length, 3, 'a rebuild promoted through the form is wired');

  editor.undo();
  poke('inner');
  editor.undo();
  poke('outer');
  assert.equal(runs.length, 5, 'undo re-renders stay wired too');
  instance.dispose();
});

/* P3.5 WO-5 (F5) — compound patches: applyAll is ONE undo entry, unwinds whole on a mid-sequence
   refusal, publishes once, and never coalesces with the keystroke merging. */

edit('applyAll: one undo entry — undo replays the inverses in reverse, redo goes forward', () => {
  const {source, instance, editor} = open();
  const input = find(source, 'nameInput');
  const details = find(source, 'details');
  const before = instance.dump();

  // adjacent siblings: each inverse re-adds at the index the member had when it was removed, so
  // only a reverse-order replay restores the original order
  editor.apply({op: 'set-prop', node: input, name: 'label', value: 'Edited'});
  editor.applyAll([{op: 'remove', node: input}, {op: 'remove', node: details}]);
  assert.deepEqual(domNames(instance, source.root), ['leaf']);

  assert.equal(editor.undo(), true);
  assert.deepEqual(domNames(instance, source.root), ['nameInput', 'details', 'leaf'],
    'one undo brings the whole compound back, in the original order');
  assert.equal(input.props.label, 'Edited', 'and only the compound — not the entry before it');
  assert.equal(editor.undo(), true);
  assert.deepEqual(instance.dump(), before);
  assert.equal(editor.canUndo.value, false);

  assert.equal(editor.redo(), true);
  assert.equal(editor.redo(), true);
  assert.deepEqual(domNames(instance, source.root), ['leaf'], 'redo replays the compound forward');
  instance.dispose();
});

edit('applyAll: a mid-sequence refusal unwinds the applied members and leaves no history', () => {
  const {source, instance, editor} = open();
  const details = find(source, 'details');
  const note = find(source, 'note');
  const before = instance.dump();

  // removing `details` detaches `note` with it, so the second member is refused mid-sequence
  assert.throws(() => editor.applyAll([{op: 'remove', node: details}, {op: 'remove', node: note}]),
    /u2 editor: /);
  assert.deepEqual(instance.dump(), before, 'the applied member was unwound');
  assert.deepEqual(domNames(instance, source.root), ['nameInput', 'details', 'leaf']);
  assert.equal(instance.nodes().has(note), true);
  assert.equal(editor.canUndo.value, false, 'a refused compound is not history');
  instance.dispose();
});

edit('applyAll: a throw out of a member\'s execution unwinds the whole compound, no history', () => {
  // a component ctor throw is contained per node (spec.ts _build renders a placeholder), so an
  // escaping throw is an engine invariant failure mid-rerender — modeled by a rerender that
  // throws once; the throwing member's mutation has landed and must be unwound with the rest
  const {source, instance, editor} = open();
  const input = find(source, 'nameInput');
  const leaf = find(source, 'leaf');
  const before = instance.dump();
  const rerender = instance.rerender.bind(instance);
  let calls = 0;
  instance.rerender = (node) => {
    if (++calls === 2)
      throw new Error('boom');
    rerender(node);
  };
  try {
    assert.throws(() => editor.applyAll([{op: 'remove', node: input}, {op: 'remove', node: leaf}]),
      /boom/, 'the error propagates');
  } finally {
    instance.rerender = rerender;
  }
  assert.deepEqual(instance.dump(), before, 'the document deep-equals the pre-call dump');
  assert.deepEqual(domNames(instance, source.root), ['nameInput', 'details', 'leaf'],
    'and the DOM was re-rendered back');
  assert.equal(editor.canUndo.value, false, 'a failed compound is not history');
  instance.dispose();
});

edit('applyAll publishes once, carrying every member and the structural OR', () => {
  const {source, instance, editor} = open();
  const input = find(source, 'nameInput');
  const leaf = find(source, 'leaf');
  const seen = [];
  const scope = new Scope();
  scope.effect(() => {
    const applied = editor.onDidApply.value;
    if (applied)
      seen.push(applied);
  });

  editor.applyAll([
    {op: 'set-prop', node: input, name: 'label', value: 'Edited'},
    {op: 'remove', node: leaf},
  ]);
  assert.equal(seen.length, 1, 'one publication per applyAll');
  assert.deepEqual(seen[0].patches.map((p) => p.op), ['set-prop', 'remove']);
  assert.equal(seen[0].structural, true, 'any structural member makes the change structural');
  assert.equal(seen[0].origin, 'apply');

  editor.applyAll([
    {op: 'set-prop', node: input, name: 'label', value: 'Twice'},
    {op: 'set-prop', node: input, name: 'size', value: 2},
  ]);
  assert.equal(seen.length, 2);
  assert.equal(seen[1].structural, false);

  editor.undo();
  assert.equal(seen.length, 3);
  assert.equal(seen[2].origin, 'undo');
  assert.deepEqual(seen[2].patches.map((p) => [p.op, p.name]),
    [['set-prop', 'size'], ['set-prop', 'label']], 'undo publishes the inverses it executed, in order');
  scope.dispose();
  instance.dispose();
});

edit('coalescing never merges into or out of a compound', () => {
  const {source, instance, editor} = open();
  const input = find(source, 'nameInput');

  editor.apply({op: 'set-prop', node: input, name: 'label', value: 'A'});
  editor.applyAll([
    {op: 'set-prop', node: input, name: 'label', value: 'B'},
    {op: 'set-prop', node: input, name: 'size', value: 1},
  ]);
  editor.apply({op: 'set-prop', node: input, name: 'label', value: 'C'});

  editor.undo();
  assert.equal(input.props.label, 'B', 'the patch after the compound was its own entry');
  editor.undo();
  assert.equal(input.props.label, 'A', 'and the compound never absorbed it');
  assert.equal(input.props.size, undefined);
  editor.undo();
  assert.equal(input.props.label, 'Name', 'nor did the compound merge into the entry before it');
  assert.equal(editor.canUndo.value, false);
  instance.dispose();
});

edit('patches inside a never-activated tab pane show fresh content when the pane first opens', async () => {
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'u2-lazy-tabs', name: 'tabs', children: [
      {tag: 'u2-fake-box', name: 'pane1', props: {title: 'First'}},
      {tag: 'u2-fake-box', name: 'pane2', props: {title: 'Second'}, children: [
        {tag: 'u2-fake-box', name: 'inner', children: [
          {tag: 'span', name: 'note', props: {text: 'old'}},
        ]},
      ]},
    ]},
  };
  const {instance, editor} = open(source);
  const note = find(source, 'note');
  const before = element(instance.nodes().get(note));
  assert.equal(before.isConnected, false, 'pane 2 was never activated — its subtree is detached');

  // the deep re-render lands inside the detached subtree; the pane-level one promotes to the tabs
  editor.apply({op: 'set-prop', node: note, name: 'text', value: 'new'});
  editor.apply({op: 'add', parent: find(source, 'inner'), index: 1,
    node: {tag: 'span', name: 'added', props: {text: 'added'}}});

  const strip = instance.nodes().get(source.root);
  strip.activeTab.value = 'tab-1';
  await flush();
  const panel = strip.root.querySelectorAll('.u2-tabs-panel')[1];
  const after = element(instance.nodes().get(find(source, 'note')));
  assert.equal(panel.contains(after), true, 'the first activation shows the re-rendered content');
  assert.equal(after.textContent, 'new');
  assert.equal(panel.contains(before), false, 'and never the disposed old element');
  assert.equal(panel.contains(element(instance.nodes().get(find(source, 'added')))), true,
    'a node added while the pane was hidden is there too');
  instance.dispose();
});

/* WO-11 — the binding substrate at the editor's boundary: path validation, sub-bind keys,
   cycle refusal, structured event entries and what rename does to them. */

edit('canApply(set-bind): the path must parse, and a dotted key needs a subBindable head', () => {
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-box', name: 'layout', children: [
      {tag: 'u2-fake-src', name: 'orders', props: {func: 'demoOrders'}},
      {tag: 'u2-fake-input', name: 'daysInput', props: {label: 'Days'}},
    ]},
  };
  const {instance, editor} = open(source);
  const orders = find(source, 'orders');
  const input = find(source, 'daysInput');

  assert.equal(editor.canApply({op: 'set-bind', node: orders, name: 'params.days',
    path: '$.daysInput'}), null, 'a dotted key on a subBindable prop');
  assert.match(editor.canApply({op: 'set-bind', node: orders, name: 'func.x', path: '$.daysInput'}),
    /prop "func" does not take sub-binds/);
  assert.match(editor.canApply({op: 'set-bind', node: orders, name: 'nope.x', path: '$.daysInput'}),
    /has no prop "nope"/);
  assert.match(editor.canApply({op: 'set-bind', node: input, name: 'value', path: 'daysInput'}),
    /must start with "\$\."/);
  for (const path of ['$.a b', '$.1x', '$.a[\'b', '$.a..b']) {
    assert.match(editor.canApply({op: 'set-bind', node: input, name: 'value', path}),
      /is not a valid bind path/, path);
  }
  assert.equal(editor.canApply({op: 'set-bind', node: input, name: 'value',
    path: '$.orders.currentRow[\'Mol Weight\']'}), null,
  'an unresolvable path is contained, not refused');

  editor.apply({op: 'set-bind', node: orders, name: 'params.days', path: '$.daysInput'});
  assert.deepEqual(instance.dump().root.children[0].bind, {'params.days': '$.daysInput'});
  editor.undo();
  assert.equal('bind' in instance.dump().root.children[0], false);
  instance.dispose();
});

edit('canApply(set-bind): a patch that would close a loop is refused, naming it', () => {
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-box', name: 'layout', children: [
      {tag: 'u2-fake-src', name: 'orders', props: {func: 'demoOrders'}},
      {tag: 'u2-fake-input', name: 'daysInput', bind: {value: '$.orders'}},
      {tag: 'u2-fake-input', name: 'other', props: {label: 'Other'}},
    ]},
  };
  const {instance, editor} = open(source);
  const orders = find(source, 'orders');
  const before = instance.dump();

  assert.equal(editor.canApply({op: 'set-bind', node: orders, name: 'params.days', path: '$.daysInput'}),
    'binding "$.daysInput" would create a cycle: orders → daysInput → orders');
  assert.throws(() => editor.apply({op: 'set-bind', node: orders, name: 'params.days',
    path: '$.daysInput'}), /would create a cycle/);
  assert.deepEqual(instance.dump(), before, 'the probe left the document exactly as it was');
  assert.equal(editor.canUndo.value, false);

  assert.equal(editor.canApply({op: 'set-bind', node: orders, name: 'params.days', path: '$.other'}),
    null, 'the same key to a node off the loop is fine');
  assert.equal(editor.canApply({op: 'set-bind', node: orders, name: 'params.days', path: '$.orders'}),
    'binding "$.orders" would create a cycle: orders → orders', 'a self-loop too');
  instance.dispose();
});

edit('canApply(set-event): both forms validated, args must be a plain record', () => {
  const {source, instance, editor} = open();
  const input = find(source, 'nameInput');
  const check = (command) => editor.canApply({op: 'set-event', node: input, event: 'input', command});

  assert.equal(check('cmd:save'), null);
  assert.equal(check({cmd: 'cmd:Pkg:Save'}), null);
  assert.equal(check({cmd: 'cmd:Pkg:Save', args: {id: '$.x', n: 1}}), null);
  assert.equal(check(undefined), null, 'clearing the wiring');
  assert.match(check({cmd: 'save'}), /must name a command as 'cmd:' followed by the command name/);
  assert.match(check({args: {}}), /got '\{"args":\{\}\}'/, 'a missing cmd reports the entry');
  assert.equal(check({cmd: 'cmd:Pkg:Save', args: ['a']}),
    'event args must be a plain object of argument values');
  assert.equal(check({cmd: 'cmd:Pkg:Save', args: 'x'}),
    'event args must be a plain object of argument values');

  editor.apply({op: 'set-event', node: input, event: 'input',
    command: {cmd: 'cmd:Pkg:Save', args: {n: 1}}});
  assert.deepEqual(instance.dump().root.children[0].on.input, {cmd: 'cmd:Pkg:Save', args: {n: 1}});
  editor.undo();
  assert.equal('on' in instance.dump().root.children[0], false);
  instance.dispose();
});

edit('rename: structured entries move too — the component tier cmd and $.old arg paths', () => {
  const source = {
    $schema: 'dg-ui/1',
    root: {tag: 'u2-fake-box', name: 'layout', children: [
      {tag: 'u2-fake-src', name: 'old', props: {func: 'f'}},
      {tag: 'u2-fake-input', name: 'consumer', bind: {value: '$.old'},
        on: {input: {cmd: 'cmd:old.refresh', args: {row: '$.old.currentRow', keep: '$.older',
          note: 'plain', n: 3}},
        change: {cmd: 'cmd:Pkg:Save', args: {id: '$.old'}}}},
    ]},
  };
  const {instance, editor} = open(source, {older: 'x'});
  const before = instance.dump();

  editor.apply({op: 'rename', node: find(source, 'old'), name: 'orders'});
  const on = instance.dump().root.children[1].on;
  assert.deepEqual(on.input, {cmd: 'cmd:orders.refresh',
    args: {row: '$.orders.currentRow', keep: '$.older', note: 'plain', n: 3}});
  assert.deepEqual(on.change, {cmd: 'cmd:Pkg:Save', args: {id: '$.orders'}},
    'a platform-tier cmd keeps its name; only its args move');
  assert.equal(instance.dump().root.children[1].bind.value, '$.orders');

  editor.undo();
  assert.deepEqual(instance.dump(), before, 'the reverse rename rewrites everything back');
  instance.dispose();
});

edit('renameRefs: structured entries, walked directly', () => {
  const node = {tag: 'div',
    on: {a: {cmd: 'cmd:old.run', args: {x: '$.old', y: '$.oldish', z: 7}},
      b: 'cmd:old.run', c: {cmd: 'cmd:old'}, d: {cmd: 'cmd:x:y'}}};
  assert.deepEqual(renameRefs(node, 'old', 'next'), [node]);
  assert.deepEqual(node.on, {
    a: {cmd: 'cmd:next.run', args: {x: '$.next', y: '$.oldish', z: 7}},
    b: 'cmd:next.run',
    c: {cmd: 'cmd:old'},
    d: {cmd: 'cmd:x:y'},
  });
  assert.deepEqual(renameRefs({tag: 'div', on: {a: {cmd: 'cmd:other.run'}}}, 'old', 'next'), []);
});
