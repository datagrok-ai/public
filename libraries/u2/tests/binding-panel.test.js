/* Phase 3.4 (UB-7, R10–R12) — the panel over universal binding: unboundOf/bindsOf/read() as the
   model (R10), the chip in place of a bound row's editor (R11), the picker button on every
   property row, the unified Bindings section, one patch per gesture through reissue, a refresh
   that echoes nothing — and the viewer tier path unmoved (UB-8). Harness idioms of
   tests/handler-edit.test.js; dg-stub.mjs serves the platform. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {Control} from '../src/core/component.js';
import {TextInput} from '../src/components/inputs/text-input.js';
import {Registry} from '../src/spec/registry.js';
import {registerAll} from '../src/spec/registrations.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {SpecEditor} from '../src/spec/editor.js';
import {dfBindings} from '../src/sources/df-bindings.js';
import {DataFrame, Property, WidgetDescriptor, platform} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {SpecNodeRef} = await import('../src/dg/designer/node-ref.js');
const {SpecNodeHandler, propsFor, disposePanel} = await import('../src/dg/designer/handler.js');
const {bindsOf, unboundOf} = await import('../src/dg/designer/prop-model.js');
const {registerPlatformComponents} = await import('../src/dg/viewers/registrations.js');
const {shell} = await import('datagrok-api/grok');

/** A parent whose childProps COLLIDE with the child's own declaration — the dedup case. */
const DECK = {
  tag: 'u2-p-deck',
  create: () => new Control(),
  description: 'Fake stack that reads a title off every child',
  props: [],
  childProps: [
    {name: 'title', type: 'string', description: 'Card title'},
    {name: 'badge', type: 'string', description: 'Corner badge'},
  ],
  acceptsChildren: true,
  example: {tag: 'u2-p-deck'},
};

const CARD = {
  tag: 'u2-p-card',
  create: (props) => new TextInput({label: typeof props.title === 'string' ? props.title : undefined}),
  description: 'Fake card that declares the title its parent also reads',
  props: [{name: 'title', type: 'string'}],
  example: {tag: 'u2-p-card'},
};

function editable() {
  return {$schema: 'dg-ui/1', root: {tag: 'u2-div-v', name: 'layout', children: [
    {tag: 'u2-text-input', name: 'plain', props: {label: 'Name', value: 'Aspirin'}},
    {tag: 'u2-text-input', name: 'boundValue', props: {label: 'Bound'}, bind: {value: '$.reagent'}},
    {tag: 'u2-button', name: 'save', props: {icon: 'save'}, bind: {text: '$.reagent'}},
    {tag: 'u2-tabs', name: 'tabs', children: [
      {tag: 'u2-panel', name: 'firstPane', props: {title: 'First'}},
    ]},
    {tag: 'u2-p-deck', name: 'deck', children: [
      {tag: 'u2-p-card', name: 'card', props: {title: 'Card'}},
    ]},
    {tag: 'div', name: 'block', props: {text: 'Hello'}},
  ]}};
}

function find(spec, name) {
  const walk = (n) => {
    if (n.name === name)
      return n;
    for (const child of n.children ?? []) {
      const found = walk(child);
      if (found)
        return found;
    }
    return null;
  };
  return [spec.root, ...spec.components ?? []].map(walk).find((node) => node !== null) ?? null;
}

function section(panel, title) {
  const kids = [...panel.children];
  const i = kids.findIndex((el) => el.tagName === 'H3' && el.textContent === title);
  if (i >= 0)
    return kids[i + 1];
  const pane = [...panel.querySelectorAll('.u2-accordion-pane')]
    .find((p) => p.querySelector('.u2-accordion-title')?.textContent === title);
  return pane?.querySelector('.u2-accordion-content') ?? null;
}

function field(panel, title, name) {
  return section(panel, title).querySelector(`[data-u2-prop="${name}"]`)
    .querySelector('input, textarea, select');
}

const bindRows = (panel) => [...section(panel, 'Bindings').querySelectorAll('[data-u2-prop]')]
  .map((el) => el.dataset.u2Prop);

function picker() {
  const dialog = document.querySelector('.u2-bind-picker');
  const list = dialog?.querySelector('.u2-list');
  if (list) {
    list.clientHeight = 400;
    fire(list, 'scroll');
  }
  return dialog;
}

const pickerRow = (dialog, label) => [...dialog.querySelectorAll('.u2-list-row')]
  .find((row) => (row.querySelector('.u2-tree-label')?.textContent ?? '').startsWith(label));

const okButton = (dialog) => [...dialog.querySelectorAll('.u2-dialog-buttons button')]
  .find((b) => b.textContent === 'OK');

const cancelButton = (dialog) => [...dialog.querySelectorAll('.u2-dialog-buttons button')]
  .find((b) => b.textContent === 'CANCEL');

function edit(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const warn = console.warn;
    const warnings = [];
    const infos = [];
    const warning = shell.warning;
    const info = shell.info;
    console.warn = () => {};
    shell.warning = (message) => warnings.push(message);
    shell.info = (message) => infos.push(message);
    const reg = new Registry();
    registerAll(reg);
    for (const meta of [DECK, CARD])
      reg.register(meta);
    const spec = editable();
    const instance = renderSpec(spec, new SpecContext({data: {reagent: 'Ethanol', tone: 'red'}}), reg);
    const editor = new SpecEditor(instance);
    const patches = [];
    const apply = editor.apply.bind(editor);
    editor.apply = (patch) => {
      patches.push(patch);
      apply(patch);
    };
    const handler = new SpecNodeHandler();
    try {
      await body({
        instance, editor, patches, warnings, infos,
        node: (n) => find(spec, n),
        ref: (n) => new SpecNodeRef(instance, find(spec, n), editor),
        panel: (n) => handler.renderProperties(new SpecNodeRef(instance, find(spec, n), editor)),
        roPanel: (n) => handler.renderProperties(new SpecNodeRef(instance, find(spec, n))),
      });
    } finally {
      disposePanel();
      instance.dispose();
      shell.warning = warning;
      shell.info = info;
      console.warn = warn;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

edit('unboundOf: every declared-but-unbound prop — component, HTML tag, minus what is bound',
  ({ref}) => {
    const plain = unboundOf(ref('plain'));
    for (const name of ['value', 'label', 'placeholder', 'color'])
      assert.ok(plain.includes(name), `${name} is declared and unbound`);
    const bound = unboundOf(ref('boundValue'));
    assert.equal(bound.includes('value'), false, 'a live-tier bind excludes its prop');
    assert.ok(bound.includes('label'));
    assert.equal(unboundOf(ref('save')).includes('text'), false, 'a re-render bind excludes its prop too');
    const block = unboundOf(ref('block'));
    for (const name of ['text', 'cls', 'color'])
      assert.ok(block.includes(name), `${name}: HTML props and the appearance group both offer`);
  });

edit('unboundOf: childProps join, deduped against the child\'s own declaration', ({ref}) => {
  const pane = unboundOf(ref('firstPane'));
  for (const name of ['title', 'icon', 'cls'])
    assert.ok(pane.includes(name), `${name} is offered on the pane`);
  const card = unboundOf(ref('card'));
  assert.equal(card.filter((name) => name === 'title').length, 1,
    'declared by the card AND read by the deck — listed once');
  assert.ok(card.includes('badge'), 'the parent-only prop still joins');
});

edit('a bound prop\'s read(): the bind path for the chip, the live value everywhere else (S2)',
  ({ref}) => {
    const props = propsFor(ref('boundValue'), true).find((s) => s.title === 'Properties');
    assert.equal(props.values.value, '$.reagent');
    assert.equal(props.read().value, '$.reagent', 'what the chip shows and refreshes to');
    assert.equal(props.props.find((p) => p.name === 'value').set, undefined,
      'the document-side truth: a bound prop has no set');
    const plain = propsFor(ref('boundValue')).find((s) => s.title === 'Properties');
    assert.equal(plain.values.value, 'Ethanol',
      'without the flag a bound prop reads its live value — the read-only tables\' truth');
    assert.deepEqual(bindsOf(ref('boundValue')), {value: '$.reagent'}, 'bindsOf is the raw document');
  });

edit('the read-only panel shows live values for bound props; paths stay in the Bindings table (S2)',
  ({roPanel}) => {
    const shown = roPanel('boundValue');
    const props = section(shown, 'Properties');
    assert.match(props.textContent, /Ethanol/, 'the live value, not the path');
    assert.doesNotMatch(props.textContent, /\$\.reagent/, 'the path never duplicates into the props table');
    assert.match(section(shown, 'Bindings').textContent, /\$\.reagent/, 'the Bindings table carries the path');
  });

edit('every property row of a plain node carries its picker button — childProps rows included',
  ({ref, panel}) => {
    const shown = panel('plain');
    for (const s of propsFor(ref('plain'))) {
      for (const p of s.props) {
        assert.notEqual(section(shown, s.title)
          .querySelector(`[data-u2-prop="${p.name}"] [data-u2-bind-pick="${p.name}"]`), null,
        `${s.title}/${p.name} has its … button (UB-7a)`);
      }
    }
    const pane = panel('firstPane');
    for (const name of ['title', 'icon']) {
      assert.notEqual(section(pane, 'Parent (u2-tabs)')
        .querySelector(`[data-u2-prop="${name}"] [data-u2-bind-pick="${name}"]`), null,
      `the parent's ${name} row has its … button`);
    }
  });

edit('a bound row renders the chip: the caption, the path, ⇄ only where the declaration says two-way',
  ({panel}) => {
    const editorCaption = section(panel('plain'), 'Properties')
      .querySelector('[data-u2-prop="value"] .u2-input-label').textContent;
    const bound = panel('boundValue');
    const row = section(bound, 'Properties').querySelector('[data-u2-prop="value"]');
    const chip = row.querySelector('.u2-bind-chip');
    assert.notEqual(chip, null, 'the chip replaces the editor (UB-7b)');
    assert.equal(row.querySelector('.u2-input-label').textContent, editorCaption,
      'the chip row keeps the caption an editor row would carry (F-4)');
    assert.equal(row.getAttribute('title'), '$.reagent', 'the full path rides the row tooltip (P1)');
    assert.equal(chip.querySelector('.u2-bind-chip-path').textContent, '$.reagent');
    assert.notEqual(chip.querySelector('.u2-bind-chip-two-way'), null, 'value is twoWay — ⇄ shows');
    assert.notEqual(row.querySelector('[data-u2-unbind="value"]'), null, 'the × carries its handle');
    assert.equal(row.querySelector('input, textarea, select'), null, 'the editor is gone');

    const save = panel('save');
    const text = section(save, 'Properties').querySelector('[data-u2-prop="text"] .u2-bind-chip');
    assert.notEqual(text, null, 'a re-render bind gets the same chip');
    assert.equal(text.querySelector('.u2-bind-chip-two-way'), null, 'strictly one-way — no ⇄ (UB-3)');
  });

edit('× unbinds with exactly one patch; the SAME panel redraws the row as an editor (F-3)',
  async ({node, panel, patches, infos}) => {
    const shown = panel('boundValue');
    const writes = shell.dart.writes.length;
    fire(section(shown, 'Properties').querySelector('[data-u2-unbind="value"]'), 'click');
    await flush();

    assert.deepEqual(patches.map((p) => [p.op, p.name, p.path]), [['set-bind', 'value', undefined]]);
    assert.equal(node('boundValue').bind, undefined);
    assert.deepEqual(infos, ['value: binding removed'], 'clearing wiring is deliberate but never silent');
    assert.equal(shell.dart.writes.length - writes, 0,
      'no shell.o write — a same-identity re-issue is a platform no-op, the panel redraws itself');

    const row = section(shown, 'Properties').querySelector('[data-u2-prop="value"]');
    assert.equal(row.querySelector('.u2-bind-chip'), null);
    assert.notEqual(row.querySelector('input'), null, 'the editor is back without re-selecting the node');
    assert.deepEqual(bindRows(shown), ['add-binding'], 'and the bound row left the Bindings section');
  });

edit('a row-button pick is ONE patch through the funnel; the SAME panel flips the row to the chip (F-1)',
  async ({node, panel, patches}) => {
    const shown = panel('plain');
    const writes = shell.dart.writes.length;
    fire(section(shown, 'Properties').querySelector('[data-u2-bind-pick="label"]'), 'click');
    const dialog = picker();
    assert.notEqual(dialog, null, 'the row button opens the picker');
    fire(pickerRow(dialog, 'reagent'), 'click');
    okButton(dialog).click();
    await flush();

    assert.deepEqual(patches.map((p) => [p.op, p.name, p.path]), [['set-bind', 'label', '$.reagent']]);
    assert.equal(node('plain').bind.label, '$.reagent');
    assert.equal(shell.dart.writes.length - writes, 0,
      'no shell.o write — a same-identity re-issue is a platform no-op, the panel redraws itself');

    const row = section(shown, 'Properties').querySelector('[data-u2-prop="label"]');
    assert.notEqual(row.querySelector('.u2-bind-chip'), null,
      'the freshly bound row renders the chip without re-selecting the node');
    assert.notEqual(row.querySelector('.u2-input-label'), null, 'with its caption (F-4)');
    assert.equal(field(shown, 'Bindings', 'label').value, '$.reagent', 'and the Bindings row is there');
  });

edit('a bound prop\'s row cannot commit a dead literal — the stale-row path is closed (F-2)',
  async ({editor, node, panel, patches}) => {
    const shown = panel('plain');
    const stale = section(shown, 'Properties').querySelector('[data-u2-prop="label"] input');
    editor.apply({op: 'set-bind', node: node('plain'), name: 'label', path: '$.reagent'});
    await flush();
    assert.notEqual(section(shown, 'Properties').querySelector('[data-u2-prop="label"] .u2-bind-chip'),
      null, 'a set-bind from any surface flips the live panel\'s row');

    stale.value = 'Stale literal';
    fire(stale, 'input');
    fire(stale, 'change');
    await flush();
    assert.deepEqual(patches.filter((p) => p.op === 'set-prop'), [],
      'the replaced row commits nothing over the bound prop');
  });

edit('Bindings: bound rows only, and "Add binding…" humanizes every unbound prop', ({panel}) => {
  const bound = panel('boundValue');
  assert.deepEqual(bindRows(bound), ['value', 'add-binding']);
  assert.equal(section(bound, 'Bindings').querySelector('[data-u2-prop="value"] .u2-input-label')
    .textContent, 'Value', 'a bound row\'s label is humanized like the Add-binding list (P5)');
  const select = section(bound, 'Bindings').querySelector('[data-u2-prop="add-binding"] select');
  const options = [...select.querySelectorAll('option')].map((o) => [o.value, o.textContent]);
  assert.deepEqual(options.find(([value]) => value === 'tooltipText'), ['tooltipText', 'Tooltip Text']);
  assert.equal(options.some(([value]) => value === 'value'), false, 'what is bound is not offered');

  assert.deepEqual(bindRows(panel('block')), ['add-binding'], 'the same shape for an HTML node');
});

edit('a chip-body click opens the picker; the × does not — unbind is one gesture, one patch',
  async ({node, panel, patches}) => {
    const shown = panel('save');
    fire(section(shown, 'Properties').querySelector('[data-u2-prop="text"] .u2-bind-chip'), 'click');
    const dialog = picker();
    assert.notEqual(dialog, null, 'the chip re-picks (R11)');
    cancelButton(dialog).click();
    await flush();
    assert.equal(patches.length, 0, 'a cancelled re-pick patches nothing');
    assert.equal(node('save').bind.text, '$.reagent');

    const after = panel('save');
    fire(section(after, 'Properties').querySelector('[data-u2-unbind="text"]'), 'click');
    await flush();
    assert.equal(document.querySelector('.u2-bind-picker'), null, 'the × never opens the picker');
    assert.deepEqual(patches.map((p) => [p.op, p.name, p.path]), [['set-bind', 'text', undefined]]);
    assert.equal(node('save').bind, undefined);
  });

edit('the fold headers stay live after a bind, an unbind, and on a fresh panel post-set-bind (B1)',
  async ({editor, node, panel}) => {
    const toggles = (host, title) => {
      const header = [...host.querySelectorAll('.u2-accordion-pane')]
        .find((p) => p.querySelector('.u2-accordion-title')?.textContent.startsWith(title))
        ?.querySelector('.u2-accordion-header');
      if (!header)
        return false;
      const before = header.getAttribute('aria-expanded');
      fire(header, 'click');
      const flipped = header.getAttribute('aria-expanded') !== before;
      fire(header, 'click');
      return flipped && header.getAttribute('aria-expanded') === before;
    };
    const shown = panel('plain');
    fire(section(shown, 'Properties').querySelector('[data-u2-bind-pick="label"]'), 'click');
    const dialog = picker();
    fire(pickerRow(dialog, 'reagent'), 'click');
    okButton(dialog).click();
    await flush();
    assert.ok(toggles(shown, 'Bindings'), 'the Bindings fold still toggles after a panel bind');
    assert.ok(toggles(shown, 'Appearance'), 'the Appearance fold too');

    fire(section(shown, 'Properties').querySelector('[data-u2-unbind="label"]'), 'click');
    await flush();
    assert.ok(toggles(shown, 'Bindings'), 'and again after the unbind');

    editor.apply({op: 'set-bind', node: node('plain'), name: 'label', path: '$.reagent'});
    await flush();
    const fresh = panel('plain');
    assert.ok(toggles(fresh, 'Bindings'),
      'a fresh panel built while the last patch was a set-bind is live from birth');
  });

edit('a re-pick names the prop in the title and arrives with the current path\'s row selected (M3)',
  async ({panel}) => {
    const shown = panel('save');
    fire(section(shown, 'Properties').querySelector('[data-u2-prop="text"] .u2-bind-chip'), 'click');
    let dialog = picker();
    assert.equal(dialog.querySelector('.u2-dialog-title-text').textContent, 'text — bind to');
    await flush();
    dialog = picker();
    const selected = dialog.querySelector('.u2-list-row-selected .u2-tree-label');
    assert.ok(selected?.textContent.startsWith('reagent'),
      `the current path's row is selected: ${selected?.textContent}`);
    cancelButton(dialog).click();
    await flush();

    fire(section(shown, 'Properties').querySelector('[data-u2-bind-pick="icon"]'), 'click');
    const unbound = picker();
    assert.equal(unbound.querySelector('.u2-dialog-title-text').textContent, 'icon — bind to');
    assert.equal(unbound.querySelector('.u2-list-row-selected'), null,
      'nothing preselected where nothing is bound');
    cancelButton(unbound).click();
    await flush();
  });

edit('typing a new path into a bound Bindings row re-points the binding — one patch, node follows',
  async ({node, panel, patches, instance}) => {
    const shown = panel('save');
    const row = field(shown, 'Bindings', 'text');
    row.value = '$.tone';
    fire(row, 'change');
    fire(row, 'blur');
    await flush();
    assert.deepEqual(patches.map((p) => [p.op, p.name, p.path]), [['set-bind', 'text', '$.tone']]);
    assert.equal(node('save').bind.text, '$.tone');
    assert.equal(instance.node('save').root.textContent, 'red',
      'the node re-rendered onto the new source\'s value');
  });

edit('a panel refresh echoes nothing — a chip is display, not an editor',
  async ({editor, node, panel, patches}) => {
    panel('boundValue');
    editor.apply({op: 'set-prop', node: node('boundValue'), name: 'label', value: 'Renamed'});
    await flush();
    assert.deepEqual(patches.map((p) => [p.op, p.name]), [['set-prop', 'label']],
      'only the patch we applied — the refresh committed no set-bind');
  });

edit('undo after a pick takes the chip away with the binding; redo brings both back',
  async ({editor, node, panel}) => {
    const shown = panel('plain');
    fire(section(shown, 'Properties').querySelector('[data-u2-bind-pick="label"]'), 'click');
    const dialog = picker();
    fire(pickerRow(dialog, 'reagent'), 'click');
    okButton(dialog).click();
    await flush();
    assert.equal(node('plain').bind.label, '$.reagent');

    editor.undo();
    assert.equal(node('plain').bind, undefined);
    const undone = panel('plain');
    assert.equal(section(undone, 'Properties').querySelector('[data-u2-prop="label"] .u2-bind-chip'),
      null, 'the rendering matches the document at each step');
    assert.deepEqual(bindRows(undone), ['add-binding']);

    editor.redo();
    assert.equal(node('plain').bind.label, '$.reagent');
    const redone = panel('plain');
    assert.notEqual(section(redone, 'Properties').querySelector('[data-u2-prop="label"] .u2-bind-chip'),
      null);
  });

test('the viewer tier is unmoved (UB-8): no row buttons, no chips — Bindings alone carries the gesture',
  async () => {
    const live = Scope.liveCount;
    const warn = console.warn;
    console.warn = () => {};
    WidgetDescriptor.registry = [new WidgetDescriptor('Scatter plot', [
      new Property('xColumnName', 'string', {category: 'Data'}),
    ])];
    const reg = new Registry();
    registerPlatformComponents(reg);
    const scope = new Scope();
    const df = new DataFrame([{name: 'city', type: 'string'}], [{city: 'Kyiv'}], 'orders');
    const spec = {$schema: 'dg-ui/1', root: {tag: 'u2-div-v', name: 'box', children: [
      {tag: 'u2-viewer-scatter-plot', name: 'plot', bind: {table: '$.orders'}, props: {xColumnName: 'city'}},
    ]}};
    const instance = renderSpec(spec,
      new SpecContext({data: {orders: dfBindings(signal(df), scope)}}), reg, {designTime: true});
    const editor = new SpecEditor(instance);
    const handler = new SpecNodeHandler();
    try {
      const plot = find(spec, 'plot');
      assert.equal(instance.nodes().get(plot).propertyTier, true);
      const shown = handler.renderProperties(new SpecNodeRef(instance, plot, editor));
      assert.equal(shown.querySelector('.u2-bind-chip'), null, 'no chips on the tier path');
      assert.equal(shown.querySelector('.u2-designer-look [data-u2-bind-pick]'), null,
        'the platform grid carries no row buttons');
      assert.equal(shown.querySelector('[data-u2-prop="xColumnName"]'), null, 'no u2 row for a look prop');
      assert.deepEqual([...shown.querySelectorAll('[data-u2-bind-pick]')].map((el) => el.dataset.u2BindPick),
        ['table', 'add-binding', 'd4-data-event'],
      'the Bindings rows\' buttons and the Events picker only');
    } finally {
      disposePanel();
      instance.dispose();
      scope.dispose();
      WidgetDescriptor.registry = [];
      platform.reset();
      console.warn = warn;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });

test('Phase 1.7 acceptance: canApply(set-bind) answers null on u2-button.text and h3.text', async () => {
  const live = Scope.liveCount;
  const reg = new Registry();
  registerAll(reg);
  const spec = {$schema: 'dg-ui/1', root: {tag: 'u2-div-v', name: 'box', children: [
    {tag: 'u2-button', name: 'go', props: {text: 'Run'}},
    {tag: 'h3', name: 'title', props: {text: 'Header'}},
  ]}};
  const instance = renderSpec(spec, new SpecContext({data: {s: 'x'}}), reg);
  const editor = new SpecEditor(instance);
  try {
    assert.equal(editor.canApply({op: 'set-bind', node: find(spec, 'go'), name: 'text', path: '$.s'}),
      null, 'a re-render bind on a component prop may apply');
    assert.equal(editor.canApply({op: 'set-bind', node: find(spec, 'title'), name: 'text', path: '$.s'}),
      null, 'and on a plain HTML node');
  } finally {
    instance.dispose();
    resetDom();
    await flush();
  }
  assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
});

test('a subBindable prop is never offered whole: params is off the Add-binding list', async () => {
  const live = Scope.liveCount;
  const warn = console.warn;
  console.warn = () => {};
  const reg = new Registry();
  registerAll(reg);
  const spec = {$schema: 'dg-ui/1',
    components: [{tag: 'u2-func-source', name: 'orders', props: {func: 'Sales:Recent'}}],
    root: {tag: 'div'}};
  const instance = renderSpec(spec, new SpecContext(), reg);
  const editor = new SpecEditor(instance);
  try {
    const offered = unboundOf(new SpecNodeRef(instance, spec.components[0], editor));
    assert.equal(offered.includes('params'), false,
      'bound through its dotted sub-bind rows only, never as a whole');
    for (const name of ['func', 'auto', 'debounce', 'designData', 'sample'])
      assert.ok(offered.includes(name), `${name} is still offered`);
  } finally {
    instance.dispose();
    console.warn = warn;
    resetDom();
    await flush();
  }
  assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
});
