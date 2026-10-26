import {signal, computed, Scope, Control, Filters, FilterBuilder} from '../../src/index.js';
import {divH, divV, span, button} from '../../src/core/elements.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-inputs-css', '../../css/inputs.css');
injectOnce('u2-choice-css', '../../css/choice.css');
injectOnce('u2-number-css', '../../css/number.css');
injectOnce('u2-icons-css', '../../css/icons.css');
injectOnce('u2-buttons-css', '../../css/buttons.css');
injectOnce('u2-tooltip-css', '../../css/tooltip.css');
injectOnce('u2-badge-css', '../../css/badge.css');
injectOnce('u2-tags-css', '../../css/tags.css');
injectOnce('u2-typeahead-css', '../../css/typeahead.css');
injectOnce('u2-date-css', '../../css/date.css');
injectOnce('u2-filter-css', '../../css/filter.css');

/** A literal schema — the platform-free proof: plain IProperty objects plus value lists. */
const SCHEMA = Filters.schema([
  {name: 'name', type: 'string', friendlyName: 'Name'},
  {name: 'age', type: 'int', friendlyName: 'Age', min: 0, max: 120},
  {name: 'mw', type: 'double', friendlyName: 'Mol weight'},
  {name: 'sex', type: 'string', friendlyName: 'Sex', choices: ['F', 'M']},
  {name: 'status', type: 'string', friendlyName: 'Status'},
  {name: 'active', type: 'bool', friendlyName: 'Active'},
  {name: 'created', type: 'datetime', friendlyName: 'Created'},
  {name: 'owner', type: 'string', friendlyName: 'Owner', ref: 'Core.users'},
], {
  name: ['Aspirin', 'Ibuprofen', 'Paracetamol', 'Naproxen', 'Diclofenac'],
  status: ['Open', 'Blocked', 'Closed', 'Draft'],
  owner: ['Alice', 'Bob', 'Carol', 'Dave'],
});

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

export async function render(main) {
  main.append(el('h1', null, 'Filter builder'));
  const intro = el('p');
  intro.innerHTML = 'Schema-driven query builder over an immutable <code>FilterGroup</code>: one ' +
    'row per condition (property → operator → value), a single <code>and</code>/<code>or</code> ' +
    'toggle, <code>+</code> adds and <code>−</code> removes. Editors come from the property kind ' +
    '(number, choice, a date-time box that also takes a span such as <code>-1w</code>, tags for <code>in</code>), the query ' +
    'string is <code>Filters.format(value)</code> and problems come from ' +
    '<code>Filters.validate</code>. Everything here runs on a literal schema — no platform.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  const parts = [];
  const section = (title, builder) => {
    main.append(el('h2', null, title));
    const component = Control.build(builder);
    parts.push(component);
    main.append(component.root);
    return component;
  };

  const tree = signal(Filters.group('and', [
    Filters.cond('age', '>', 30),
    Filters.cond('sex', '=', 'F'),
  ]));
  section('Vertical (sidebar), live query readout', () => {
    const fb = new FilterBuilder({label: 'Criteria', schema: SCHEMA, bind: tree, showQuery: true,
      name: 'criteria'});
    return [fb, readout('problems', computed(() => fb.problems.value.length === 0 ? 'none' :
      fb.problems.value.map((p) => p.message).join('; ')))];
  });
  main.append(el('p', 'u2-gallery-status',
    'Pick a property and the operator list follows it; switching kinds swaps the value editor. ' +
    'The footer is the canonical string — paste it into any `filter()` call.'));

  section('Horizontal (filters on top of a view) — the same tree', () => {
    const fb = new FilterBuilder({schema: SCHEMA, bind: tree, orientation: 'horizontal', inline: true,
      name: 'criteriaH'});
    const width = signal(720);
    const host = divV([fb]);
    const scope = Scope.ambient;
    scope.effect(() => host.style.width = `${width.value}px`);
    return [host, divH([
      button('720px', () => width.value = 720), button('400px', () => width.value = 400),
    ]), readout('query', fb.query)];
  });
  main.append(el('p', 'u2-gallery-status',
    'Both builders are bound to one signal: edit either and the other follows. Click the ' +
    'and/or chip between two rows to flip the connector; narrow the host to see the rows wrap.'));

  section('Every editor kind', () => {
    const fb = new FilterBuilder({schema: SCHEMA, showQuery: true, name: 'kinds', value: Filters.group('or', [
      Filters.cond('name', 'like', 'pro'),
      Filters.cond('mw', 'between', [200, 500]),
      Filters.cond('status', 'in', ['Open', 'Blocked']),
      Filters.cond('active', '=', true),
      Filters.cond('created', '>', {span: '-1w'}),
      Filters.cond('name', 'is null'),
    ])});
    return [fb];
  });
  main.append(el('p', 'u2-gallery-status',
    'Contains → text with suggestions, between → two numbers, in → tags over the schema values, ' +
    'bool → true/false, datetime → calendar or a typed span (-1w, 2d, now), is empty → no editor.'));

  section('Advanced (nested, drag to reorder)', () => {
    const fb = new FilterBuilder({schema: SCHEMA, mode: 'advanced', showQuery: true, name: 'advanced',
      value: Filters.group('and', [
        Filters.cond('age', '>', 30),
        Filters.group('or', [
          Filters.cond('status', 'in', ['Open', 'Blocked']),
          Filters.group('and', [Filters.cond('sex', '=', 'F'), Filters.cond('mw', 'between', [200, 500])]),
        ]),
        Filters.cond('created', '>', {span: '-1w'}),
      ])});
    return [fb, readout('status', computed(() => JSON.stringify({mode: fb.mode.value, problems: fb.problems.value.length})))];
  });
  main.append(el('p', 'u2-gallery-status',
    'Each group has its own and/or, a `not` chip, + condition, + group and −; the left rule marks ' +
    'the depth. Drag a row by its grip onto a group header to move it inside, or between rows to ' +
    'reorder; Esc cancels. "Simple" is refused while groups exist unless they can be flattened.'));

  section('Locked template — only the value of the first row may change, nothing added', () => {
    const template = {
      root: Filters.group('and', [
        Filters.cond('status', 'in', ['Open'], {lock: 'value'}),
        Filters.cond('created', '>', {span: '-30d'}, {lock: 'all'}),
      ]),
      allowedProperties: ['status', 'created', 'name'],
      allowAdd: false,
    };
    const fb = new FilterBuilder({schema: SCHEMA, template, value: Filters.applyTemplate(template),
      showQuery: true, name: 'template'});
    return [fb, readout('diff', computed(() => JSON.stringify(Filters.diff(template.root, fb.value.value))))];
  });

  section('Locked template, advanced — a locked group keeps its shape, the free group is yours', () => {
    const template = {
      root: Filters.group('and', [
        Filters.cond('age', '>', 18),
        Filters.group('or', [
          Filters.cond('status', '=', 'Open'), Filters.cond('status', '=', 'Blocked'),
        ], {lock: 'value'}),
        Filters.group('or', [Filters.cond('name', 'like', 'pro')]),
      ]),
      allowedOperators: {age: ['>', '>=', 'between']},
    };
    const fb = new FilterBuilder({schema: SCHEMA, mode: 'advanced', template, value: Filters.applyTemplate(template),
      showQuery: true, name: 'templateAdvanced'});
    return [fb, readout('diff', computed(() => JSON.stringify(Filters.diff(template.root, fb.value.value)))),
      readout('problems', computed(() => fb.problems.value.map((p) => p.message).join('; ') || 'none'))];
  });
  main.append(el('p', 'u2-gallery-status',
    'The locked group inherits its lock to its rows (values editable, structure fixed, no grip); ' +
    'the free group can be edited, negated and dragged. The diff is what a saved template records.'));

  section('Simple-only template — allowAdvanced: false hides the mode switch', () => {
    const template = {root: Filters.group('and', [Filters.cond('active', '=', true)]), allowAdvanced: false};
    const fb = new FilterBuilder({schema: SCHEMA, template, value: Filters.applyTemplate(template), name: 'simpleOnly'});
    return [fb];
  });

  section('Empty, then add', () => {
    const fb = new FilterBuilder({label: 'Filter', schema: SCHEMA, showQuery: true, name: 'empty'});
    return [fb, readout('status', computed(() => JSON.stringify({
      query: fb.query.value, problems: fb.problems.value.length, mode: fb.mode.value})))];
  });
  main.append(el('p', 'u2-gallery-status',
    'A fresh row is a problem until it has a value (the red rule and the line under the row say which); ' +
    'the builder\'s validity is the first problem.'));

  section('Nested tree in simple mode', () => {
    const fb = new FilterBuilder({schema: SCHEMA, showQuery: true, name: 'nested', value: Filters.group('and', [
      Filters.cond('age', '>', 30),
      Filters.group('or', [Filters.cond('sex', '=', 'F'), Filters.cond('name', 'like', 'an')]),
    ])});
    return [fb, el('p', 'u2-gallery-status',
      'A sub-group shows as a read-only summary row with the hint above; Flatten is offered only ' +
      'when the connectors allow it (here they do not), Advanced always.')];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
