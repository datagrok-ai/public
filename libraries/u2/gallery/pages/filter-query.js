import {signal, computed, Scope, Control, Filters, FilterBuilder, FilterQueryInput} from '../../src/index.js';
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
injectOnce('u2-filter-query-css', '../../css/filter-query.css');

/** The same literal schema as the builder page — plain IProperty objects plus value lists. */
const SCHEMA = Filters.schema([
  {name: 'name', type: 'string', friendlyName: 'Name'},
  {name: 'age', type: 'int', friendlyName: 'Age', min: 0, max: 120},
  {name: 'mw', type: 'double', friendlyName: 'Mol weight'},
  {name: 'sex', type: 'string', friendlyName: 'Sex', choices: ['F', 'M']},
  {name: 'status', type: 'string', friendlyName: 'Status'},
  {name: 'active', type: 'bool', friendlyName: 'Active'},
  {name: 'created', type: 'datetime', friendlyName: 'Created'},
], {
  name: ['Aspirin', 'Ibuprofen', 'Paracetamol', 'Naproxen', 'Diclofenac'],
  status: ['Open', 'Blocked', 'Closed', 'Draft'],
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
  main.append(el('h1', null, 'Filter query input'));
  const intro = el('p');
  intro.innerHTML = 'One-line query box over the filter grammar with completion at the caret: ' +
    'properties, then the property\'s operators, then its values (from the schema, or literal ' +
    'hints per kind), then <code>and</code>/<code>or</code>. Enter or blur parses the text into ' +
    'the same immutable <code>FilterGroup</code> the builder edits; a tree written elsewhere ' +
    're-formats the text. Everything here runs on a literal schema — no platform.';
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
  section('Query input and builder on one tree', () => {
    const q = new FilterQueryInput({label: 'Query', schema: SCHEMA, bind: tree, name: 'query',
      placeholder: 'name like "asp" and created > -1w'});
    const fb = new FilterBuilder({label: 'Criteria', schema: SCHEMA, bind: tree, name: 'criteria'});
    return [q, readout('problems', computed(() => q.problems.value.length === 0 ? 'none' :
      q.problems.value.map((p) => p.message).join('; '))), fb];
  });
  main.append(el('p', 'u2-gallery-status',
    'Type `age > 30 and (sex = "F" or name like "an")` then Enter: the builder shows the nested ' +
    'group. Edit a row in the builder and the text follows; a text with a problem keeps the old ' +
    'tree and says what is wrong.'));

  section('Completion at every position', () => {
    const q = new FilterQueryInput({schema: SCHEMA, inline: true, name: 'complete',
      placeholder: 'start typing a property…'});
    const host = divV([q]);
    host.style.width = '480px';
    return [host, readout('query', q.query)];
  });
  main.append(el('p', 'u2-gallery-status',
    'Down arrow opens the list; Enter or a click inserts the row and moves on to the next ' +
    'expectation. Strings are quoted for you; `name =` offers the schema values, `active =` ' +
    'true/false, `created >` relative spans such as -1w and now.'));

  main.append(el('h2', null, 'Disposal'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
