/* FunctionsBrowser: the smart-search grammar (plain terms AND-ed, #tag, @role, internal-tagged
   items always out), tag counts over the term-filtered set, the check↔query round-trip in both
   directions, the sole-checked-tag toggle, the live show* signals, and setItems over a signal.
   Mock items only — the platform feed is the dg factory's, tested through e2e. Same contract as
   the other UI suites: every test leaves the live-scope count where it found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {signal, Scope} from '../src/index.js';
import {FunctionsBrowser, filterFuncItems, tagCounts}
  from '../src/components/collections/functions-browser.js';

function ui(name, body) {
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

function item(name, tags = [], roles = [], extra = {}) {
  const label = extra.label ?? name;
  return {name, label, tags, roles, signature: '(x) : int',
    search: `${name} ${label} ${extra.description ?? ''}`.toLowerCase(), ...extra};
}

const ITEMS = [
  item('Chem:Sketch', ['chem'], [], {label: 'Sketch molecule'}),
  item('Chem:FindSimilar', ['chem', 'search'], []),
  item('Sin', ['math'], []),
  item('Cos', ['math'], [], {description: 'The adjacent over the hypotenuse'}),
  item('Bio:Align', ['bio'], ['app']),
  item('SecretFn', ['internal'], []),
];

const NONE = new Set();
const names = (items) => items.map((x) => x.name);

function mount(options) {
  const fb = new FunctionsBrowser(options);
  document.body.append(fb.root);
  const list = fb.root.querySelector('.u2-list');
  list.clientHeight = 220;
  fire(list, 'scroll');
  return fb;
}

const tagRow = (fb, tag) =>
  fb.root.querySelector(`[data-u2="fb-tags"] .u2-fb-pane-row[data-value="${tag}"]`);

test('filterFuncItems: plain terms AND over the search string, case-insensitive', () => {
  assert.deepEqual(names(filterFuncItems(ITEMS, 'chem find', NONE, NONE)), ['Chem:FindSimilar']);
  assert.deepEqual(names(filterFuncItems(ITEMS, 'SIN', NONE, NONE)), ['Sin']);
  assert.equal(filterFuncItems(ITEMS, '', NONE, NONE).length, 5, 'an empty query filters nothing');
  assert.deepEqual(filterFuncItems(ITEMS, 'sketch cos', NONE, NONE), [], 'terms are AND-ed');
  assert.deepEqual(names(filterFuncItems(ITEMS, 'hypotenuse', NONE, NONE)), ['Cos'],
    'a term matches the description too');
});

test('filterFuncItems: #tag and @role terms filter, checked sets do the same', () => {
  assert.deepEqual(names(filterFuncItems(ITEMS, '#math', NONE, NONE)), ['Sin', 'Cos']);
  assert.deepEqual(names(filterFuncItems(ITEMS, '', new Set(['math']), NONE)), ['Sin', 'Cos']);
  assert.deepEqual(names(filterFuncItems(ITEMS, '@app', NONE, NONE)), ['Bio:Align']);
  assert.deepEqual(names(filterFuncItems(ITEMS, '', NONE, new Set(['app']))), ['Bio:Align']);
  assert.deepEqual(names(filterFuncItems(ITEMS, '#chem find', NONE, NONE)), ['Chem:FindSimilar'],
    'terms and tags combine');
  assert.deepEqual(names(filterFuncItems(ITEMS, '#chem #math', NONE, NONE)),
    ['Chem:Sketch', 'Chem:FindSimilar', 'Sin', 'Cos'], 'tags are any-of');
});

test('filterFuncItems: internal-tagged items never show, whatever the query', () => {
  assert.equal(names(filterFuncItems(ITEMS, '', NONE, NONE)).includes('SecretFn'), false);
  assert.deepEqual(filterFuncItems(ITEMS, 'secretfn', NONE, NONE), []);
  assert.deepEqual(filterFuncItems(ITEMS, '#internal', NONE, NONE), []);
});

test('tagCounts: every tag is a key, counted against the term-filtered set', () => {
  const all = tagCounts(ITEMS, '');
  assert.equal(all.get('chem'), 2);
  assert.equal(all.get('math'), 2);
  assert.equal(all.get('internal'), 0, 'internal items never count, but the tag stays a key');

  const filtered = tagCounts(ITEMS, 'sin');
  assert.equal(filtered.get('math'), 1);
  assert.equal(filtered.get('chem'), 0, 'a zero-count tag keeps its row');
  assert.equal(filtered.get('search'), 0);
});

ui('typing #tag / @role checks the panes; the typed query is never rewritten under the caret', () => {
  const fb = mount({items: ITEMS, roles: [{role: 'app', description: 'An app.'}]});
  fb.query.value = '#chem foo';
  assert.deepEqual([...fb.checkedTags.value], ['chem']);
  assert.equal(fb.query.value, '#chem foo', 'the query keeps the user\'s own order');
  assert.equal(tagRow(fb, 'chem').querySelector('input').checked, true);

  fb.query.value = '@app';
  assert.deepEqual([...fb.checkedRoles.value], ['app']);
  assert.deepEqual([...fb.checkedTags.value], [], 'a dropped #tag unchecks its box');
  fb.dispose();
});

ui('writing the checked sets rewrites the query, terms preserved', () => {
  const fb = mount({items: ITEMS});
  fb.query.value = 'foo #chem';
  fb.checkedTags.value = new Set(['math']);
  assert.equal(fb.query.value, 'foo #math');
  fb.checkedRoles.value = new Set(['app']);
  assert.equal(fb.query.value, 'foo #math @app');
  fb.checkedTags.value = new Set();
  assert.equal(fb.query.value, 'foo @app');
  fb.dispose();
});

ui('pane clicks: the checkbox toggles membership, the label checks exclusively, and a click on ' +
  'the sole checked tag unchecks it', () => {
  const fb = mount({items: ITEMS});
  fire(tagRow(fb, 'chem').querySelector('.u2-fb-pane-label'), 'click');
  assert.deepEqual([...fb.checkedTags.value], ['chem']);
  assert.equal(fb.query.value, '#chem');

  fire(tagRow(fb, 'math').querySelector('.u2-fb-pane-label'), 'click');
  assert.deepEqual([...fb.checkedTags.value], ['math'], 'a label click checks exclusively');

  const box = tagRow(fb, 'chem').querySelector('input');
  box.checked = true;
  fire(box, 'change');
  assert.deepEqual([...fb.checkedTags.value].sort(), ['chem', 'math'],
    'the checkbox itself adds to the set');
  assert.equal(tagRow(fb, 'chem').querySelector('input'), box,
    'the activated checkbox survives the update — no rebuild, no focus loss');

  fb.checkedTags.value = new Set(['math']);
  fire(tagRow(fb, 'math').querySelector('.u2-fb-pane-label'), 'click');
  assert.deepEqual([...fb.checkedTags.value], [], 'the sole checked tag toggles off');
  assert.equal(fb.query.value, '');
  fb.dispose();
});

ui('panes: counts follow the terms, ignoreTags and visibleTags shape the list, roles only when given',
  () => {
    const fb = mount({items: ITEMS});
    assert.equal(fb.root.querySelector('[data-u2="fb-roles"]'), null, 'no roles option, no pane');
    assert.equal(tagRow(fb, 'internal'), null, 'ignored tags stay out of the pane');
    assert.equal(tagRow(fb, 'chem').querySelector('.u2-fb-count').textContent, '2');
    fb.query.value = 'sin';
    assert.equal(tagRow(fb, 'chem').querySelector('.u2-fb-count').textContent, '0');
    assert.equal(tagRow(fb, 'math').querySelector('.u2-fb-count').textContent, '1');
    fb.dispose();

    const shaped = mount({items: ITEMS, visibleTags: ['math'],
      roles: [{role: 'app', description: 'An app.'}, {role: 'Transform', description: ''}]});
    assert.notEqual(shaped.root.querySelector('[data-u2="fb-roles"]'), null);
    assert.equal(shaped.root
      .querySelector('[data-u2="fb-roles"] [data-value="app"] .u2-fb-pane-label').title, 'An app.');
    assert.equal(shaped.root
      .querySelector('[data-u2="fb-roles"] [data-value="transform"] .u2-fb-pane-label').title,
    'Transform', 'a role with no description answers its name, never an empty title');
    assert.equal(tagRow(shaped, 'math').querySelector('.u2-fb-pane-label').title, 'math',
      'a truncated tag label answers its full name on hover');
    assert.notEqual(tagRow(shaped, 'math'), null);
    assert.equal(tagRow(shaped, 'chem'), null, 'visibleTags restricts the pane');
    shaped.dispose();
  });

ui('live show* toggles hide the search row, the panes, the signature and the play icon', () => {
  const sig = signal(true);
  const fb = mount({items: ITEMS, runAction: () => {}, showSignature: sig});
  const searchRow = fb.root.querySelector('[data-u2="fb-search"]');
  const panes = fb.root.querySelector('[data-u2="fb-panes"]');

  fb.showSearch.value = false;
  assert.equal(searchRow.hidden, true);
  fb.showSearch.value = true;
  assert.equal(searchRow.hidden, false);

  fb.showTags.value = false;
  assert.equal(panes.hidden, true);
  fire(fb.root.querySelector('[data-u2="fb-search"] [data-u2="icon-button"]'), 'click');
  assert.equal(fb.showTags.value, true, 'the filter icon toggles the panes back');
  assert.equal(panes.hidden, false);

  sig.value = false;
  assert.equal(fb.root.classList.contains('u2-fb-no-signature'), true,
    'a followed option signal drives the toggle');
  sig.value = true;
  assert.equal(fb.root.classList.contains('u2-fb-no-signature'), false);
  fb.showSignature.value = false;
  assert.equal(sig.value, true, 'a local toggle never writes back into the option signal');

  fb.showRunButton.value = false;
  assert.equal(fb.root.classList.contains('u2-fb-no-run'), true);
  fb.dispose();
});

ui('rows: play icon and signature; selection feeds `selected` and onChanged; dblclick/Enter activate',
  () => {
    const runs = [];
    const changes = [];
    const activated = [];
    const fb = mount({items: ITEMS, runAction: (x) => runs.push(x.name),
      onChanged: (x) => changes.push(x?.name ?? null), onActivate: (x) => activated.push(x.name)});
    const events = [];
    fb.onEvent('activate').subscribe((x) => events.push(x.name));

    const row = fb.root.querySelector('.u2-list-row[data-index="0"]');
    assert.equal(row.getAttribute('data-u2-func'), 'Chem:Sketch');
    assert.equal(row.querySelector('.u2-fb-label').textContent, 'Sketch molecule');
    assert.equal(row.querySelector('.u2-fb-sig').textContent, '(x) : int');

    fire(row, 'click');
    assert.equal(fb.selected.value.name, 'Chem:Sketch');
    assert.deepEqual(changes, ['Chem:Sketch']);

    fire(row.querySelector('.u2-fb-run'), 'click');
    assert.deepEqual(runs, ['Chem:Sketch'], 'the play icon runs the item');

    fire(row, 'dblclick');
    fire(fb.root.querySelector('.u2-list'), 'keydown', {key: 'Enter'});
    assert.deepEqual(activated, ['Chem:Sketch', 'Chem:Sketch']);
    assert.deepEqual(events, ['Chem:Sketch', 'Chem:Sketch'], 'the component event fires too');
    fb.dispose();
  });

ui('setItems over a signal: the list, the counts and the status line follow it', () => {
  const items = signal([item('Sin', ['math'])]);
  const fb = mount({});
  fb.setItems(items);
  assert.equal(fb.root.querySelector('[data-u2="fb-status"]').textContent, '1 of 1');
  assert.equal(fb.root.querySelectorAll('.u2-list-row').length, 1);

  items.value = [item('Sin', ['math']), item('Cos', ['math']), item('Chem:Sketch', ['chem'])];
  assert.equal(fb.root.querySelector('[data-u2="fb-status"]').textContent, '3 of 3');
  assert.equal(tagRow(fb, 'math').querySelector('.u2-fb-count').textContent, '2');

  fb.query.value = 'nothing-matches';
  assert.equal(fb.root.querySelector('.u2-list').hidden, true);
  const empty = fb.root.querySelector('[data-u2="fb-empty"]');
  assert.equal(empty.hidden, false);
  assert.equal(empty.querySelector('.u2-fb-empty-message').textContent,
    'No functions match "nothing-matches".');
  const clear = empty.querySelector('[data-u2="fb-clear"]');
  assert.equal(clear.hidden, false);
  fire(clear, 'click');
  assert.equal(fb.query.value, '', 'Clear search empties the query');
  assert.equal(fb.root.querySelector('.u2-list').hidden, false);
  fb.dispose();
});

ui('options: initial search + pre-checked tags merge, ignoreTags overrides the default pane set',
  () => {
    const fb = mount({items: ITEMS, search: 'foo', tags: ['chem'], ignoreTags: ['math']});
    assert.equal(fb.query.value, 'foo #chem');
    assert.deepEqual([...fb.checkedTags.value], ['chem']);
    assert.equal(tagRow(fb, 'math'), null, 'an explicit ignoreTags replaces the default');
    assert.notEqual(tagRow(fb, 'internal'), null, '…entirely — internal only hides by default');
    fb.dispose();
  });

ui('selection: a filter that drops the selected key clears it — never a phantom emission', () => {
  const changes = [];
  const fb = mount({items: ITEMS, onChanged: (x) => changes.push(x?.name ?? null)});
  fire(fb.root.querySelector('.u2-list-row[data-index="2"]'), 'click');
  assert.equal(fb.selected.value.name, 'Sin');

  fb.query.value = 'chem';
  assert.equal(fb.selected.value, null, 'the key left the filtered set');
  assert.deepEqual(changes, ['Sin', null],
    'no emission for the item that took over the stale index');
  fb.dispose();
});

ui('selection: a filter the selected key survives shifts the index silently — no re-emission', () => {
  const changes = [];
  const fb = mount({items: ITEMS, onChanged: (x) => changes.push(x?.name ?? null)});
  fire(fb.root.querySelector('.u2-list-row[data-index="3"]'), 'click');
  assert.equal(fb.selected.value.name, 'Cos');

  fb.query.value = 'cos';
  assert.equal(fb.selected.value.name, 'Cos', 'the selection survived the filter');
  assert.deepEqual(changes, ['Cos'], 'the surviving selection is not re-emitted');
  const row = fb.root.querySelector('.u2-list-row[data-index="0"]');
  assert.equal(row.getAttribute('data-u2-func'), 'Cos');
  assert.equal(row.getAttribute('aria-selected'), 'true', 'the list re-selected it at its new index');

  fb.query.value = '';
  assert.equal(fb.selected.value.name, 'Cos');
  assert.deepEqual(changes, ['Cos'], 'clearing the filter re-emits nothing either');
  fb.dispose();
});

ui('bind surface: $.fb.search aliases the query, $.fb.selected answers the qualified name', () => {
  const fb = mount({items: ITEMS});
  assert.equal(fb.bindStep('search'), fb.query, 'search is the query signal itself');
  const step = fb.bindStep('selected');
  assert.equal(step.value, null);
  fire(fb.root.querySelector('.u2-list-row[data-index="1"]'), 'click');
  assert.equal(step.value, 'Chem:FindSimilar');
  const props = fb.bindProps();
  assert.deepEqual(props.find((p) => p.name === 'search'),
    {name: 'search', type: 'string', writable: true});
  assert.deepEqual(props.find((p) => p.name === 'selected'),
    {name: 'selected', type: 'string', writable: false});
  fb.dispose();
});

ui('smoke: construct, interact, dispose — no leaked scopes', () => {
  const fb = mount({items: ITEMS, roles: [{role: 'app', description: 'App'}], runAction: () => {}});
  fb.query.value = '#chem';
  fire(fb.root.querySelector('.u2-list-row[data-index="0"]'), 'click');
  fb.dispose();
});
