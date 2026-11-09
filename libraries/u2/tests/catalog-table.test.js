/* CatalogTable — and the BasicTable it stands on: live columns, sortable headers, foldable sections —
   on the node DOM shim. Every test leaves the live-scope count where it found it. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {layout} from './layout.mjs';
import {signal, computed, Scope, AsyncValue, SortOrder} from '../src/index.js';
import {BasicTable} from '../src/components/collections/table.js';
import {CatalogTable} from '../src/components/collections/catalog-table.js';

function scoped(name, body) {
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

// ---- BasicTable ------------------------------------------------------------------------------------

const DRUGS = [
  {name: 'Aspirin', mw: 180.16, kind: 'Analgesic'},
  {name: 'Caffeine', mw: 194.19, kind: 'Stimulant'},
  {name: 'Ibuprofen', mw: 206.28, kind: 'Analgesic'},
];

const nameCol = {header: 'Name', name: 'name', render: (d) => d.name, sortable: true};
const mwCol = {header: 'MW', name: 'mw', render: (d) => String(d.mw), align: 'right', sortable: true};

const headers = (t) => [...t.root.querySelectorAll('th[data-u2="column-header"]')].map((th) => th.textContent);
const rowKeys = (t) => [...t.root.querySelectorAll('.u2-table-row')].map((r) => r.dataset.u2Item);

scoped('BasicTable: live columns build the header and the rows again; the rows keep their keys', () => {
  const columns = signal([nameCol]);
  const t = new BasicTable({columns, items: DRUGS, keyOf: (d) => d.name});
  assert.deepEqual(headers(t), ['Name']);
  assert.deepEqual([...t.root.querySelectorAll('.u2-table-row')[0].querySelectorAll('td')].map((td) => td.textContent),
    ['Aspirin']);
  columns.value = [nameCol, mwCol];
  assert.deepEqual(headers(t), ['Name', 'MW']);
  const cells = [...t.root.querySelectorAll('.u2-table-row')[0].querySelectorAll('td')];
  assert.deepEqual(cells.map((td) => td.textContent), ['Aspirin', '180.16'], 'the rows are built for the new columns');
  assert.equal(cells[1].dataset.column, 'mw');
  assert.equal(cells[1].classList.contains('u2-table-align-right'), true);
  assert.deepEqual(rowKeys(t), ['Aspirin', 'Caffeine', 'Ibuprofen']);
  t.columns = [mwCol];
  assert.deepEqual(headers(t), ['MW'], 'the setter takes a list as well');
  t.dispose();
});

scoped('BasicTable: a sortable header cycles ascending, descending and none, with the arrow and aria-sort', () => {
  const sort = signal('');
  const t = new BasicTable({columns: [nameCol, {...mwCol, sortable: false}], items: DRUGS, sort});
  const [name, mw] = t.root.querySelectorAll('th');
  assert.equal(name.classList.contains('u2-table-sortable'), true);
  assert.equal(mw.classList.contains('u2-table-sortable'), false, 'not sortable: no gesture');
  assert.equal(name.tabIndex, 0);
  fire(name, 'click');
  assert.equal(sort.value, 'name');
  assert.equal(name.getAttribute('aria-sort'), 'ascending');
  assert.equal(name.querySelectorAll('.u2-table-sort').length, 1);
  assert.equal(name.classList.contains('u2-table-sorted'), true);
  fire(name, 'keydown', {key: 'Enter'});
  assert.equal(sort.value, '!name', 'Enter sorts too');
  assert.equal(name.getAttribute('aria-sort'), 'descending');
  fire(name, 'click');
  assert.equal(sort.value, '');
  assert.equal(name.getAttribute('aria-sort'), 'none');
  assert.equal(name.querySelectorAll('.u2-table-sort').length, 0);
  fire(mw, 'click');
  assert.equal(sort.value, '', 'a header that does not sort does nothing');
  sort.value = '!name';
  assert.equal(name.getAttribute('aria-sort'), 'descending', 'the header follows the signal');
  t.dispose();
});

scoped('BasicTable: sortGesture dblclick sorts on a double click only', () => {
  const sort = signal('');
  const t = new BasicTable({columns: [nameCol], items: DRUGS, sort, sortGesture: 'dblclick'});
  const th = t.root.querySelector('th');
  fire(th, 'click');
  assert.equal(sort.value, '');
  fire(th, 'dblclick');
  assert.equal(sort.value, 'name');
  t.dispose();
});

test('SortOrder: one string for every table — parse, next, of', () => {
  assert.equal(SortOrder.next('name', 'name'), '!name');
  assert.equal(SortOrder.next('!name', 'name'), '');
  assert.equal(SortOrder.next('mw', 'name'), 'name');
  assert.equal(SortOrder.next('name,!mw', 'name'), '!name', 'the first key is the one a header sorts by');
  assert.deepEqual(SortOrder.parse('!Created,Code'), {column: 'Created', descending: true});
  assert.equal(SortOrder.parse(''), null);
  assert.equal(SortOrder.of('Code', true), '!Code');
});

scoped('BasicTable: onSort and onCollapse take what the headers ask for, in place of writing the signals', () => {
  const asked = [];
  const sort = signal('');
  const collapsed = signal([]);
  const t = new BasicTable({columns: [nameCol], items: [{section: 'Acids'}, ...DRUGS], sort, collapsed,
    onSort: (s) => asked.push(['sort', s]), onCollapse: (c) => asked.push(['collapse', c])});
  fire(t.root.querySelector('th'), 'click');
  fire(t.root.querySelector('.u2-table-section th'), 'click');
  assert.deepEqual(asked, [['sort', 'name'], ['collapse', ['Acids']]]);
  assert.equal(sort.value, '', 'the signals are the owner\'s to write');
  assert.deepEqual(collapsed.value, []);
  t.dispose();
});

scoped('BasicTable: contextActions — a row\'s actions as its context menu, by the pointer or the keyboard', () => {
  const ran = [];
  const t = new BasicTable({columns: [nameCol], items: DRUGS,
    contextActions: (item) => [{text: `Open ${item.name}`, run: () => ran.push(item.name)}]});
  document.body.append(t.root);
  const row = t.root.querySelectorAll('.u2-table-row')[1];
  fire(row.querySelector('td'), 'contextmenu', {clientX: 5, clientY: 5});
  assert.deepEqual([...document.querySelectorAll('.u2-menu-item')].map((i) => i.textContent.trim()),
    [`Open ${DRUGS[1].name}`]);
  fire(document.querySelector('.u2-menu-item'), 'click');
  assert.deepEqual(ran, [DRUGS[1].name]);
  fire(row.querySelector('td'), 'keydown', {key: 'F10', shiftKey: true});
  assert.equal(document.querySelectorAll('.u2-menu').length, 1, 'Shift+F10 opens it under the row');
  t.dispose();
});

scoped('BasicTable: a row under a folded section is revealed by its first cell\'s text too', async () => {
  const collapsed = signal(['Acids']);
  const t = new BasicTable({columns: [nameCol], items: [{section: 'Acids'}, ...DRUGS], collapsed});
  document.body.append(t.root);
  const found = await t.revealItem(DRUGS[0].name.toUpperCase());
  assert.deepEqual(collapsed.value, [], 'the section unfolds');
  assert.equal(found?.dataset.u2Item, DRUGS[0].name);
  t.dispose();
});

scoped('BasicTable: foldable sections chevron, count, aria-expanded; a fold hides their rows; reveal unfolds', async () => {
  const collapsed = signal([]);
  const t = new BasicTable({columns: [nameCol, mwCol], keyOf: (d) => d.name, collapsed, sectionCounts: true,
    sectionLabel: (title) => `${title}s`,
    items: [{section: 'Analgesic'}, DRUGS[0], DRUGS[2], {section: 'Stimulant'}, DRUGS[1]]});
  const section = t.root.querySelectorAll('.u2-table-section')[0];
  assert.equal(section.getAttribute('aria-expanded'), 'true');
  assert.equal(section.querySelector('[data-u2-part="item-label"]').textContent, 'Analgesics');
  assert.equal(section.querySelector('.u2-badge').textContent, '2');
  assert.equal(section.querySelector('.u2-table-chevron').classList.contains('fa-angle-down'), true);
  assert.equal(section.querySelector('th').colSpan, 2);
  fire(section.querySelector('th'), 'click');
  assert.deepEqual(collapsed.value, ['Analgesic']);
  assert.deepEqual(rowKeys(t), ['Caffeine'], 'the folded rows are not rendered');
  const again = t.root.querySelectorAll('.u2-table-section')[0];
  assert.equal(again, section, 'the header is patched in place, so it keeps the focus');
  assert.equal(section.getAttribute('aria-expanded'), 'false');
  assert.equal(section.querySelector('.u2-table-chevron').classList.contains('fa-angle-right'), true);
  assert.equal(section.querySelector('.u2-badge').textContent, '2', 'the count still counts the folded rows');
  fire(t.root.querySelectorAll('.u2-table-section')[1].querySelector('th'), 'keydown', {key: ' '});
  assert.deepEqual(collapsed.value, ['Analgesic', 'Stimulant'], 'Space folds the focused header');
  const row = await t.revealItem('Ibuprofen');
  assert.equal(row?.dataset.u2Item, 'Ibuprofen', 'revealing a folded row unfolds its section');
  assert.deepEqual(collapsed.value, ['Stimulant']);
  t.dispose();
});

scoped('BasicTable: pinned and capped columns, and the form row height', () => {
  const t = new BasicTable({columns: [{...nameCol, pinned: 'start', maxWidth: '200px'}, {...mwCol, pinned: 'end'}],
    items: DRUGS, rowHeight: 'control'});
  assert.equal(t.root.classList.contains('u2-table-control'), true);
  const [first, last] = t.root.querySelectorAll('.u2-table-row')[0].querySelectorAll('td');
  assert.equal(first.classList.contains('u2-table-pin-start'), true);
  assert.equal(first.classList.contains('u2-table-clip'), true);
  assert.equal(first.style.maxWidth, '200px');
  assert.equal(last.classList.contains('u2-table-pin-end'), true);
  layout(first, {scrollWidth: 300, clientWidth: 200});
  fire(first, 'pointerover');
  assert.equal(first.title, 'Aspirin', 'a cut cell shows its whole text');
  layout(first, {scrollWidth: 100});
  fire(first, 'pointerover');
  assert.equal(first.hasAttribute('title'), false);
  t.dispose();
});

// ---- CatalogTable ----------------------------------------------------------------------------------

function campaigns() {
  const c = (code, name, extra) => ({code, name, template: code.split('-')[0], status: 'In Progress',
    created: '2026/1/1', molecules: 0, modified: 0, fields: {}, ...extra});
  return [
    c('KHL-2', 'JAK2 selectivity', {created: '2026/3/12', molecules: 420, modified: 5, fields: {Chemist: 'Ann', MW: 350}}),
    c('KHL-12', 'PI3K quick wins', {created: '2026/10/19', molecules: 3, modified: 9, fields: {Chemist: 'Cy', MW: 90}}),
    c('KHL-13', 'BTK follow-up', {created: '2026/3/30', molecules: 260, modified: 1, status: 'Submitted'}),
    c('GPCR-1', 'A2A antagonists', {created: '2026/2/20', molecules: 750, modified: 7, status: ''}),
    c('FRAG-3', 'Mpro merges', {created: '2026/2/27', molecules: undefined, modified: 3, status: 'Completed',
      fields: {Chemist: 'Ann', MW: 1000}}),
  ];
}

const TEMPLATES = {KHL: 'Kinase hit-to-lead', GPCR: 'GPCR series', FRAG: 'Fragment expansion'};

function catalogColumns(extra = []) {
  return [
    {name: 'Name', value: (c) => c.name, pinned: true, link: 'Continue Campaign', maxWidth: '320px'},
    {name: 'Code', value: (c) => c.code},
    {name: 'Created', value: (c) => c.created},
    {name: 'Molecules', value: (c) => c.molecules, align: 'right'},
    {name: 'Status', value: (c) => c.status},
    ...extra,
  ];
}

function groupings() {
  return [
    {name: 'Template', key: (c) => c.template, groupLabel: (key) => TEMPLATES[key] ?? key},
    {name: 'Status', key: (c) => c.status},
    {name: 'fields.Chemist', label: 'Chemist', key: (c) => c.fields.Chemist, menuGroup: 'Campaign Fields'},
    {name: 'fields.MW', label: 'MW', key: (c) => c.fields.MW, menuGroup: 'Campaign Fields',
      compare: (a, b) => Number(a) - Number(b)},
  ];
}

function catalog(options = {}) {
  const t = new CatalogTable({name: 'campaigns', label: 'Campaigns', items: campaigns(), keyOf: (c) => c.code,
    columns: catalogColumns(), groupings: groupings(), defaultColumns: ['Code', 'Created', 'Molecules', 'Status'],
    defaultOrder: (a, b) => b.modified - a.modified, ...options});
  document.body.append(t.root);
  return t;
}

const codes = (t) => [...t.root.querySelectorAll('.u2-table-row')].map((r) => r.dataset.u2Item);
const header = (t, name) => [...t.root.querySelectorAll('th')].find((th) => th.dataset.column === name);

scoped('CatalogTable: columns recomputed as new objects of the same shape keep the header and the rows', () => {
  const items = signal(campaigns());
  const t = catalog({items, columns: computed(() => {
    items.value;
    return catalogColumns();
  })});
  const head = t.root.querySelector('thead');
  const rows = t.root.querySelectorAll('.u2-table-row');
  items.value = [...items.peek()];
  assert.ok(t.root.querySelector('thead') === head, 'the header is not built again');
  assert.ok([...t.root.querySelectorAll('.u2-table-row')].every((r, i) => r === rows[i]), 'nor the rows');
  assert.equal(t.getState().extra.hidden, '', 'the shown columns are told by name');
  t.dispose();
});
const menuItems = () => [...document.querySelectorAll('.u2-menu-item')].map((el) => ({
  label: el.querySelector('.u2-menu-label').textContent, checked: el.getAttribute('aria-checked'),
  detail: el.querySelector('.u2-menu-detail')?.textContent}));

scoped('CatalogTable: the pinned name first, the default columns, the default order, the kind and its name', () => {
  const t = catalog();
  assert.equal(t.root.dataset.u2, 'catalog-table');
  assert.equal(t.root.dataset.u2Name, 'campaigns');
  assert.equal(t.root.getAttribute('aria-label'), 'Campaigns');
  assert.deepEqual(headers(t), ['Name', 'Code', 'Created', 'Molecules', 'Status']);
  assert.deepEqual(codes(t), ['KHL-12', 'GPCR-1', 'KHL-2', 'FRAG-3', 'KHL-13'], 'last modified first');
  const name = t.root.querySelector('.u2-table-row td');
  assert.equal(name.classList.contains('u2-table-pin-start'), true);
  assert.equal(name.querySelector('a').textContent, 'PI3K quick wins');
  assert.equal(t.root.querySelectorAll('.u2-table-row')[0].classList.contains('u2-row'), true);
  assert.equal(t.root.querySelector('.u2-table').classList.contains('u2-table-control'), true, '28px rows');
  t.dispose();
});

scoped('CatalogTable: a double-click on a header sorts up, down, then back to the default order', () => {
  const t = catalog();
  const molecules = header(t, 'Molecules');
  fire(molecules, 'click');
  assert.equal(t.state.value.sort, '', 'a single click does nothing');
  fire(molecules, 'dblclick');
  assert.equal(t.state.value.sort, 'Molecules');
  assert.deepEqual(codes(t), ['KHL-12', 'KHL-13', 'KHL-2', 'GPCR-1', 'FRAG-3'], 'numbers, the empty count last');
  fire(header(t, 'Molecules'), 'dblclick');
  assert.deepEqual(codes(t), ['GPCR-1', 'KHL-2', 'KHL-13', 'KHL-12', 'FRAG-3'], 'descending, the empty count last too');
  assert.equal(header(t, 'Molecules').getAttribute('aria-sort'), 'descending');
  fire(header(t, 'Molecules'), 'dblclick');
  assert.equal(t.state.value.sort, '');
  assert.deepEqual(codes(t), ['KHL-12', 'GPCR-1', 'KHL-2', 'FRAG-3', 'KHL-13']);
  t.sortBy('Created', false);
  assert.deepEqual(codes(t), ['KHL-12', 'KHL-13', 'KHL-2', 'FRAG-3', 'GPCR-1'], 'YYYY/M/D sorts as dates');
  assert.equal(header(t, 'Created').getAttribute('aria-sort'), 'descending', 'code writes reach the header');
  t.dispose();
});

scoped('CatalogTable: codes sort as text in natural order — KHL-12 is no date', () => {
  const t = catalog();
  t.sortBy('Code');
  assert.deepEqual(codes(t), ['FRAG-3', 'GPCR-1', 'KHL-2', 'KHL-12', 'KHL-13']);
  const sorted = (values) => CatalogTable.sortItems(values.map((v) => ({v})), {name: 'x', value: (r) => r.v}, true)
    .map((r) => r.v);
  assert.deepEqual(sorted(['KHL-12', 'KHL-1', 'KHL-2']), ['KHL-1', 'KHL-2', 'KHL-12'], 'text, in natural order');
  assert.deepEqual(sorted(['2026/3/30', null, '2026-10-01', '', '2025/1/2']), ['2025/1/2', '2026/3/30', '2026-10-01',
    null, ''], 'dates in either spelling, the empty values last');
  assert.deepEqual(sorted([new Date(2026, 1, 1), undefined, new Date(2025, 1, 1)]).map((d) => d?.getFullYear()),
    [2025, 2026, undefined]);
  assert.deepEqual(sorted(['12', 3.5, '', '1e3']), [3.5, '12', '1e3', ''], 'numbers');
  t.sortBy('Status');
  assert.deepEqual(codes(t), ['FRAG-3', 'KHL-12', 'KHL-2', 'KHL-13', 'GPCR-1'],
    'text; equal values keep the default order; the empty status last');
  t.dispose();
});

scoped('CatalogTable: a sorter replaces the comparator and gets the items in the default order', () => {
  const calls = [];
  const t = catalog({sorter: (items, column, ascending) => {
    calls.push([items.map((c) => c.code).join(','), column.name, ascending]);
    return [...items].reverse();
  }});
  t.sortBy('Code', false);
  assert.deepEqual(calls, [['KHL-12,GPCR-1,KHL-2,FRAG-3,KHL-13', 'Code', false]]);
  assert.deepEqual(codes(t), ['KHL-13', 'FRAG-3', 'KHL-2', 'GPCR-1', 'KHL-12']);
  t.dispose();
});

scoped('CatalogTable: groups over every item, labelled, counted, in order of appearance, the empty one last', () => {
  const t = catalog();
  t.groupBy('Template');
  const sections = [...t.root.querySelectorAll('.u2-table-section')];
  assert.deepEqual(sections.map((s) => s.querySelector('[data-u2-part="item-label"]').textContent),
    ['Kinase hit-to-lead', 'GPCR series', 'Fragment expansion'], 'labelled by the template name, keyed by its key');
  assert.deepEqual(sections.map((s) => s.dataset.u2Item), ['section:KHL', 'section:GPCR', 'section:FRAG']);
  assert.deepEqual(sections.map((s) => s.querySelector('.u2-badge').textContent), ['3', '1', '1']);
  assert.deepEqual(sections.map((s) => s.getAttribute('aria-expanded')), ['true', 'true', 'true'], 'groups start open');
  assert.deepEqual(codes(t), ['KHL-12', 'KHL-2', 'KHL-13', 'GPCR-1', 'FRAG-3']);
  t.groupBy('Status');
  assert.deepEqual([...t.root.querySelectorAll('.u2-table-section [data-u2-part="item-label"]')].map((s) => s.textContent),
    ['In Progress', 'Completed', 'Submitted', 'No information']);
  t.groupBy('fields.MW');
  assert.deepEqual([...t.root.querySelectorAll('.u2-table-section')].map((s) => s.dataset.u2Item),
    ['section:90', 'section:350', 'section:1000', 'section:No information'], 'a grouping orders its groups its way');
  t.dispose();
});

scoped('CatalogTable: one group shows no header unless asked; with no grouping there are none', () => {
  const items = campaigns().filter((c) => c.template === 'KHL');
  const t = catalog({items});
  t.groupBy('Template');
  assert.equal(t.root.querySelectorAll('.u2-table-section').length, 0);
  assert.equal(codes(t).length, 3);
  t.dispose();
  const single = catalog({items, singleGroupHeader: true, emptyGroupLabel: 'No Information'});
  single.groupBy('Template');
  assert.equal(single.root.querySelectorAll('.u2-table-section').length, 1);
  single.groupBy('fields.Chemist');
  assert.equal(single.root.querySelectorAll('.u2-table-section')[2].textContent.startsWith('No Information'), true);
  single.dispose();
});

scoped('CatalogTable: a folded group keeps its rows out, per grouping, and says so in the state', async () => {
  const t = catalog();
  t.groupBy('Template');
  fire(t.root.querySelectorAll('.u2-table-section')[0].querySelector('th'), 'click');
  assert.deepEqual(t.state.value.collapsed, {Template: ['KHL']});
  assert.deepEqual(codes(t), ['GPCR-1', 'FRAG-3']);
  const items = t.getItems().list;
  assert.deepEqual(items.map((i) => [i.key, i.rendered, i.expanded]).slice(0, 4),
    [['section:KHL', true, false], ['KHL-12', false, undefined], ['KHL-2', false, undefined], ['KHL-13', false, undefined]]);
  assert.equal(items[1].label, 'PI3K quick wins', 'a row reads as its first column');
  t.groupBy('Status');
  assert.equal(codes(t).length, 5, 'another grouping has its own folds');
  t.groupBy('Template');
  assert.deepEqual(codes(t), ['GPCR-1', 'FRAG-3'], 'and the first one comes back with its own');
  const row = await t.revealItem('BTK follow-up');
  assert.equal(row?.dataset.u2Item, 'KHL-13', 'a row in a folded group is revealed by its label, unfolding the group');
  assert.deepEqual(t.state.value.collapsed, {Template: []});
  assert.equal(await t.revealItem('Nothing like it'), null);
  t.dispose();
});

scoped('CatalogTable: the state block reads sort, grouping, columns, folds, groups and counts', () => {
  const t = catalog({columns: catalogColumns([{name: 'Chemist', value: (c) => c.fields.Chemist, menuGroup: 'Fields'}])});
  t.sortBy('Molecules', false);
  t.groupBy('Template');
  t.foldGroup('GPCR');
  const state = t.getState();
  assert.equal(state.kind, 'catalog-table');
  assert.equal(state.name, 'campaigns');
  assert.equal(state.label, 'Campaigns');
  assert.deepEqual(state.extra, {sort: 'Molecules desc', groupBy: 'Template', columns: 'Name,Code,Created,Molecules,Status',
    hidden: 'Chemist', collapsed: 'GPCR', groups: 'GPCR series (1), Kinase hit-to-lead (3), Fragment expansion (1)',
    rows: 4, total: 5}, 'the groups come in the order their first row has in the sort');
  assert.equal(state.items.count, 8, 'three group headers and five rows');
  t.dispose();
});

scoped('CatalogTable: the columns menu toggles in toggle order, a menuGroup is a submenu with details', () => {
  const t = catalog({columns: catalogColumns([
    {name: 'fields.Chemist', header: 'Chemist', value: (c) => c.fields.Chemist, menuGroup: 'Campaign Fields',
      menuDetail: 'Kinase hit-to-lead, Fragment expansion'}])});
  const menu = t.build(() => t.columnsMenu());
  menu.show({x: 0, y: 0});
  assert.deepEqual(menuItems().slice(0, 5), [
    {label: 'Code', checked: 'true', detail: undefined}, {label: 'Created', checked: 'true', detail: undefined},
    {label: 'Molecules', checked: 'true', detail: undefined}, {label: 'Status', checked: 'true', detail: undefined},
    {label: 'Campaign Fields', checked: null, detail: undefined}]);
  fire(document.querySelectorAll('.u2-menu-item')[1], 'click');
  assert.deepEqual(t.state.value.columns, ['Code', 'Molecules', 'Status']);
  assert.equal(menuItems()[1].checked, 'false', 'the menu stays open, Created unchecked in it');
  fire(document.querySelectorAll('.u2-menu-item')[1], 'click');
  assert.deepEqual(t.state.value.columns, ['Code', 'Molecules', 'Status', 'Created'], 'picked again: shown, at the end');
  assert.equal(menuItems()[1].checked, 'true');
  menu.close();
  t.showColumn('Created', false);
  t.showColumn('fields.Chemist');
  t.showColumn('Created');
  assert.deepEqual(t.state.value.columns, ['Code', 'Molecules', 'Status', 'fields.Chemist', 'Created'],
    'shown columns come in the order they were turned on');
  assert.deepEqual([...t.root.querySelectorAll('thead th')].map((th) => th.textContent),
    ['Name', 'Code', 'Molecules', 'Status', 'Chemist', 'Created']);
  t.dispose();
});

scoped('CatalogTable: the group menu checks the current grouping and picks another', () => {
  const t = catalog();
  t.build(() => t.groupMenu()).show({x: 0, y: 0});
  assert.deepEqual(menuItems().map((i) => [i.label, i.checked]),
    [['None', 'true'], ['Template', 'false'], ['Status', 'false'], ['Campaign Fields', null]]);
  fire(document.querySelectorAll('.u2-menu-item')[2], 'click');
  assert.equal(t.state.value.groupBy, 'Status');
  t.dispose();
});

scoped('CatalogTable: choices are loaded once, saved on each change, and junk from storage is ignored', () => {
  const saved = [];
  const persist = {
    load: () => ({columns: ['Status', 'fields.Gone', 7], sort: '!Code', groupBy: 'Status',
      collapsed: {Status: ['Completed'], Template: 'oops'}}),
    save: (s) => saved.push(JSON.parse(JSON.stringify(s))),
  };
  const t = catalog({persist});
  assert.deepEqual(t.state.value, {columns: ['Code', 'Created', 'Molecules', 'Status'],
    sort: '!Code', groupBy: 'Status', collapsed: {Status: ['Completed']}},
  'a column list with a non-string falls back to the defaults');
  assert.deepEqual(saved, [], 'loading saves nothing');
  t.groupBy(null);
  assert.equal(saved.length, 1);
  assert.equal(saved[0].groupBy, null);
  t.dispose();

  const kept = catalog({persist: {load: () => ({columns: ['Status', 'fields.Later']}), save: () => {}}});
  assert.deepEqual([...kept.root.querySelectorAll('thead th')].map((th) => th.textContent), ['Name', 'Status'],
    'a column that does not exist yet is kept, not shown');
  kept.dispose();

  for (const [old, now] of [[{column: 'Code', ascending: false}, '!Code'], [{column: 'Code', ascending: true}, 'Code'],
    [null, '']]) {
    const legacy = catalog({persist: {load: () => ({sort: old}), save: () => {}}});
    assert.equal(legacy.state.value.sort, now, `a sort kept as ${JSON.stringify(old)} still reads`);
    legacy.dispose();
  }
});

scoped('CatalogTable: a kept column shows once live columns bring it', () => {
  const columns = signal(catalogColumns());
  const t = catalog({columns, state: {columns: ['Status', 'fields.Chemist']}});
  assert.deepEqual([...t.root.querySelectorAll('thead th')].map((th) => th.textContent), ['Name', 'Status']);
  columns.value = catalogColumns([{name: 'fields.Chemist', header: 'Chemist', value: (c) => c.fields.Chemist}]);
  assert.deepEqual([...t.root.querySelectorAll('thead th')].map((th) => th.textContent), ['Name', 'Status', 'Chemist']);
  assert.equal(t.root.querySelectorAll('.u2-table-row')[0].querySelectorAll('td')[2].textContent, 'Cy');
  t.dispose();
});

scoped('CatalogTable: storage(key) keeps the whole state as JSON in localStorage', () => {
  const store = new Map();
  const previous = Object.getOwnPropertyDescriptor(globalThis, 'localStorage');
  Object.defineProperty(globalThis, 'localStorage', {configurable: true, value: {
    getItem: (k) => store.get(k) ?? null, setItem: (k, v) => store.set(k, v)}});
  try {
    const t = catalog({persist: 'demo'});
    t.sortBy('Code');
    assert.equal(JSON.parse(store.get('u2.catalog.demo')).sort, 'Code');
    t.dispose();
    const back = catalog({persist: 'demo'});
    assert.equal(back.state.value.sort, 'Code');
    back.dispose();
  }
  finally {
    if (previous)
      Object.defineProperty(globalThis, 'localStorage', previous);
    else
      delete globalThis.localStorage;
  }
});

function actions(log) {
  return (c) => [
    {text: 'Manage campaign permissions', icon: 'users', requires: 'edit', run: () => log.push(`share ${c.code}`)},
    {text: 'Delete campaign', icon: 'trash-alt', requires: 'edit', destructive: true, run: () => log.push(`delete ${c.code}`)},
  ];
}

scoped('CatalogTable: row actions follow each row\'s access: edit, view only, still asked', () => {
  const log = [];
  const access = signal('pending');
  const t = catalog({actions: actions(log), actionsVisible: 'always',
    viewOnly: {icon: 'eye', tooltip: 'View only'},
    access: (c) => c.code === 'GPCR-1' ? 'view' : access});
  assert.equal(t.root.classList.contains('u2-catalog-actions-always'), true);
  const cellOf = (code) => [...t.root.querySelectorAll('.u2-table-row')].find((r) => r.dataset.u2Item === code)
    .querySelector('.u2-catalog-actions');
  const pending = cellOf('KHL-12').querySelector('.u2-row-actions');
  assert.equal(pending.classList.contains('u2-catalog-actions-pending'), true, 'the room is kept');
  assert.equal(pending.querySelectorAll('[data-u2="icon-button"]').length, 2);
  const view = cellOf('GPCR-1');
  assert.equal(view.querySelectorAll('[data-u2="icon-button"]').length, 0, 'view only: no edit actions');
  assert.equal(view.querySelector('.u2-catalog-view-only').getAttribute('aria-label'), 'View only');
  access.value = 'edit';
  const buttons = [...cellOf('KHL-12').querySelectorAll('[data-u2="icon-button"]')];
  assert.deepEqual(buttons.map((b) => b.getAttribute('aria-label')), ['Manage campaign permissions', 'Delete campaign']);
  assert.equal(cellOf('KHL-12').querySelector('.u2-catalog-actions-pending'), null);
  fire(buttons[1], 'click');
  assert.deepEqual(log, ['delete KHL-12']);
  const cell = cellOf('KHL-12').parentElement;
  assert.equal(cell.classList.contains('u2-table-pin-end'), true, 'the actions column stays at the end');
  t.dispose();
});

scoped('CatalogTable: a row action runs as the table\'s: what it opens outlives the row\'s icons', () => {
  const access = signal('edit');
  let owner;
  const t = catalog({access: () => access,
    actions: () => [{text: 'Delete campaign', icon: 'trash-alt', run: () => owner = Scope.ambient}]});
  fire(t.root.querySelector('.u2-table-row [data-u2="icon-button"]'), 'click');
  assert.ok(owner !== undefined, 'the action ran with an owner');
  let closed = false;
  owner.own(() => closed = true);
  access.value = 'view';
  assert.equal(closed, false, 'the icons were rebuilt: what the action opened stays open');
  t.dispose();
  assert.equal(closed, true, 'it goes with the table');
});

scoped('CatalogTable: a right click opens the row\'s actions as a menu; a view-only row has none', () => {
  const log = [];
  const t = catalog({actions: actions(log), access: (c) => c.code === 'GPCR-1' ? 'view' : 'edit'});
  const rows = t.root.querySelectorAll('.u2-table-row');
  const event = fire(rows[1].querySelector('td'), 'contextmenu');
  assert.equal(event, true, 'nothing to offer: the browser keeps its menu');
  assert.equal(document.querySelectorAll('.u2-menu-item').length, 0);
  fire(rows[0].querySelectorAll('td')[1], 'contextmenu');
  assert.deepEqual(menuItems().map((i) => i.label), ['Manage campaign permissions', 'Delete campaign']);
  fire(document.querySelectorAll('.u2-menu-item')[0], 'click');
  assert.deepEqual(log, ['share KHL-12']);
  t.dispose();
});

scoped('CatalogTable: a link column activates its row', () => {
  const opened = [];
  const t = catalog({onActivate: (c) => opened.push(c.code)});
  fire(t.root.querySelector('.u2-table-row a'), 'click');
  fire(t.root.querySelectorAll('.u2-table-row a')[1], 'keydown', {key: 'Enter'});
  assert.deepEqual(opened, ['KHL-12', 'GPCR-1']);
  t.dispose();
});

scoped('CatalogTable: loading shows a skeleton, an error says why with Retry, no items the empty state', async () => {
  const items = signal([]);
  const loading = signal(true);
  const error = signal(null);
  let retried = 0;
  const t = catalog({items, loading, error, refresh: async () => retried++, empty: 'No campaigns yet.'});
  assert.notEqual(t.root.querySelector('.u2-skeleton'), null);
  assert.equal(t.root.getAttribute('aria-busy'), 'true');
  error.value = 'Could not load campaigns: no access';
  loading.value = false;
  const failure = t.root.querySelector('.u2-empty-state-error');
  assert.equal(failure.querySelector('.u2-empty-state-text').textContent, 'Could not load campaigns: no access');
  fire(failure.querySelector('button'), 'click');
  await flush();
  assert.equal(retried, 1);
  error.value = null;
  assert.equal(t.root.querySelector('.u2-empty-state-text').textContent, 'No campaigns yet.');
  assert.equal(t.root.querySelector('table'), null);
  items.value = campaigns();
  assert.equal(codes(t).length, 5);
  assert.equal(t.root.querySelector('.u2-empty-state'), null);
  t.dispose();
});

scoped('CatalogTable: from an AsyncValue — skeleton first, then rows that stay through a refresh', async () => {
  let answer = campaigns().slice(0, 2);
  let release;
  const owner = new Scope(undefined, {label: 'test'});
  const source = owner.build(() => new AsyncValue(() => new Promise((resolve) => release = () => resolve(answer)),
    {debounceMs: 0}));
  const t = catalog({items: source, footer: (shown) => `${shown.length} campaigns`});
  const [, , sync] = t.chrome();
  assert.notEqual(t.root.querySelector('.u2-skeleton'), null);
  assert.equal(sync.getAttribute('aria-disabled'), 'true', 'the icons wait for the first load');
  await flush();
  release();
  await flush();
  assert.equal(codes(t).length, 2);
  assert.equal(t.root.querySelector('[data-u2-part="footer"]').textContent, '2 campaigns');
  answer = campaigns();
  const done = t.refresh();
  assert.equal(codes(t).length, 2, 'the rows stay while the newer answer comes');
  assert.equal(t.root.getAttribute('aria-busy'), 'true');
  assert.equal(sync.classList.contains('u2-catalog-refreshing'), true);
  await flush();
  release();
  await done;
  assert.equal(codes(t).length, 5);
  assert.equal(t.root.hasAttribute('aria-busy'), false);
  t.dispose();
  owner.dispose();
});

scoped('CatalogTable: chrome — the icons, their names and tooltips, the grouping tooltip following the state', () => {
  const t = catalog({refresh: async () => {}});
  const [columns, group, sync] = t.chrome({columns: 'Edit visible columns in campaigns table',
    group: (g) => `Group Campaigns. Current: ${g}`, refresh: 'Refresh campaigns table'});
  assert.deepEqual([columns, group, sync].map((b) => b.dataset.u2Name), ['campaigns-columns', 'campaigns-group', 'campaigns-refresh']);
  assert.equal(columns.getAttribute('aria-label'), 'Edit visible columns in campaigns table');
  assert.equal(group.getAttribute('aria-label'), 'Group Campaigns. Current: None');
  t.groupBy('fields.Chemist');
  assert.equal(group.getAttribute('aria-label'), 'Group Campaigns. Current: Chemist');
  document.body.append(columns, group, sync);
  fire(group, 'click');
  assert.deepEqual(menuItems().map((i) => i.label), ['None', 'Template', 'Status', 'Campaign Fields']);
  assert.equal(menuItems().some((i) => i.checked === 'true'), false, 'the current grouping sits in its submenu');
  t.dispose();
});

scoped('CatalogTable: a sort moves the rows it keeps, a grouping change builds none of them again', () => {
  let renders = 0;
  const many = Array.from({length: 300}, (_, i) => ({code: `C-${i}`, name: `Campaign ${i}`, template: ['A', 'B', 'C'][i % 3],
    status: 'In Progress', created: '2026/1/1', molecules: (i * 37) % 101, modified: i, fields: {}}));
  const columns = catalogColumns();
  columns[1] = {...columns[1], render: (c) => {
    renders++;
    return c.code;
  }};
  const t = catalog({items: many, columns});
  assert.equal(codes(t).length, 300);
  assert.equal(renders, 300);
  const first = t.root.querySelector('.u2-table-row');
  t.sortBy('Molecules');
  t.groupBy('Template');
  t.sortBy('Molecules', false);
  assert.equal(renders, 300, 'no row was built again');
  assert.equal(codes(t).length, 300, 'the full list, no cap');
  assert.equal([...t.root.querySelectorAll('.u2-table-row')].includes(first), true);
  t.dispose();
});
