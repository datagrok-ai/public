import {signal, Scope} from '../../src/index.js';
import {div, divH, span, button} from '../../src/core/elements.js';
import {DataTable} from '../../src/components/collections/data-table.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

for (const name of ['elements', 'buttons', 'icons', 'menu', 'data-table'])
  injectOnce(`u2-${name}-css`, `../../css/${name}.css`);

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

/** A line that follows a signal, owned by the control that survives the page. */
function readout(owner, label, read) {
  const line = el('p', 'u2-gallery-status');
  owner.effect(() => line.textContent = `${label} = ${read()}`);
  return line;
}

const ROW_COUNT = 200000;
const STATUS = ['Approved', 'Phase II', 'Retired'];

const rowAt = (i) => ({
  id: `c${i}`,
  name: `Compound ${i.toLocaleString()}`,
  formula: `C${9 + (i % 14)}H${8 + (i % 23)}O${1 + (i % 5)}`,
  mw: 120 + ((i * 37) % 4200) / 10,
  status: STATUS[i % STATUS.length],
});

const COLUMNS = [
  {name: 'name', header: 'Name'},
  {name: 'formula', header: 'Formula', width: '140px'},
  {name: 'mw', header: 'MW', align: 'right', width: '90px',
    render: (item) => item.mw.toFixed(1)},
  {name: 'status', header: 'Status', width: '110px'},
];

export async function render(main) {
  main.append(el('h1', null, 'Data table'));
  const intro = el('p');
  intro.innerHTML = 'Columns over the virtual list\'s render-range diff: pooled ' +
    '<code>role="row"</code> rows under a sticky header, cells pooled with their row, and only ' +
    'the visible window in the DOM. The small-data control with real <code>&lt;table&gt;</code> ' +
    'semantics is <code>BasicTable</code>; this one is for collections that do not fit.';
  main.append(intro);

  const items = Array.from({length: ROW_COUNT}, (_, i) => rowAt(i));
  const log = signal('nothing activated yet');
  const table = new DataTable({
    columns: COLUMNS,
    rowHeight: 24,
    keyOf: (item) => item.id,
    onActivate: (item, index) => log.value = `${item.name} (row ${index})`,
    contextActions: (item) => [
      {name: 'Copy name', icon: 'copy', run: () => log.value = `copied “${item.name}”`},
      {name: 'Open', run: () => log.value = `opened ${item.name}`},
    ],
  });
  table.root.style.height = '360px';
  table.setItems(items);
  main.append(table.root);

  const status = el('p', 'u2-gallery-status');
  const update = () => status.textContent =
    `${ROW_COUNT.toLocaleString()} rows · selectedIndex = ${table.selectedIndex.value}` +
    ` · selected = ${table.selectedIndices.value.size} · rows in the DOM = ${table.renderedCount}` +
    ` · live scopes = ${Scope.liveCount}`;
  table.effect(update);
  table.root.addEventListener('scroll', update);
  table.own(() => table.root.removeEventListener('scroll', update));
  main.append(status);
  main.append(readout(table, 'onActivate', () => log.value));
  main.append(div([span('Click a row (Ctrl to add, Shift for a range), right-click for its ' +
    'actions, Enter or double-click to activate, ↑ ↓ Home End PageUp PageDown to move.')],
  'u2-gallery-status'));

  const actions = divH([
    button('Scroll to 100,000', () => {
      table.scrollToIndex(100000);
      update();
    }),
    button('Replace items', () => {
      table.setItems(items.slice(1));
      update();
    }),
  ]);
  actions.style.gap = 'var(--dg-space-m)';
  main.append(actions);

  main.append(el('h2', null, 'Editor colours'));
  const editIntro = el('p');
  editIntro.innerHTML = 'A <code>cellState</code> answers per ROW KEY and column — never per frame ' +
    'index, since a row view is a keyed proxy. An unsaved edit reads as the platform grid\'s amber ' +
    '<code>DIRTY_CELL_COLOR</code>, a refusal as its red <code>INVALID_CELL_COLOR</code> with the ' +
    'message as the cell\'s tooltip, and a cell the caller may not write is muted.';
  main.append(editIntro);

  const changed = new Set(['c1.name', 'c3.mw']);
  const errors = new Map([['c2.formula', {message: 'Formula does not parse', kind: 'validation'}]]);
  const edited = new DataTable({
    columns: COLUMNS,
    rowHeight: 24,
    keyOf: (item) => item.id,
    cellState: {
      isChanged: (key, column) => changed.has(`${key}.${column}`),
      canEdit: (key, column) => column !== 'status',
      errorOf: (key, column) => errors.get(`${key}.${column}`) ?? null,
    },
  });
  edited.root.style.height = '180px';
  edited.setItems(items.slice(0, 6));
  main.append(edited.root);
  main.append(divH([
    button('Edit another cell', () => {
      changed.add('c4.formula');
      edited.refresh();
    }),
    button('Clear the refusal', () => {
      errors.clear();
      edited.refresh();
    }),
  ]));

  main.append(el('h2', null, 'Disposal'));
  const disposed = el('p', 'u2-gallery-status');
  main.append(button('Dispose', () => {
    table.dispose();
    edited.dispose();
    disposed.textContent = `Disposed: live scopes = ${Scope.liveCount}. Scrolling, keys and the ` +
      'context menu are dead — every listener and effect was owned by the component scope.';
  }));
  main.append(disposed);
  update();
}
