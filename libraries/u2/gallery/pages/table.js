import {signal, computed, Scope, Control} from '../../src/index.js';
import {div, divH, span, button} from '../../src/core/elements.js';
import {BasicTable, tableFromMap} from '../../src/components/collections/table.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-table-css', '../../css/table.css');

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

const STATUS_COLORS = {Approved: 'var(--dg-green-2)', Phase: 'var(--dg-orange-2)', Retired: 'var(--dg-grey-4)'};

function status(text) {
  const chip = span(text);
  chip.style.color = STATUS_COLORS[text];
  return chip;
}

const COMPOUNDS = [
  {name: 'Aspirin', formula: 'C9H8O4', mw: 180.16, status: 'Approved'},
  {name: 'Caffeine', formula: 'C8H10N4O2', mw: 194.19, status: 'Approved'},
  {name: 'Ibuprofen', formula: 'C13H18O2', mw: 206.28, status: 'Approved'},
  {name: 'Paracetamol', formula: 'C8H9NO2', mw: 151.16, status: 'Approved'},
  {name: 'Sildenafil', formula: 'C22H30N6O4S', mw: 474.58, status: 'Phase'},
  {name: 'Rofecoxib', formula: 'C17H14O4S', mw: 314.36, status: 'Retired'},
];

const COLUMNS = [
  {header: 'Name', render: (c) => c.name, width: '160px'},
  {header: 'Formula', render: (c) => c.formula, width: '140px'},
  {header: 'MW', render: (c) => c.mw.toFixed(2), align: 'right', width: '80px'},
  {header: 'Status', render: (c) => status(c.status), width: '90px'},
];

export async function render(main) {
  main.append(el('h1', null, 'Table'));
  const intro = el('p');
  intro.innerHTML = 'A real <code>&lt;table&gt;</code> for small data — the ' +
    '<code>ui.table</code>/<code>ui.tableFromMap</code> parity control. Cells are whatever the ' +
    'renderer returns (a string becomes text, an element is appended), and each re-render runs the ' +
    'renderers under a fresh child scope, so signal-bound cells are released instead of accumulating. ' +
    'Thousands of rows still belong in the virtual list.';
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

  const log = signal('nothing clicked yet');
  section('Compounds — selectable', () => {
    const table = new BasicTable({
      columns: COLUMNS,
      items: COMPOUNDS,
      selectable: true,
      onRowClick: (item, index) => log.value = `${item.name} (row ${index})`,
    });
    return [
      table,
      readout('selectedIndex', computed(() => String(table.selectedIndex.value))),
      readout('onRowClick', log),
      div([span('Click a row, or focus the table and use ↑ ↓ Home End, then Enter.')], 'u2-gallery-status'),
    ];
  });

  const rows = signal(COMPOUNDS.slice(0, 3));
  section('Live items', () => {
    const table = new BasicTable({columns: COLUMNS, items: rows});
    const buttons = divH([
      button('Add row', () => {
        const next = COMPOUNDS[rows.value.length % COMPOUNDS.length];
        rows.value = [...rows.value, {...next, name: `${next.name} ${rows.value.length + 1}`}];
      }),
      button('Shuffle', () => rows.value = rows.value.slice().sort(() => Math.random() - 0.5)),
      button('Reset', () => rows.value = COMPOUNDS.slice(0, 3)),
    ]);
    buttons.style.gap = 'var(--dg-space-m)';
    return [table, buttons, readout('rows', computed(() => String(rows.value.length)))];
  });

  const ticks = signal(0);
  section('Signal-bound cells', () => {
    const table = new BasicTable({
      columns: [
        {header: 'Compound', render: (c) => c.name},
        {header: 'Clicks', render: () => span(computed(() => String(ticks.value))), align: 'right'},
      ],
      items: COMPOUNDS.slice(0, 3),
    });
    const row = divH([
      button('Tick', () => ticks.value++),
      button('Re-render rows', () => table.setItems(COMPOUNDS.slice(0, 3))),
    ]);
    row.style.gap = 'var(--dg-space-m)';
    return [
      table,
      row,
      div([span('Every cell binds the same signal; re-rendering the rows releases the old bindings, ' +
        'so the live-scope count above stays put.')], 'u2-gallery-status'),
    ];
  });

  main.append(el('h2', null, 'tableFromMap'));
  main.append(tableFromMap({
    'Name': 'Aspirin',
    'Formula': 'C9H8O4',
    'Molecular weight': 180.16,
    'Status': status('Approved'),
    'Source': 'DrugBank',
  }));

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Every table above is owned by its section component: disposing the sections releases the row ' +
    'scopes, the effects and the listeners, and live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
