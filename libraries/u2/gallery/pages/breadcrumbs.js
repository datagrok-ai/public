import {computed, Scope, bindText} from '../../src/index.js';
import {Breadcrumbs} from '../../src/components/navigation/breadcrumbs.js';

function injectOnce(id, make) {
  if (document.getElementById(id)) return;
  const e = make();
  e.id = id;
  document.head.append(e);
}

injectOnce('u2-breadcrumbs-css', () => {
  const l = document.createElement('link');
  l.rel = 'stylesheet';
  l.href = new URL('../../css/breadcrumbs.css', import.meta.url).href;
  return l;
});

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function button(text, onClick) {
  const b = el('button', null, text);
  b.addEventListener('click', onClick);
  return b;
}

function frame(width) {
  const box = el('div');
  box.style.width = width;
  box.style.padding = 'var(--dg-space-s) var(--dg-space-m)';
  box.style.border = 'var(--dg-border-width) solid var(--dg-border)';
  box.style.boxSizing = 'border-box';
  return box;
}

const PATHS = {
  'Home / Projects / Chemistry / demog.csv': ['Home', 'Projects', 'Chemistry', 'demog.csv'],
  'Two segments': ['Home', 'Settings'],
  'One segment': ['Home'],
  'Long names': ['Home', 'Shared with everyone', 'Discovery chemistry 2026', 'sar-small-molecules.csv'],
};

const DEEP = ['Home', 'Files', 'Demo', 'Chem', 'Screening', 'Runs', '2026', 'August', 'plate-17.csv'];

export async function render(main) {
  main.append(el('h1', null, 'Breadcrumbs'));
  const intro = el('p');
  intro.innerHTML = 'Real buttons in a <code>nav</code>, so Tab reaches them and Enter or Space ' +
    'activates. The last segment is the current one: not clickable, ' +
    '<code>aria-current="page"</code>. When the row does not fit, the middle collapses into a ' +
    '<code>…</code> segment that expands the full path inline; a ResizeObserver re-measures, so ' +
    'widening the container brings the segments back.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  const log = el('div', 'u2-gallery-status');
  const logLine = (text) => log.append(el('div', null, text));

  const path = new Breadcrumbs({
    items: PATHS['Home / Projects / Chemistry / demog.csv'],
    onClick: (index, item) => logLine(`clicked [${index}] "${item}"`),
  });
  const pathFrame = frame('420px');
  pathFrame.append(path.root);
  main.append(pathFrame);

  const current = el('code');
  bindText(path.scope, current, computed(() => path.items.value.join(' / ')));
  const readout = el('p', 'u2-gallery-status');
  readout.append('items = ', current);
  main.append(readout);

  main.append(el('h2', null, 'setItems'));
  const setNote = el('p', 'u2-gallery-status',
    'setItems replaces the path and re-collapses it — an expanded "…" does not survive a new path.');
  const row = el('div');
  for (const [label, items] of Object.entries(PATHS))
    row.append(button(label, () => path.setItems(items)), ' ');
  main.append(setNote, row);

  main.append(el('h2', null, 'Overflow'));
  const deepNote = el('p', 'u2-gallery-status',
    'Nine segments in a narrow box: the middle collapses. Click "…" to expand it inline, or widen ' +
    'the box — the ResizeObserver re-measures and shows whatever fits.');
  const deep = new Breadcrumbs({
    items: DEEP,
    onClick: (index, item) => logLine(`deep: clicked [${index}] "${item}"`),
  });
  const deepFrame = frame('300px');
  deepFrame.append(deep.root);
  const widths = el('div');
  widths.style.paddingTop = 'var(--dg-space-m)';
  for (const width of ['240px', '300px', '460px', '760px'])
    widths.append(button(width, () => deepFrame.style.width = width), ' ');
  main.append(deepNote, deepFrame, widths);

  main.append(el('h2', null, 'Click log'));
  main.append(log);

  main.append(el('h2', null, 'Disposal'));
  const note = el('p', 'u2-gallery-status',
    'Dispose drops the render effect, the click listener and the ResizeObserver — live scopes drop back.');
  main.append(note, button('Dispose both', () => {
    path.dispose();
    deep.dispose();
    refresh();
  }));
  refresh();
}
