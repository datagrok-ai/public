import {signal, Scope, Control} from '../../src/index.js';
import {Tooltip} from '../../src/core/tooltip.js';

function injectOnce(id, make) {
  if (document.getElementById(id)) return;
  const e = make();
  e.id = id;
  document.head.append(e);
}

injectOnce('u2-tooltip-css', () => {
  const l = document.createElement('link');
  l.rel = 'stylesheet';
  l.href = new URL('../../css/tooltip.css', import.meta.url).href;
  return l;
});

const LONG_TEXT = 'Molecular weight is computed from the canonical structure: isotopes are ' +
  'resolved to their standard atomic weights, explicit and implicit hydrogens are counted, ' +
  'and salts are kept unless the structure was desalted upstream. Wraps at 400 px.';

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

function target(text) {
  const e = el('span', null, text);
  e.style.cssText = 'padding: var(--dg-space-s) var(--dg-space-m); border: var(--dg-border-width) ' +
    'dashed var(--dg-border); border-radius: var(--dg-radius); cursor: default;';
  return e;
}

export async function render(main) {
  main.append(el('h1', null, 'Tooltip'));
  const intro = el('p');
  intro.innerHTML = 'One singleton element in the shared <code>Overlay</code> layer, positioned ' +
    'with Floating UI. <code>Tooltip.bind(el, content, scope)</code> shows after a 300 ms hover ' +
    'delay; moving to another bound element within 100 ms switches instantly. Pointer-leave, ' +
    'pointer-down, scroll and Esc hide it. The element is <code>pointer-events: none</code>, so ' +
    'it never steals the hover, and it is reused — never disposed.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  const page = new Control();
  const now = signal(new Date().toLocaleTimeString());
  const timer = setInterval(() => now.value = new Date().toLocaleTimeString(), 1000);
  page.own(() => clearInterval(timer));

  main.append(el('h2', null, 'String, callback and long content'));
  const row = el('div');
  row.style.cssText = 'display: flex; gap: var(--dg-space-xl); align-items: center;';

  const plain = target('String tooltip');
  Tooltip.bind(plain, 'Molecular weight, g/mol', page.scope);

  const live = target('HTML callback');
  Tooltip.bind(live, () => {
    const box = el('div');
    box.append(el('b', null, 'Evaluated at show time'));
    box.append(el('div', null, `signal read on open: ${now.value}`));
    return box;
  }, page.scope);

  const long = target('Long text');
  Tooltip.bind(long, LONG_TEXT, page.scope);

  row.append(plain, live, long);
  main.append(row);
  main.append(el('p', 'u2-gallery-status',
    'Hover one, then slide across the others: the first wait is 300 ms, the switches are instant. ' +
    'The callback tooltip re-reads the clock signal every time it opens.'));

  main.append(el('h2', null, 'Disposal'));
  const note = el('p', 'u2-gallery-status',
    'Disposing the page component unbinds all three: hovering shows nothing, and live scopes drop. ' +
    'The tooltip element itself stays in the overlay host — it is shared and reused.');
  main.append(note, button('Dispose page component', () => {
    page.dispose();
    refresh();
  }));
  refresh();
}
