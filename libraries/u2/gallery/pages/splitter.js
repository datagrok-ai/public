import {Scope} from '../../src/index.js';
import {Splitter} from '../../src/components/splitter.js';

const OUTER_SIZES = [0.35, 0.65];
const INNER_SIZES = [0.5, 0.5];

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

function ensureCss() {
  if (document.getElementById('u2-splitter-css')) return;
  const link = el('link');
  link.id = 'u2-splitter-css';
  link.rel = 'stylesheet';
  link.href = new URL('../../css/splitter.css', import.meta.url).href;
  document.head.append(link);
}

function pane(title, text) {
  const p = el('div');
  p.style.padding = '8px 12px';
  p.style.height = '100%';
  p.style.boxSizing = 'border-box';
  p.append(el('b', null, title), el('p', 'u2-gallery-status', text));
  return p;
}

function fmt(sizes) {
  return sizes.map((s) => `${(s * 100).toFixed(1)}%`).join(' / ');
}

export async function render(main) {
  ensureCss();
  main.append(el('h1', null, 'Splitter'));
  const intro = el('p');
  intro.innerHTML = 'Panels laid out by flex-basis percentages driven off a <code>sizes</code> signal. ' +
    'The sash is a 4px hit area over a 1px line (VS Code\'s model): drag it with pointer capture, ' +
    'double-click to split the two neighbours evenly, or focus it and use the arrow keys ' +
    '(2% per press, 10% with Shift). The right panel is a second splitter — nesting is just a panel ' +
    'whose content happens to be another component.';
  main.append(intro);

  const inner = new Splitter([
    pane('Nested top', 'Vertical splitter inside the right panel of the horizontal one.'),
    pane('Nested bottom', 'Its sash and the outer one drag independently; both keep their own sizes signal.'),
  ], {direction: 'vertical', sizes: INNER_SIZES});

  const outer = new Splitter([
    pane('Left', 'Drag the sash on the right edge. Panels never shrink below minSize (60px by default).'),
    inner.root,
  ], {direction: 'horizontal', sizes: OUTER_SIZES});

  const frame = el('div');
  frame.style.height = '360px';
  frame.style.border = 'var(--dg-border-width) solid var(--dg-border)';
  frame.append(outer.root);
  main.append(frame);

  const readout = el('p', 'u2-gallery-status');
  outer.effect(() => readout.textContent =
    `outer = ${fmt(outer.sizes.value)} · nested = ${fmt(inner.sizes.value)} · live scopes = ${Scope.liveCount}`);
  main.append(readout);

  const disposed = el('p', 'u2-gallery-status');
  const row = el('div');
  row.append(
    button('Reset', () => {
      outer.sizes.value = OUTER_SIZES.slice();
      inner.sizes.value = INNER_SIZES.slice();
    }),
    ' ',
    button('Dispose', () => {
      outer.dispose();
      inner.dispose();
      disposed.textContent = `Disposed: live scopes = ${Scope.liveCount}. ` +
        'Both splitters are inert — every sash listener, hover timer and layout effect belonged to ' +
        'the component scopes.';
    }));
  main.append(row, disposed);
}
