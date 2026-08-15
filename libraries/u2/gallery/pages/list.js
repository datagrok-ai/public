import {Scope} from '../../src/index.js';
import {VirtualList} from '../../src/components/list.js';

const ITEM_COUNT = 1000000;

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
  if (document.getElementById('u2-list-css')) return;
  const link = el('link');
  link.id = 'u2-list-css';
  link.rel = 'stylesheet';
  link.href = new URL('../../css/list.css', import.meta.url).href;
  document.head.append(link);
}

export async function render(main) {
  ensureCss();
  main.append(el('h1', null, 'Virtual list'));
  const intro = el('p');
  intro.innerHTML = 'A port of VS Code\'s <code>listView</code> virtualization onto the u2 core: ' +
    'only the visible rows plus a small overscan exist in the DOM, row elements are recycled, and ' +
    'selection is a signal the row highlight is an effect of. Click a row, or focus the list and ' +
    'use ↑ ↓ Home End PageUp PageDown.';
  main.append(intro);

  const items = Array.from({length: ITEM_COUNT}, (_, i) => i);
  const list = new VirtualList({
    render: (item) => el('span', null, `Row ${item.toLocaleString()} — generated, never all in the DOM`),
  });
  list.root.style.height = '360px';
  list.setItems(items);
  main.append(list.root);

  const readout = el('p', 'u2-gallery-status');
  const update = () => readout.textContent =
    `${ITEM_COUNT.toLocaleString()} items · selectedIndex = ${list.selectedIndex.value}` +
    ` · rows in the DOM = ${list.renderedCount} · live scopes = ${Scope.liveCount}`;
  list.effect(update);
  list.root.addEventListener('scroll', update);
  list.own(() => list.root.removeEventListener('scroll', update));
  main.append(readout);

  const disposed = el('p', 'u2-gallery-status');
  const row = el('div');
  row.append(
    button('Scroll to 500,000', () => {
      list.scrollToIndex(500000);
      update();
    }),
    ' ',
    button('Dispose', () => {
      list.dispose();
      disposed.textContent = `Disposed: live scopes = ${Scope.liveCount}. ` +
        'Scrolling and keys are dead — every listener and effect was owned by the component scope.';
    }));
  main.append(row, disposed);
}
