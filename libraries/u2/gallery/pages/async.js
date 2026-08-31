import {AsyncSource, Scope} from '../../src/index.js';
import {AsyncView, loader, skeleton} from '../../src/components/display/async-view.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-async-css', '../../css/async.css');

const COMPOUNDS = [
  'Acetic acid', 'Acetone', 'Adenine', 'Aniline', 'Ascorbic acid', 'Benzene', 'Caffeine',
  'Cholesterol', 'Citric acid', 'Dopamine', 'Ethanol', 'Glucose', 'Glycine', 'Ibuprofen',
  'Morphine', 'Nicotine', 'Paracetamol', 'Phenol', 'Salicylic acid', 'Serotonin', 'Urea',
];

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

function fakeFetch(query, abort) {
  return new Promise((resolve, reject) => {
    const timer = setTimeout(() => {
      const q = query.trim().toLowerCase();
      if (q === 'err') {
        reject(new Error('Compound service unavailable (HTTP 503)'));
        return;
      }
      if (q === 'none') {
        resolve([]);
        return;
      }
      resolve(q ? COMPOUNDS.filter((c) => c.toLowerCase().includes(q)) : COMPOUNDS);
    }, 600);
    abort.addEventListener('abort', () => {
      clearTimeout(timer);
      reject(new DOMException('Aborted', 'AbortError'));
    });
  });
}

function renderList(items) {
  const list = el('div');
  list.style.display = 'flex';
  list.style.flexWrap = 'wrap';
  list.style.gap = 'var(--dg-space-m)';
  for (const item of items)
    list.append(el('span', 'u2-gallery-status', item));
  return list;
}

function frame(width) {
  const box = el('div');
  box.style.width = width;
  box.style.minHeight = '90px';
  box.style.padding = 'var(--dg-space-m)';
  box.style.border = 'var(--dg-border-width) solid var(--dg-border)';
  box.style.borderRadius = 'var(--dg-radius)';
  return box;
}

export async function render(main) {
  main.append(el('h1', null, 'Async states'));
  const intro = el('p');
  intro.innerHTML = 'One <code>AsyncSource</code>, one <code>AsyncView</code>: idle, loading ' +
    '(spinner or skeleton), ready, empty and error-with-Retry, swapped atomically. Controls ' +
    'stop hand-rolling spinners — <code>loader()</code> and <code>skeleton()</code> are also ' +
    'usable standalone. Type <code>err</code> to fail the provider, <code>none</code> to get ' +
    'no rows; the provider takes 600 ms and every new query aborts the in-flight one.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refreshCount = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  main.append(el('h2', null, 'Own source — query on Enter'));
  const search = el('input');
  search.placeholder = 'Compound, "err" or "none"…';
  search.style.width = '260px';
  const owned = AsyncView.owned(fakeFetch, renderList, {empty: 'No compounds match'});
  search.addEventListener('keydown', (e) => {
    if (e.key === 'Enter') owned.refresh(search.value);
  });
  const ownedBox = frame('520px');
  ownedBox.append(owned.root);
  main.append(search, ownedBox);

  main.append(el('h2', null, 'Shared source — spinner vs skeleton'));
  main.append(el('p', 'u2-gallery-status',
    'Both views render the same AsyncSource; neither owns it, so the page disposes it.'));
  const shared = new AsyncSource(fakeFetch);
  const spinnerView = new AsyncView(shared, renderList, {empty: 'No compounds'});
  const skeletonView = new AsyncView(shared, renderList, {empty: 'No compounds', skeleton: true});
  const spinnerBox = frame('auto');
  const skeletonBox = frame('auto');
  spinnerBox.style.flex = '1';
  skeletonBox.style.flex = '1';
  spinnerBox.append(el('div', 'u2-gallery-status', 'spinner'), spinnerView.root);
  skeletonBox.append(el('div', 'u2-gallery-status', 'skeleton'), skeletonView.root);
  const pair = el('div');
  pair.style.display = 'flex';
  pair.style.gap = 'var(--dg-space-l)';
  pair.append(spinnerBox, skeletonBox);
  const controls = el('div');
  controls.style.display = 'flex';
  controls.style.gap = 'var(--dg-space-m)';
  controls.style.margin = '0 0 var(--dg-space-m)';
  controls.append(
    button('Load', () => shared.query('')),
    button('Load empty', () => shared.query('none')),
    button('Fail', () => shared.query('err'))
  );
  main.append(controls, pair);

  main.append(el('h2', null, 'Specimens'));
  const specimens = el('div');
  specimens.style.display = 'flex';
  specimens.style.alignItems = 'flex-start';
  specimens.style.gap = 'var(--dg-space-xxl)';
  const spinners = el('div');
  spinners.style.display = 'flex';
  spinners.style.flexDirection = 'column';
  spinners.style.gap = 'var(--dg-space-l)';
  spinners.append(loader(), loader('Loading compounds…'));
  const skeletons = frame('320px');
  skeletons.append(skeleton(), skeleton(5));
  specimens.append(spinners, skeletons);
  main.append(specimens);

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Disposing a view releases its state effect and the scope of whatever it last rendered; ' +
    'the owned view also disposes its source, the shared one is disposed here by hand.'));
  main.append(button('Dispose views and shared source', () => {
    owned.dispose();
    spinnerView.dispose();
    skeletonView.dispose();
    shared.dispose();
    refreshCount();
  }));

  owned.refresh('');
  shared.query('');
  refreshCount();
}
