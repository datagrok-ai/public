import {Scope, signal} from '../../src/index.js';
import {AsyncPager} from '../../src/core/pager.js';
import {VirtualList} from '../../src/components/list.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-list-css', '../../css/list.css');
injectOnce('u2-async-css', '../../css/async.css');

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

// The "server": 120 reports, served in pages of 20 with a latency no user would call instant.
const KINDS = ['bug', 'idea', 'question'];
const REPORTS = Array.from({length: 120}, (_, i) => ({
  id: i,
  title: `Report #${String(i).padStart(3, '0')} — ${KINDS[i % 3]} in the ${i % 2 ? 'grid' : 'scatter plot'}`,
  kind: KINDS[i % 3],
}));
const PAGE_SIZE = 20;
const ITEM_HEIGHT = 26;
const LATENCY_MS = 300;

const delay = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

function matching(query) {
  const q = query.trim().toLowerCase();
  return q ? REPORTS.filter((r) => r.title.toLowerCase().includes(q)) : REPORTS;
}

export function render(main) {
  main.append(el('h1', null, 'Async pager'));
  const intro = el('p');
  intro.innerHTML = '<code>AsyncPager&lt;T&gt;</code> is the paged twin of <code>AsyncSource</code>: ' +
    'pages accumulate into one array signal a <code>VirtualList</code> renders directly, a page ' +
    'shorter than <code>pageSize</code> marks the collection <code>done</code>, and a ' +
    '<code>reset()</code> invalidates whatever is in flight before re-reading the total. The ' +
    'filter below is a <em>thunk</em>, re-read on every reset — which is how one pager serves a ' +
    'query box for its whole lifetime. This page pages a mock array; in the platform, ' +
    '<code>u2/dg</code>\'s <code>dapiPager(() =&gt; grok.dapi.reports)</code> pages a server ' +
    'collection through the same seam.';
  main.append(intro);

  const query = signal('');
  const pager = new AsyncPager({
    pageSize: PAGE_SIZE,
    count: async () => {
      await delay(LATENCY_MS);
      return matching(query.value).length;
    },
    fetchPage: async (page) => {
      await delay(LATENCY_MS);
      return matching(query.value).slice(page * PAGE_SIZE, (page + 1) * PAGE_SIZE);
    },
  });
  Scope.ambient?.own(() => pager.dispose());

  const search = el('input');
  search.placeholder = 'Filter, then watch the pager restart…';
  search.style.width = '320px';
  search.addEventListener('input', () => {
    query.value = search.value;
    pager.reset();
  });

  const controls = el('p');
  controls.append(search, ' ', button('Reset', () => pager.reset()), ' ',
    button('Load more', () => pager.loadMore()));
  main.append(controls);

  const list = new VirtualList({
    itemHeight: ITEM_HEIGHT,
    keyOf: (r) => String(r.id),
    render: (r) => el('span', null, `${r.title}  ·  ${r.kind}`),
  });
  list.root.style.height = '320px';
  list.setItems(pager.items);
  const onScroll = () => {
    const root = list.root;
    if (root.scrollTop + root.clientHeight > root.scrollHeight - 5 * ITEM_HEIGHT)
      pager.loadMore();
  };
  list.root.addEventListener('scroll', onScroll);
  list.own(() => list.root.removeEventListener('scroll', onScroll));
  main.append(list.root);

  const readout = el('p', 'u2-gallery-status');
  list.effect(() => readout.textContent =
    `state = ${pager.state.value} · loaded = ${pager.items.value.length} · ` +
    `total = ${pager.total.value ?? '…'} · pages of ${PAGE_SIZE}, ${LATENCY_MS} ms each`);
  main.append(readout);
  main.append(el('p', 'u2-gallery-status',
    'Scroll to within five rows of the end and the next page is requested; the list keeps its ' +
    'scroll position because the items signal only ever grows.'));

  pager.reset();
}
