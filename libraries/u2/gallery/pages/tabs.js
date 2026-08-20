import {computed, Scope, bindText} from '../../src/index.js';
import {TabStrip} from '../../src/components/containers/tabs.js';

function injectOnce(id, make) {
  if (document.getElementById(id)) return;
  const e = make();
  e.id = id;
  document.head.append(e);
}

injectOnce('u2-tabs-css', () => {
  const l = document.createElement('link');
  l.rel = 'stylesheet';
  l.href = new URL('../../css/tabs.css', import.meta.url).href;
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

function paragraph(text) {
  const p = el('div');
  p.style.padding = '4px 0';
  p.textContent = text;
  return p;
}

const FILES = [
  'overview.md', 'package.ts', 'detectors.js', 'tokens.css', 'combobox.ts', 'list.ts',
  'scope.ts', 'signals.ts', 'queries.sql', 'demog.csv',
];

export async function render(main) {
  main.append(el('h1', null, 'Tabs'));
  const intro = el('p');
  intro.innerHTML = 'IDE document tabs: manual activation (focus moves with ←/→/Home/End, ' +
    'Enter or Space activates), Delete or middle-click closes a closable tab, panels are hidden ' +
    'with <code>display:none</code> so their state survives switching, and lazy content is built ' +
    'once on first activation. Twelve tabs force the overflow chevron at normal widths.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  const log = el('div', 'u2-gallery-status');
  const logLine = (text) => log.append(el('div', null, text));

  const lazy = () => {
    const started = performance.now();
    const box = el('div');
    box.append(paragraph('Built lazily on first activation — never rebuilt afterwards.'));
    for (let i = 0; i < 300; i++)
      box.append(paragraph(`row ${i}: ${'█'.repeat(1 + (i % 12))}`));
    logLine(`lazy.json content built in ${(performance.now() - started).toFixed(1)} ms ` +
      `at ${new Date().toLocaleTimeString()}`);
    return box;
  };

  const stateful = el('div');
  stateful.append(paragraph('Type below, switch tabs, come back — the value is still here ' +
    '(panels are hidden, not removed).'));
  const memo = el('input');
  memo.placeholder = 'Scratch…';
  stateful.append(memo);

  const defs = FILES.map((name) => ({id: name, label: name, content: paragraph(`Contents of ${name}`)}));
  defs.splice(1, 0, {id: 'draft.md', label: 'draft.md', closable: true, content: stateful});
  defs.push({id: 'lazy.json', label: 'lazy.json', closable: true, content: lazy});

  const tabs = new TabStrip({tabs: defs});
  tabs.root.style.height = '260px';
  tabs.root.style.maxWidth = '620px';
  tabs.root.style.border = 'var(--dg-border-width) solid var(--dg-border)';
  main.append(tabs.root);

  const active = el('b');
  const closed = el('code');
  bindText(tabs.scope, active, computed(() => tabs.activeTab.value ?? '(none)'));
  bindText(tabs.scope, closed, computed(() => tabs.onTabClosed.value ?? '(none)'));
  const readout = el('p', 'u2-gallery-status');
  readout.append('activeTab = ', active, '  ·  onTabClosed = ', closed);
  main.append(readout);

  main.append(el('h2', null, 'Dirty marker'));
  const hint = el('p', 'u2-gallery-status',
    'A ● dot sits where the ✕ is; hovering the tab swaps it back to ✕. A pinned (non-closable) ' +
    'tab shows the dot alone.');
  let draftDirty = false;
  let pinnedDirty = false;
  const row = el('div');
  row.append(
    button('Toggle dirty — draft.md (closable)', () => {
      draftDirty = !draftDirty;
      tabs.setDirty('draft.md', draftDirty);
    }),
    ' ',
    button('Toggle dirty — overview.md (pinned)', () => {
      pinnedDirty = !pinnedDirty;
      tabs.setDirty('overview.md', pinnedDirty);
    })
  );
  main.append(hint, row);

  main.append(el('h2', null, 'Lazy build log'));
  main.append(log);

  main.append(el('h2', null, 'Disposal'));
  const note = el('p', 'u2-gallery-status',
    'Dispose drops the active/menu effects, the header and document listeners and the ' +
    'ResizeObserver — live scopes drop back.');
  main.append(note, button('Dispose tab strip', () => {
    tabs.dispose();
    refresh();
  }));
  refresh();
}
