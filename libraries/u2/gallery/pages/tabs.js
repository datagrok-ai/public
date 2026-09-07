import {computed, Scope, bindText} from '../../src/index.js';
import {TabStrip} from '../../src/components/containers/tabs.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  l.id = id;
  document.head.append(l);
}

injectOnce('u2-icons-css', '../../css/icons.css');
injectOnce('u2-tabs-css', '../../css/tabs.css');

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function html(tag, cls, markup) {
  const e = el(tag, cls);
  e.innerHTML = markup;
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
  intro.innerHTML = 'The platform\'s <code>ui.tabControl</code> look by default: manual activation ' +
    '(focus moves with ←/→/Home/End, Enter or Space activates), Delete or middle-click closes a ' +
    'closable tab, panels are hidden with <code>display:none</code> so their state survives ' +
    'switching, and lazy content is built once on first activation. Twelve tabs force the ' +
    'overflow chevron at normal widths.';
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

  main.append(el('h2', null, 'Document variant'));
  main.append(html('p', 'u2-gallery-status',
    '<code>variant: \'document\'</code> — the IDE editor-group skin: secondary stripe, accent top border, ' +
    'active tab on the page background.'));
  const documents = new TabStrip({variant: 'document', tabs: [
    {id: 'package.ts', label: 'package.ts', closable: true, content: paragraph('Contents of package.ts')},
    {id: 'README.md', label: 'README.md', closable: true, content: paragraph('Contents of README.md')},
    {id: 'tokens.css', label: 'tokens.css', content: paragraph('Contents of tokens.css')},
  ]});
  documents.root.style.height = '120px';
  documents.root.style.maxWidth = '620px';
  documents.root.style.border = 'var(--dg-border-width) solid var(--dg-border)';
  main.append(documents.root);

  main.append(el('h2', null, 'Vertical, with icons'));
  main.append(html('p', 'u2-gallery-status',
    '<code>orientation: \'vertical\'</code> — the header is a column on the left and ↑/↓ move focus; ' +
    '<code>icon</code> is a Font Awesome name or an element, and a tab with an empty label is ' +
    'icon-only (its <code>tooltip</code> carries the name).'));
  const vertical = new TabStrip({orientation: 'vertical', tabs: [
    {id: 'home', label: 'Home', icon: 'home', content: paragraph('Home')},
    {id: 'settings', label: 'Settings', icon: 'cog', content: paragraph('Settings')},
    {id: 'charts', label: 'Charts', icon: 'chart-line', content: paragraph('Charts')},
    {id: 'help', label: '', icon: 'question-circle', tooltip: 'Help', content: paragraph('Help (icon-only tab)')},
  ]});
  vertical.root.style.height = '180px';
  vertical.root.style.maxWidth = '620px';
  vertical.root.style.border = 'var(--dg-border-width) solid var(--dg-border)';
  main.append(vertical.root);

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
