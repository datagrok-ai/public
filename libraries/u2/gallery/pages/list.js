import {Scope, timestamp} from '../../src/index.js';
import {VirtualList} from '../../src/components/collections/list.js';
import {rowActions} from '../../src/components/actions/actions.js';

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
  for (const name of ['list', 'elements', 'buttons', 'icons', 'menu', 'adaptive']) {
    if (document.getElementById(`u2-${name}-css`)) continue;
    const link = el('link');
    link.id = `u2-${name}-css`;
    link.rel = 'stylesheet';
    link.href = new URL(`../../css/${name}.css`, import.meta.url).href;
    document.head.append(link);
  }
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

  main.append(el('h2', null, 'Rich rows'));
  const richIntro = el('p');
  richIntro.innerHTML = 'The list-item-rendering recipe as primitives: <code>rowActions</code> ' +
    'reveals icon shortcuts on hover or focus, right-click opens the FULL action list ' +
    '(<code>contextActions</code>), <code>timestamp</code> shows the short date with the full ' +
    'date-time on hover, and <code>u2-adaptive</code> + <code>u2-p2</code> hide the assignee ' +
    'as the pane narrows — drag the resize handle of your browser window to see it degrade.';
  main.append(richIntro);

  const log = el('p', 'u2-gallery-status');
  const day = 86400000;
  const issues = Array.from({length: 200}, (_, i) => ({
    id: i,
    title: `Issue #${i + 1}: something worth a two-line row`,
    assignee: ['ada', 'grace', 'edsger'][i % 3],
    created: new Date(Date.now() - i * 11 * day),
  }));
  const actionsOf = (issue) => [
    {name: 'Resolve', icon: 'check', run: () => log.textContent = `Resolved #${issue.id + 1}`},
    {name: 'Copy title', icon: 'copy', run: () => log.textContent = `Copied "${issue.title}"`},
    {name: 'Assign to me', run: () => log.textContent = `Assigned #${issue.id + 1} to you`},
  ];
  const rich = new VirtualList({
    itemHeight: 42,
    keyOf: (issue) => String(issue.id),
    contextActions: actionsOf,
    render: (issue) => {
      const root = el('div');
      root.style.cssText = 'display:flex;flex-direction:column;justify-content:center;flex:1;min-width:0;gap:2px';
      const top = el('div');
      top.style.cssText = 'display:flex;align-items:baseline;gap:8px;min-width:0';
      const title = el('span', null, issue.title);
      title.style.cssText = 'overflow:hidden;text-overflow:ellipsis;white-space:nowrap;min-width:0';
      const assignee = el('span', 'u2-p2', `@${issue.assignee}`);
      assignee.style.cssText = 'margin-left:auto;color:var(--dg-text-color-light);flex-shrink:0';
      top.append(title, assignee);
      const bottom = el('div');
      bottom.style.cssText = 'display:flex;align-items:center;gap:8px';
      bottom.append(timestamp(issue.created), rowActions(actionsOf(issue)));
      root.append(top, bottom);
      return root;
    },
  });
  rich.root.classList.add('u2-adaptive');
  rich.root.style.height = '240px';
  rich.setItems(issues);
  main.append(rich.root, log);
}
