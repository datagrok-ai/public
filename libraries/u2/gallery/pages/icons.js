import {signal, computed, Scope, Component} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {icon} from '../../src/components/icon.js';

// The gallery is a standalone host: without this marker css/icons.css contributes no font at all.
document.documentElement.classList.add('u2-standalone');

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-tooltip-css', '../../css/tooltip.css');
injectOnce('u2-icons-css', '../../css/icons.css');

if (!document.getElementById('u2-icons-page-css')) {
  const style = document.createElement('style');
  style.id = 'u2-icons-page-css';
  style.textContent = `
    .u2-icon-grid {
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(104px, 1fr));
      gap: var(--dg-space-s);
    }
    .u2-icon-cell {
      display: flex;
      flex-direction: column;
      align-items: center;
      gap: var(--dg-space-s);
      padding: var(--dg-space-m) var(--dg-space-s);
      border: var(--dg-border-width) solid var(--dg-border-subtle);
      border-radius: var(--dg-radius);
      overflow: hidden;
    }
    .u2-icon-cell .name {
      max-width: 100%;
      font-size: var(--dg-font-size-xs);
      color: var(--dg-text-color-light);
      text-overflow: ellipsis;
      overflow: hidden;
      white-space: nowrap;
    }`;
  document.head.append(style);
}

const NAMES = [
  'plus', 'minus', 'times', 'check', 'search', 'filter', 'sync', 'undo', 'redo', 'clone',
  'save', 'edit', 'pen', 'trash-alt', 'copy', 'paste', 'cut', 'tag', 'tags',
  'folder', 'folder-open', 'file', 'file-alt', 'file-import', 'file-export',
  'download', 'upload', 'print', 'share-alt', 'link', 'external-link-alt',
  'chevron-up', 'chevron-down', 'chevron-left', 'chevron-right',
  'angle-double-left', 'angle-double-right',
  'arrow-up', 'arrow-down', 'arrow-left', 'arrow-right',
  'sort', 'sort-up', 'sort-down', 'bars', 'ellipsis-h', 'ellipsis-v',
  'cog', 'cogs', 'wrench', 'sliders-h', 'tools', 'key', 'lock', 'unlock',
  'user', 'users', 'user-plus', 'sign-in-alt', 'sign-out-alt',
  'home', 'star', 'heart', 'bookmark', 'flag', 'bell', 'envelope', 'comment', 'comments',
  'question-circle', 'info-circle', 'exclamation-triangle', 'exclamation-circle',
  'check-circle', 'times-circle', 'ban', 'clock', 'history', 'calendar', 'calendar-alt',
  'table', 'th', 'th-list', 'list', 'columns',
  'chart-bar', 'chart-line', 'chart-pie', 'chart-area', 'project-diagram',
  'database', 'server', 'code', 'terminal', 'play', 'pause', 'stop', 'spinner',
  'eye', 'eye-slash', 'expand', 'compress',
  'flask', 'atom', 'microscope', 'dna', 'vial',
];

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

function row(children) {
  const r = divH(children);
  r.style.gap = 'var(--dg-space-l)';
  r.style.alignItems = 'center';
  r.style.flexWrap = 'wrap';
  return r;
}

function labelled(text, node) {
  const cell = el('div');
  cell.style.cssText = 'display: flex; flex-direction: column; align-items: center; gap: 4px;';
  cell.append(node, el('div', 'u2-gallery-status', text));
  return cell;
}

export async function render(main) {
  main.append(el('h1', null, 'Icons'));
  const intro = el('p');
  intro.innerHTML = '<code>icon(name, options)</code> builds the platform\'s own markup — ' +
    '<code>&lt;i class="grok-icon u2-icon fal fa-{name}"&gt;</code> — so in the app it renders ' +
    'with the Font Awesome 5 Pro sheet Datagrok already loads, untouched. Standalone hosts like ' +
    'this gallery set <code>u2-standalone</code> on <code>&lt;html&gt;</code> and load ' +
    '<code>css/icons.css</code>, which maps the same classes onto the vendored Font Awesome ' +
    'Pro 5.15.4 webfonts (copied from the platform\'s own assets) — so <code>fal</code> is the ' +
    'true Light face here too, and light/regular/solid genuinely differ.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  const parts = [];
  const section = (title, builder) => {
    main.append(el('h2', null, title));
    const component = Component.build(builder);
    parts.push(component);
    main.append(component.root);
    return component;
  };

  const query = signal('');
  const matches = computed(() => NAMES.filter((n) => n.includes(query.value.trim().toLowerCase())));
  section(`Browse (${NAMES.length} free names)`, () => {
    const search = el('input');
    search.type = 'search';
    search.name = 'icon-filter';
    search.setAttribute('aria-label', 'Filter icons by name');
    search.placeholder = 'Filter icons…';
    search.style.width = '240px';
    search.addEventListener('input', () => query.value = search.value);

    const grid = el('div', 'u2-icon-grid');
    const cells = new Map();
    for (const name of NAMES) {
      const cell = el('div', 'u2-icon-cell');
      cell.append(icon(name, {size: 'large', tooltip: name}), el('div', 'name', name));
      cells.set(name, cell);
      grid.append(cell);
    }
    const empty = el('div', 'u2-gallery-status', 'No icon matches that name.');
    Scope.ambient.effect(() => {
      const visible = new Set(matches.value);
      for (const [name, cell] of cells)
        cell.style.display = visible.has(name) ? '' : 'none';
      empty.style.display = visible.size === 0 ? '' : 'none';
    });
    return [
      row([search, readout('showing', computed(() => `${matches.value.length} of ${NAMES.length}`))]),
      grid,
      empty,
    ];
  });

  section('Variants', () => [
    row([
      labelled('light (fal)', icon('star', {size: 'large'})),
      labelled('regular (far)', icon('star', {size: 'large', variant: 'regular'})),
      labelled('solid (fas)', icon('star', {size: 'large', variant: 'solid'})),
    ]),
    el('div', 'u2-gallery-status',
      'Three genuine Pro weights — the same faces in the app and standalone.'),
  ]);

  section('Sizes', () => [
    row([
      labelled('small — 12px', icon('cog', {size: 'small'})),
      labelled('medium — 16px', icon('cog', {size: 'medium'})),
      labelled('large — 20px', icon('cog', {size: 'large'})),
      labelled('inherit — no size class', icon('cog')),
    ]),
    el('div', 'u2-gallery-status',
      'Without a size the icon takes the font size of whatever it sits in.'),
  ]);

  section('Color', () => {
    const inherit = (color, fontSize) => {
      const host = el('div');
      host.style.cssText = `color: ${color}; font-size: ${fontSize};`;
      host.append(icon('flask'), ' inherited from the parent');
      return host;
    };
    return [
      inherit('var(--dg-text-color)', 'var(--dg-font-size)'),
      inherit('var(--dg-accent)', 'var(--dg-font-size-h2)'),
      inherit('var(--dg-red-3)', 'var(--dg-font-size-small)'),
      el('div', 'u2-gallery-status',
        'u2-icon sets color: currentColor and no size, so an icon takes both from its context.'),
    ];
  });

  const clicks = signal(0);
  section('Tooltip and actions', () => [
    row([
      icon('info-circle', {size: 'large', tooltip: 'Bound through the ambient Scope — hover me'}),
      icon('sync', {size: 'large', cls: 'u2-icon-action', tooltip: 'Refresh'}),
      icon('trash-alt', {size: 'large', cls: 'u2-icon-action'}),
      button('Count a click', () => clicks.value++),
    ]),
    readout('clicks', computed(() => String(clicks.value))),
    el('div', 'u2-gallery-status',
      'u2-icon-action adds the pointer cursor and the muted→accent hover for clickable icons. ' +
      'A tooltip binds to the ambient Scope when there is one, and falls back to the title ' +
      'attribute when the icon is built outside a component.'),
  ]);

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Every section is a Component.build(...) builder, so the grid filter effect and every ' +
    'tooltip binding are owned by it: disposing the sections stops the filter and releases the ' +
    'hover listeners, and live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
