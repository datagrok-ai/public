import {signal, computed, Scope, Control, VirtualGrid, IconInput, ICON_NAMES, BRAND_ICON_NAMES} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {icon} from '../../src/components/display/icon.js';

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
injectOnce('u2-inputs-css', '../../css/inputs.css');
injectOnce('u2-grid-css', '../../css/grid.css');
injectOnce('u2-icon-input-css', '../../css/icon-input.css');

if (!document.getElementById('u2-icons-page-css')) {
  const style = document.createElement('style');
  style.id = 'u2-icons-page-css';
  style.textContent = `
    .u2-icon-grid {
      height: 360px;
    }
    .u2-icon-cell {
      display: flex;
      flex-direction: column;
      align-items: center;
      gap: var(--dg-space-s);
      width: 100%;
      padding: var(--dg-space-s);
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

const NAMES = [...ICON_NAMES, ...BRAND_ICON_NAMES];

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
    const component = Control.build(builder);
    parts.push(component);
    main.append(component.root);
    return component;
  };

  const query = signal('');
  const matches = computed(() => {
    const q = query.value.trim().toLowerCase().replace(/\s+/g, '-');
    return q ? NAMES.filter((n) => n.includes(q)) : NAMES;
  });
  section(`Browse (${ICON_NAMES.length} icons + ${BRAND_ICON_NAMES.length} brands, virtualized)`, () => {
    const search = el('input');
    search.type = 'search';
    search.name = 'icon-filter';
    search.setAttribute('aria-label', 'Filter icons by name');
    search.placeholder = 'Filter icons…';
    search.style.width = '240px';
    search.addEventListener('input', () => query.value = search.value);

    const grid = new VirtualGrid({cellWidth: 104, cellHeight: 56, render: (name) => {
      const cell = el('div', 'u2-icon-cell');
      cell.append(icon(name, {size: 'large', tooltip: name}), el('div', 'name', name));
      return cell;
    }});
    grid.root.classList.add('u2-icon-grid');
    grid.setItems(matches);
    const empty = el('div', 'u2-gallery-status', 'No icon matches that name.');
    const inDom = signal(0);
    const count = () => inDom.value = grid.renderedCount;
    grid.root.addEventListener('scroll', count);
    Scope.ambient.effect(() => empty.style.display = matches.value.length === 0 ? '' : 'none');
    Scope.ambient.effect(() => { matches.value; queueMicrotask(count); });
    return [
      row([search, readout('showing', computed(() => `${matches.value.length} of ${NAMES.length}`)),
        readout('in DOM', computed(() => String(inDom.value)))]),
      grid,
      empty,
    ];
  });

  section('IconInput', () => {
    const picker = new IconInput({label: 'Icon', value: 'flask'});
    const brand = new IconInput({label: 'Brand', value: 'github'});
    const inline = new IconInput({inline: true, names: ['star', 'heart', 'flag', 'bookmark']});
    return [
      row([picker, brand, inline]),
      readout('picked', computed(() => `${picker.value.value} · ${brand.value.value} · ${inline.value.value || '(none)'}`)),
      el('div', 'u2-gallery-status',
        'A click (or Enter) opens the grid; type to filter; click or Enter picks, Esc cancels, Backspace clears. ' +
        'A brand name renders in its own face without being told.'),
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
    'Every section is a Control.build(...) builder, so the grid filter effect and every ' +
    'tooltip binding are owned by it: disposing the sections stops the filter and releases the ' +
    'hover listeners, and live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
