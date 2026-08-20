import {signal, computed, Scope, bindText} from '../../src/index.js';
import {Toolbar} from '../../src/components/toolbar.js';
import {TextInput} from '../../src/components/text-input.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-toolbar-css', '../../css/toolbar.css');
injectOnce('u2-menu-css', '../../css/menu.css');
injectOnce('u2-tooltip-css', '../../css/tooltip.css');
injectOnce('u2-inputs-css', '../../css/inputs.css');

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

export async function render(main) {
  main.append(el('h1', null, 'Toolbar'));
  const intro = el('p');
  intro.innerHTML = 'A single compact row of flat icon-style buttons: toggles reflect a ' +
    '<code>Signal&lt;boolean&gt;</code> through an owned effect, a dropdown button opens a fresh ' +
    '<code>Menu</code> on every click (so checked states are current), and an inline input rides ' +
    'along in an element slot. One tab stop — Tab enters the toolbar, ←/→ move between items, ' +
    'Enter or Space activates. Narrow the row below to watch whole items collapse into the » menu.';
  main.append(intro);

  const wrap = signal(true);
  const grid = signal(false);
  const last = signal('(none)');
  const page = new Scope();
  Scope.ambient?.own(() => page.dispose());

  const bar = new Toolbar({ariaLabel: 'Document actions'});
  const filter = new TextInput({inline: true, search: true, placeholder: 'Filter…'});
  filter.root.style.width = '130px';

  const act = (name) => last.value = name;
  bar
    .addButton('New', () => act('New'), {icon: '✚', tooltip: 'New document (Ctrl+N)'})
    .addButton('Open', () => act('Open'), {icon: '▤', tooltip: 'Open… (Ctrl+O)'})
    .addButton('Save', () => act('Save'), {icon: '⭳', tooltip: 'Save (Ctrl+S)'})
    .addSeparator()
    .addButton('Undo', () => act('Undo'), {icon: '↶', tooltip: 'Undo (Ctrl+Z)'})
    .addButton('Redo', () => act('Redo'), {icon: '↷', tooltip: 'Redo (Ctrl+Y)'})
    .addSeparator()
    .addToggle('Word wrap', wrap, {tooltip: 'Wrap long lines (Alt+Z)'})
    .addToggle('Grid', grid, {tooltip: 'Show the grid overlay'})
    .addSeparator()
    .addMenu('View', (menu) => menu
      .item('Word wrap', () => wrap.value = !wrap.peek(), {check: wrap.peek(), shortcut: 'Alt+Z'})
      .item('Grid', () => grid.value = !grid.peek(), {check: grid.peek()})
      .separator()
      .group('Zoom')
        .item('Zoom in', () => act('Zoom in'), {shortcut: 'Ctrl++'})
        .item('Zoom out', () => act('Zoom out'), {shortcut: 'Ctrl+-'})
      .endGroup()
      .item('Reset layout', () => act('Reset layout')))
    .addSeparator()
    .add(filter.root);

  const host = el('div');
  host.style.cssText = 'width: 720px; max-width: 100%; border: var(--dg-border-width) solid var(--dg-border);' +
    'border-radius: var(--dg-radius); resize: none; overflow: hidden;';
  host.append(bar.root);
  main.append(host);

  const action = el('b');
  const wrapText = el('code');
  const gridText = el('code');
  const filterText = el('code');
  bindText(page, action, last);
  bindText(page, wrapText, computed(() => String(wrap.value)));
  bindText(page, gridText, computed(() => String(grid.value)));
  bindText(page, filterText, computed(() => filter.value.value || '(empty)'));
  const readout = el('p', 'u2-gallery-status');
  readout.append('last action = ', action, '  ·  word wrap = ', wrapText, '  ·  grid = ', gridText,
    '  ·  filter = ', filterText);
  main.append(readout);

  main.append(el('h2', null, 'Overflow collapse'));
  const hint = el('p', 'u2-gallery-status',
    'Items are dropped whole, from the right, and a separator left dangling at the edge goes with ' +
    'them. The » menu mirrors what was dropped — buttons become items, toggles become check items ' +
    'bound to the same signals, the dropdown becomes a submenu; it is rebuilt on every open.');
  const slider = el('input');
  slider.type = 'range';
  slider.min = '150';
  slider.max = '760';
  slider.value = '720';
  slider.style.cssText = 'width: 320px; height: auto; border: none; vertical-align: middle;';
  const widthText = el('code', null, '720px');
  slider.addEventListener('input', () => {
    host.style.width = `${slider.value}px`;
    widthText.textContent = `${slider.value}px`;
  });
  const row = el('p');
  row.append('Container width: ', slider, ' ', widthText);
  main.append(hint, row);

  main.append(el('h2', null, 'Disposal'));
  const live = el('code', null, `${Scope.liveCount} live scopes`);
  const note = el('p', 'u2-gallery-status');
  note.append('Dispose drops the toggle effects, every button listener, the tooltip bindings, ' +
    'the open menu and the ResizeObserver: ', live, '.');
  main.append(note, button('Dispose toolbar', () => {
    bar.dispose();
    filter.dispose();
    page.dispose();
    live.textContent = `${Scope.liveCount} live scopes after disposal`;
  }));
}
