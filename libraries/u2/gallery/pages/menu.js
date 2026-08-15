import {signal, computed, Scope, bindText} from '../../src/index.js';
import {Menu} from '../../src/components/menu.js';

function injectOnce(id, make) {
  if (document.getElementById(id)) return;
  const e = make();
  e.id = id;
  document.head.append(e);
}

injectOnce('u2-menu-css', () => {
  const l = document.createElement('link');
  l.rel = 'stylesheet';
  l.href = new URL('../../css/menu.css', import.meta.url).href;
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

function readout(scope, menu) {
  const state = el('code');
  const open = el('code');
  bindText(scope, state, menu.machineState);
  bindText(scope, open, computed(() => String(menu.isOpen.value)));
  const row = el('p', 'u2-gallery-status');
  row.append('machine = ', state, '  ·  isOpen = ', open);
  return row;
}

export async function render(main) {
  main.append(el('h1', null, 'Menu'));
  const intro = el('p');
  intro.innerHTML = 'Hand-written menu machine (<code>closed / open / submenu-open</code>, exposed as a ' +
    'signal), rendered into the shared <code>Overlay</code> layer — anchored via <code>Overlay.show</code>, ' +
    'at the cursor via a fixed-strategy panel with the same dismissal (outside pointerdown, Esc). ' +
    'Keys: ↑↓ move (wrapping, skipping separators and disabled items), → or Enter opens a submenu, ' +
    '← closes back to the parent, Esc closes the chain, a letter jumps to the next item starting with it.';
  main.append(intro);

  const page = new Scope();
  const wrap = signal(true);
  const last = signal('(none)');

  const action = el('b');
  const wrapText = el('b');
  bindText(page, action, last);
  bindText(page, wrapText, computed(() => String(wrap.value)));
  const values = el('p', 'u2-gallery-status');
  values.append('last action = ', action, '  ·  word wrap signal = ', wrapText);

  main.append(el('h2', null, 'Anchored — shortcuts, a checked toggle, a disabled item, nested submenus'));

  function build() {
    return Menu.popup()
      .item('Cut', () => last.value = 'Cut', {icon: '✂', shortcut: 'Ctrl+X'})
      .item('Copy', () => last.value = 'Copy', {shortcut: 'Ctrl+C'})
      .item('Paste', () => last.value = 'Paste', {shortcut: 'Ctrl+V', enabled: false})
      .separator()
      .item('Word wrap', () => wrap.value = !wrap.value, {check: wrap.value, shortcut: 'Alt+Z'})
      .separator()
      .group('Export')
        .item('CSV', () => last.value = 'Export → CSV')
        .item('Excel', () => last.value = 'Export → Excel')
        .group('Image')
          .item('PNG', () => last.value = 'Export → Image → PNG')
          .item('SVG', () => last.value = 'Export → Image → SVG')
        .endGroup()
      .endGroup()
      .item('Rename…', () => last.value = 'Rename', {shortcut: 'F2'});
  }

  const slot = el('div');
  const live = el('code', null, `${Scope.liveCount} live scopes`);
  let watch = null;
  let anchored = null;
  const open = button('Edit ▾', () => {
    if (anchored) anchored.close();
    if (watch) watch.dispose();
    watch = new Scope();
    anchored = build();
    slot.replaceChildren(readout(watch, anchored));
    bindText(watch, live, computed(() =>
      `${Scope.liveCount} live scopes while the menu is ${anchored.isOpen.value ? 'open' : 'closed'}`));
    anchored.show({anchor: open});
  });
  main.append(open, slot, values);

  main.append(el('h2', null, 'Context — right-click at the cursor'));
  const region = el('div', null, 'Right-click anywhere in this box');
  region.style.cssText = 'display: flex; align-items: center; justify-content: center; height: 120px;' +
    'border: var(--dg-border-width) dashed var(--dg-border); border-radius: var(--dg-radius);' +
    'color: var(--dg-text-color-light);';
  const context = Menu.popup()
    .item('Open', () => last.value = 'Open', {shortcut: 'Enter'})
    .item('Pin', () => last.value = 'Pin', {check: false})
    .separator()
    .group('Share')
      .item('Copy link', () => last.value = 'Share → Copy link')
      .item('Email…', () => last.value = 'Share → Email')
    .endGroup()
    .item('Delete', () => last.value = 'Delete', {shortcut: 'Del'});
  context.bindContext(region, page);
  main.append(region, readout(page, context));

  main.append(el('h2', null, 'Disposal'));
  const note = el('p', 'u2-gallery-status');
  note.append('Every show() owns a Scope, released on close: ', live, '. ',
    'Disposing the page scope drops the right-click binding and every readout effect.');
  main.append(note, button('Dispose page scope', () => {
    context.close();
    if (anchored) anchored.close();
    if (watch) watch.dispose();
    page.dispose();
    live.textContent = `${Scope.liveCount} live scopes after disposal`;
  }));

  injectOnce('u2-menu-bar-css', () => {
    const l = document.createElement('link');
    l.rel = 'stylesheet';
    l.href = new URL('../../css/menu-bar.css', import.meta.url).href;
    return l;
  });
  const {MenuBar} = await import('../../src/components/menu-bar.js');

  main.append(el('h1', null, 'Menu bar'));
  const barIntro = el('p');
  barIntro.innerHTML = 'The classic top menu: one compact row of text items, each opening a fresh ' +
    '<code>Menu</code> anchored under it. Click opens, clicking the open item closes, clicking ' +
    'another switches. While a menu is open, <b>hovering another item switches to it instantly</b> ' +
    '— no intent delay, unlike submenus. One tab stop: ←/→ move between items (and switch the open ' +
    'menu), Home/End jump to the ends, ↓ / Enter / Space open, Esc closes and returns focus to the ' +
    'item. The bar paints no background of its own, so it drops into a ribbon unchanged.';
  main.append(barIntro);

  const barScope = new Scope();
  const barWrap = signal(true);
  const layout = signal('Grid');
  const view = (name) => ({check: layout.peek() === name});

  const bar = new MenuBar()
    .menu('File', (m) => m
      .item('New', () => last.value = 'File → New', {shortcut: 'Ctrl+N'})
      .item('Open…', () => last.value = 'File → Open', {shortcut: 'Ctrl+O'})
      .separator()
      .group('Recent')
        .item('demog.csv', () => last.value = 'Recent → demog.csv')
        .item('cars.csv', () => last.value = 'Recent → cars.csv')
        .item('sales.csv', () => last.value = 'Recent → sales.csv')
      .endGroup())
    .menu('Edit', (m) => m
      .item('Undo', () => last.value = 'Edit → Undo', {shortcut: 'Ctrl+Z'})
      .item('Redo', () => last.value = 'Edit → Redo', {shortcut: 'Ctrl+Y'})
      .separator()
      .item('Word wrap', () => barWrap.value = !barWrap.peek(), {check: barWrap.peek(), shortcut: 'Alt+Z'}))
    .menu('View', (m) => m
      .item('Grid', () => layout.value = 'Grid', view('Grid'))
      .item('Cards', () => layout.value = 'Cards', view('Cards'))
      .item('Chart', () => layout.value = 'Chart', view('Chart')))
    .menu('Tools', (m) => m.item('Options…', () => last.value = 'Tools → Options'))
    .menu('Help', (m) => m
      .item('Documentation', () => last.value = 'Help → Documentation', {shortcut: 'F1'})
      .separator()
      .item('About', () => last.value = 'Help → About'));
  bar.setEnabled('Tools', false);

  const ribbon = el('div');
  ribbon.style.cssText = 'width: 720px; max-width: 100%; background: var(--dg-bg-secondary);' +
    'border: var(--dg-border-width) solid var(--dg-border); border-radius: var(--dg-radius);';
  ribbon.append(bar.root);
  main.append(ribbon);

  const openText = el('code');
  const barWrapText = el('code');
  const layoutText = el('code');
  bindText(barScope, openText, computed(() => bar.openMenu.value ?? '(closed)'));
  bindText(barScope, barWrapText, computed(() => String(barWrap.value)));
  bindText(barScope, layoutText, computed(() => layout.value));
  const barStatus = el('p', 'u2-gallery-status');
  barStatus.append('openMenu = ', openText, '  ·  word wrap = ', barWrapText, '  ·  view = ', layoutText);
  main.append(barStatus);

  let toolsEnabled = false;
  const toolsToggle = button('Enable Tools', () => {
    toolsEnabled = !toolsEnabled;
    bar.setEnabled('Tools', toolsEnabled);
    toolsToggle.textContent = toolsEnabled ? 'Disable Tools' : 'Enable Tools';
  });
  main.append(toolsToggle);
  main.append(el('p', 'u2-gallery-status',
    'Tools starts disabled: it never opens, ←/→ skip over it, and hovering it while another menu ' +
    'is open changes nothing. setEnabled() puts it back, and disabling the open item closes it.'));

  main.append(el('h2', null, 'Disposal'));
  const barLive = el('code', null, `${Scope.liveCount} live scopes`);
  const barNote = el('p', 'u2-gallery-status');
  barNote.append('Dispose closes the open menu and releases every item listener and the readouts: ',
    barLive, '.');
  main.append(barNote, button('Dispose menu bar', () => {
    bar.dispose();
    barScope.dispose();
    barLive.textContent = `${Scope.liveCount} live scopes after disposal`;
  }));
}
