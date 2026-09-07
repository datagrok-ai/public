import {
  signal, computed, Signal, Scope,
  divV, divH, panel, span, h3, button,
  TextInput, Splitter, Accordion, Breadcrumbs, Toolbar,
  Menu, MenuBar, Dialog, Tooltip, Form,
} from '@datagrok-libraries/u2';
import {readout} from './common';

export function containersPage(): HTMLElement {
  const left = panel([h3('Navigator'), span('30% wide; drag the sash, or focus it and use arrow keys.')]);
  const right = panel([h3('Content'), span('Sizes live in a signal — watch them change:')]);
  const split = new Splitter([left, right], {direction: 'horizontal', sizes: [30, 70]});
  split.root.classList.add('u2demo-split');
  right.append(readout('sizes', computed(() =>
    split.sizes.value.map((s) => s.toFixed(1)).join(' / '))));

  const acc = new Accordion();
  acc.addPane('General', divV([span('Eagerly built pane.')]), true);
  acc.addPane('Advanced', () => divV([span(`Lazily built on first expand at ${new Date().toLocaleTimeString()}.`)]));

  const path = signal(['Home', 'Projects', 'Demo', 'Tables', 'demog']);
  const crumbs = new Breadcrumbs({
    items: path.value,
    onClick: (index) => {
      path.value = path.value.slice(0, index + 1);
      crumbs.setItems(path.value);
    },
  });

  const compact = signal(false);
  const action = signal('');
  const bar = new Toolbar({ariaLabel: 'Panel toolbar'});
  bar.addButton('Run', () => action.value = 'Run clicked', {tooltip: 'Pretend to run'});
  bar.addButton('Save', () => action.value = 'Saved', {tooltip: 'Pretend to save'});
  bar.addSeparator();
  bar.addToggle('Compact', compact, {tooltip: 'Toolbar toggles bind to a signal'});
  bar.addMenu('Help', (m) => m
    .item('About u2', () => action.value = 'u2: next-gen Datagrok UI library')
    .item('Gallery', () => window.open('https://github.com/datagrok-ai/public/tree/master/libraries/u2')));

  return divV([
    h3('Splitter'), split,
    h3('Accordion'), acc,
    h3('Breadcrumbs'), crumbs,
    span('Click a segment to truncate the path.', 'u2demo-hint'),
    readout('path', computed(() => crumbs.items.value.join(' > '))),
    h3('Toolbar (panel-local)'),
    panel([bar, span('Inline toolbars belong to a panel inside the view; the view\'s own commands ' +
      'live in the ribbon. Narrow the window — items collapse into the » menu.', 'u2demo-hint')]),
    readout('toolbar', computed(() => `${action.value || '(nothing yet)'} · compact ${compact.value}`)),
  ], 'u2demo-page');
}

export function popupsPage(status: Signal<string>): HTMLElement {
  const scope = Scope.ambient!;
  const log = (action: string) => status.value = action;

  const autoRefresh = signal(false);
  const bar = new MenuBar()
    .menu('File', (m) => m
      .item('New session', () => log('New session'), {shortcut: 'Ctrl+N'})
      .item('Save layout', () => log('Layout saved'))
      .separator()
      .group('Export')
      .item('As JSON', () => log('Exported as JSON'))
      .item('As CSV', () => log('Exported as CSV'))
      .endGroup())
    .menu('View', (m) => m
      .item('Auto-refresh', () => {
        autoRefresh.value = !autoRefresh.value;
        log(`Auto-refresh ${autoRefresh.value ? 'on' : 'off'}`);
      }, {check: autoRefresh.value})
      .item('Reset panels', () => log('Panels reset')))
    .menu('Help', (m) => m
      .item('About u2', () => log('u2: next-gen Datagrok UI library'))
      .item('Sources', () => window.open('https://github.com/datagrok-ai/public/tree/master/libraries/u2')));

  const menuButton = button('Open menu', () => Menu.popup()
    .item('Add to favorites', () => log('Added to favorites'), {shortcut: 'Ctrl+D'})
    .item('Rename', () => log('Renamed'))
    .separator()
    .group('Sort by')
    .item('Name', () => log('Sorted by name'))
    .item('Size', () => log('Sorted by size'))
    .endGroup()
    .item('Delete', () => log('Deleted'), {enabled: false})
    .show({anchor: menuButton}));

  const zone = panel([span('Right-click anywhere in this panel for a context menu.')]);
  zone.classList.add('u2demo-card');
  const context = new Menu()
    .item('Copy', () => log('Copied'))
    .item('Paste', () => log('Pasted'))
    .separator()
    .item('Select all', () => log('Selected all'));
  context.bindContext(zone, scope);

  const hoverMe = span('Hover me for a tooltip.', 'u2demo-card');
  Tooltip.bind(hoverMe, () => `One shared tooltip element, evaluated at show time: ${new Date().toLocaleTimeString()}`,
    scope);

  const dlgName = new TextInput({label: 'Name', value: 'u2'});
  dlgName.addValidator((v) => v.trim() ? null : 'Required');
  const dlgForm = new Form().add(dlgName);
  const dialog = Dialog.create('u2 dialog')
    .add(span('Modal, draggable by the title bar, Esc cancels, Enter accepts.'))
    .add(dlgForm)
    .onOK(() => log(`Dialog OK, name = ${dlgName.value.value}`))
    .onCancel(() => log('Dialog cancelled'));

  return divV([
    h3('Menu bar (panel-local)'),
    panel([bar, span('A menu bar belongs to the panel it commands; a view\'s own commands live ' +
      'in the ribbon. Every item here writes the shell status bar.', 'u2demo-hint')]),
    h3('Menus, dialog and tooltip'),
    divH([menuButton, button('Open dialog', () => dialog.show({modal: true}))], 'u2demo-row'),
    zone, hoverMe,
  ], 'u2demo-page');
}
