import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, divV, span, button} from '../../src/core/elements.js';
import {iconButton, buttonWithIcon, dropDownButton} from '../../src/components/actions/buttons.js';
import {ButtonGroup} from '../../src/components/actions/button-group.js';
import {Toolbar} from '../../src/components/navigation/toolbar.js';
import {Form} from '../../src/components/forms/form.js';
import {TextInput} from '../../src/components/inputs/text-input.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

// vendored Font Awesome free maps onto the platform's icon classes for standalone hosts
document.documentElement.classList.add('u2-standalone');
injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-icons-css', '../../css/icons.css');
injectOnce('u2-buttons-css', '../../css/buttons.css');
injectOnce('u2-toolbar-css', '../../css/toolbar.css');
injectOnce('u2-menu-css', '../../css/menu.css');
injectOnce('u2-tooltip-css', '../../css/tooltip.css');
injectOnce('u2-inputs-css', '../../css/inputs.css');
injectOnce('u2-form-css', '../../css/form.css');

const ALIGN = [
  {id: 'left', label: 'Left', icon: 'align-left'},
  {id: 'center', label: 'Center', icon: 'align-center'},
  {id: 'right', label: 'Right', icon: 'align-right'},
];

const VIEWS = [
  {id: 'grid', label: 'Grid', icon: 'table'},
  {id: 'cards', label: 'Cards', icon: 'th-large'},
  {id: 'chart', label: 'Chart', icon: 'chart-bar'},
];

const ACTIONS = [
  {id: 'add', label: 'Add row', icon: 'plus'},
  {id: 'copy', label: 'Copy', icon: 'copy'},
  {id: 'delete', label: 'Delete', icon: 'trash-alt'},
  {id: 'refresh', label: 'Refresh', icon: 'sync'},
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
  r.style.gap = 'var(--dg-space-m)';
  r.style.flexWrap = 'wrap';
  return r;
}

function labelled(text, node) {
  const box = divV([el('div', 'u2-gallery-status', text), node]);
  box.style.gap = 'var(--dg-space-s)';
  box.style.alignItems = 'flex-start';
  return box;
}

export async function render(main) {
  main.append(el('h1', null, 'Buttons & groups'));
  const intro = el('p');
  intro.innerHTML = 'Icon-capable buttons compose the core <code>button()</code> — the icon layer ' +
    'lives in <code>components/actions/buttons.ts</code> because core may not import components. ' +
    '<code>ButtonGroup</code> is one control for three jobs: an action group, a radio group ' +
    '(<code>single</code> — clicking the active item never deselects it) and independent toggles ' +
    '(<code>multi</code>). Icon groups are the same control in <code>iconOnly</code> mode, and the ' +
    'density axis (<code>form</code> / <code>toolbar</code> / <code>ribbon</code>) is what makes it ' +
    'fit a form row, a toolbar or a ribbon.';
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

  const last = signal('(none)');
  const act = (name) => () => last.value = name;

  section('Buttons with icons', () => [
    row([
      buttonWithIcon('Save', 'save', act('Save')),
      buttonWithIcon('Open', 'folder-open', act('Open')),
      buttonWithIcon('Next', 'arrow-right', act('Next'), {iconRight: true}),
      buttonWithIcon('Add', 'plus', act('Add'), {primary: true}),
      buttonWithIcon('Continue', 'chevron-right', act('Continue'), {primary: true, iconRight: true}),
      buttonWithIcon('Delete', 'trash-alt', act('Delete'), {tooltip: 'Delete the selected rows'}),
      button('No icon', act('No icon')),
    ]),
    readout('last action', last),
  ]);

  const starred = signal(false);
  section('Icon buttons', () => [
    row([
      iconButton('pencil-alt', act('Edit'), {tooltip: 'Edit'}),
      iconButton('trash-alt', act('Delete'), {tooltip: 'Delete'}),
      iconButton('sync', act('Refresh'), {tooltip: 'Refresh'}),
      iconButton('ellipsis-v', act('More'), {tooltip: 'More actions'}),
      iconButton('star', act('Star'), {tooltip: 'Star this item', variant: 'solid', toggle: starred}),
    ]),
    readout('starred', computed(() => String(starred.value))),
  ]);

  section('Segmented — single (radio semantics)', () => {
    const align = new ButtonGroup({items: ALIGN, toggle: 'single'});
    return [
      align,
      readout('selected', computed(() => align.selected.value.join(', ') || '(none)')),
      el('div', 'u2-gallery-status',
        'role=radiogroup with role=radio cells and aria-checked. Tab enters once, ←/→ move, ' +
        'Enter or Space activates; clicking the active cell keeps it selected.'),
    ];
  });

  section('Segmented — multi (icon + label)', () => {
    const views = new ButtonGroup({items: VIEWS, toggle: 'multi'});
    return [
      views,
      readout('selected', computed(() => views.selected.value.join(', ') || '(none)')),
      el('div', 'u2-gallery-status', 'Plain buttons with aria-pressed; every write is a fresh array.'),
    ];
  });

  section('Icon-only action group', () => [
    new ButtonGroup({items: ACTIONS.map((i) => ({...i, onClick: act(i.label)})), iconOnly: true}),
    readout('last action', last),
    el('div', 'u2-gallery-status', 'toggle: none — nothing is ever selected, and each label ' +
      'becomes the tooltip and the accessible name.'),
  ]);

  section('The same group at three densities', () => {
    const at = (density) => labelled(density, new ButtonGroup({
      items: ACTIONS.map((i) => ({...i, onClick: act(`${i.label} (${density})`)})),
      iconOnly: true,
      density,
    }).root);
    const strip = divH([at('form'), at('toolbar'), at('ribbon')]);
    strip.style.gap = 'var(--dg-space-xxl)';
    strip.style.alignItems = 'flex-start';
    return [strip, el('div', 'u2-gallery-status',
      'form keeps the segmented frame at the button metric; toolbar and ribbon go flat at the ' +
      'inline metric — the frame appears on hover, and the ribbon packs tighter still.')];
  });

  section('Inside a toolbar', () => {
    const bar = new Toolbar({ariaLabel: 'Document actions'});
    const zoom = new ButtonGroup({
      items: [
        {id: 'out', icon: 'search-minus', label: 'Zoom out', onClick: act('Zoom out')},
        {id: 'in', icon: 'search-plus', label: 'Zoom in', onClick: act('Zoom in')},
      ],
      iconOnly: true,
      density: 'toolbar',
    });
    const mode = new ButtonGroup({items: ALIGN, toggle: 'single', iconOnly: true, density: 'toolbar'});
    bar
      .addButton('New', act('New'))
      .addButton('Open', act('Open'))
      .addSeparator()
      .add(zoom.root)
      .addSeparator()
      .add(mode.root)
      .addSeparator()
      .add(dropDownButton('Export', (m) => m
        .item('CSV', act('Export CSV'))
        .item('Excel', act('Export Excel'))));
    const host = el('div');
    host.style.cssText = 'width: 520px; max-width: 100%; border: var(--dg-border-width) solid ' +
      'var(--dg-border); border-radius: var(--dg-radius); overflow: hidden;';
    host.append(bar.root);
    return [host, readout('alignment', computed(() => mode.selected.value.join(', ') || '(none)'))];
  });

  const saved = signal('nothing saved yet');
  section('Inside a form', () => {
    const form = new Form();
    const name = new TextInput({label: 'Compound name', value: 'Aspirin'});
    const notes = new TextInput({label: 'Notes', placeholder: 'Optional'});
    const layout = new ButtonGroup({items: ALIGN, toggle: 'single'});
    form.addAll([name, notes]);
    form.addButtons((buttons) => {
      buttons.append(
        layout.root,
        buttonWithIcon('Save', 'save', () => saved.value =
          `${name.value.value} · ${layout.selected.value.join(', ') || 'no alignment'}`, {primary: true})
      );
    });
    return [form, readout('saved', saved)];
  });

  section('Disabled items', () => {
    const guarded = new ButtonGroup({
      items: [
        {id: 'run', label: 'Run', icon: 'play', onClick: act('Run')},
        {id: 'stop', label: 'Stop', icon: 'stop', onClick: act('Stop'), enabled: false},
        {id: 'reset', label: 'Reset', icon: 'undo', onClick: act('Reset')},
      ],
    });
    const stopEnabled = signal(false);
    return [
      guarded,
      row([button(computed(() => stopEnabled.value ? 'Disable Stop' : 'Enable Stop'), () => {
        stopEnabled.value = !stopEnabled.peek();
        guarded.setEnabled('stop', stopEnabled.peek());
      })]),
      el('div', 'u2-gallery-status',
        'A disabled cell drops out of the roving tab order, so ←/→ skip straight over it.'),
    ];
  });

  section('Drop-down buttons', () => [
    row([
      dropDownButton('Export', (m) => m
        .item('CSV', act('Export CSV'), {shortcut: 'Ctrl+E'})
        .item('Excel', act('Export Excel'))
        .item('JSON', act('Export JSON'))
        .separator()
        .group('Advanced')
        .item('Template…', act('Export template'))
        .item('Selected rows only', act('Export selection'))
        .endGroup()),
      dropDownButton('Add', (m) => m
        .item('Row', act('Add row'))
        .item('Column', act('Add column')), {icon: 'plus', tooltip: 'Add a row or a column'}),
      dropDownButton('Save', (m) => m
        .item('Save', act('Save'))
        .item('Save as…', act('Save as')), {icon: 'save', primary: true}),
    ]),
    readout('last action', last),
    el('div', 'u2-gallery-status', 'The whole button opens the menu, and the callback runs on ' +
      'every open — a fresh menu always shows current state. ArrowDown opens it from the keyboard.'),
  ]);

  section('Split buttons', () => [
    row([
      dropDownButton('Run', (m) => m
        .item('Run with parameters…', act('Run with parameters'))
        .item('Run in background', act('Run in background'))
        .separator()
        .item('Edit configuration…', act('Edit configuration')), {icon: 'play', split: act('Run')}),
      dropDownButton('Publish', (m) => m
        .item('Publish and share', act('Publish and share'))
        .item('Publish as draft', act('Publish as draft')), {primary: true, split: act('Publish')}),
    ]),
    readout('last action', last),
    el('div', 'u2-gallery-status', 'Two attached segments: the main part fires its own action, ' +
      'only the arrow opens the menu, and primary skins the whole control.'),
  ]);

  section('Open state', () => {
    const expanded = signal('false');
    const drop = dropDownButton('Watch me', (m) => m
      .item('One', act('One'))
      .item('Two', act('Two')));
    const observer = new MutationObserver(() => expanded.value = drop.getAttribute('aria-expanded'));
    observer.observe(drop, {attributes: true, attributeFilter: ['aria-expanded']});
    Scope.ambient.own(() => observer.disconnect());
    return [
      row([drop]),
      readout('aria-expanded', expanded),
      el('div', 'u2-gallery-status', 'The attribute is written by an effect on the menu instance ' +
        'that lives exactly as long as the open menu.'),
    ];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Every section is a Control.build(...) builder that owns its groups, toggle effects, click ' +
    'listeners and tooltip bindings: disposing the sections releases all of them and live scopes ' +
    'drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
