import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {
  signal, computed, Signal, ReadonlySignal, Component, Scope,
  divV, divH, panel, span, h3, button, bigButton, link,
  TextInput, TextArea, BoolInput, NumberInput, ChoiceInput, MultiChoiceInput,
  Combobox, TypeAhead, AsyncView, VirtualList, VirtualTree, TreeNode,
  TabStrip, Splitter, Accordion, Breadcrumbs, Toolbar,
  Menu, Dialog, Tooltip, Form, PropertyGrid,
} from '@datagrok-libraries/u2';
import * as DG from 'datagrok-api/dg';
import {asDartInput, tableInput, columnInput, leakReport, propertyForm, PropertyLike,
  userInput, entityInput, chip, EntityChip, entityCard, handlerRenderer,
  moleculeInput, moleculeRenderer}
  from '@datagrok-libraries/u2/src/dg/index.js';

const CITIES = [
  'Amsterdam', 'Athens', 'Barcelona', 'Berlin', 'Boston', 'Brussels', 'Budapest', 'Chicago',
  'Copenhagen', 'Dublin', 'Geneva', 'Helsinki', 'Kyiv', 'Lisbon', 'London', 'Madrid', 'Milan',
  'Munich', 'New York', 'Oslo', 'Paris', 'Prague', 'Rome', 'Seattle', 'Stockholm', 'Tokyo',
  'Toronto', 'Vienna', 'Warsaw', 'Zurich',
];

function readout(name: string, source: string | ReadonlySignal<unknown>): HTMLElement {
  return divH([span(`${name} = `), span(source)], 'u2demo-status');
}

function inputsPage(): HTMLElement {
  const query = signal('');
  const notify = signal(true);

  const name = new TextInput({label: 'Name', value: 'Aspirin', tooltipText: 'Compound name'});
  const search = new TextInput({label: 'Search', search: true, bind: query, placeholder: 'Filter compounds…'});
  const notes = new TextArea({label: 'Notes', placeholder: 'Grows as you type…', autoGrow: true});
  const enabled = new BoolInput({label: 'Enabled', value: true});
  const notifications = new BoolInput({label: 'Notifications', switch: true, bind: notify});
  const dose = new NumberInput({label: 'Dose, mg', mode: 'float', min: 0, max: 1000, step: 0.5, value: 250});
  const replicates = new NumberInput({label: 'Replicates', mode: 'int', min: 1, max: 10, value: 3});
  const stage = new ChoiceInput({label: 'Stage', items: ['Discovery', 'Preclinical', 'Phase I', 'Phase II'],
    value: 'Discovery'});
  const targets = new MultiChoiceInput({label: 'Targets', items: ['GPCR', 'Kinase', 'Ion channel', 'Protease'],
    value: ['Kinase']});

  const code = new TextInput({label: 'Code', placeholder: 'required, up to 10 characters'});
  code.addValidator((v) => v.trim().length > 0 ? null : 'Value is required');
  code.addValidator((v) => v.length > 10 ? 'At most 10 characters' : null);

  return divV([
    span('Every input is a signal owner: bind to it, or pass your own signal via `bind`.', 'u2demo-hint'),
    name, search, notes, enabled, notifications, dose, replicates, stage, targets, code,
    readout('search', computed(() => query.value || '(empty)')),
    readout('notifications', computed(() => String(notify.value))),
    readout('dose * replicates', computed(() =>
      String((dose.value.value ?? 0) * (replicates.value.value ?? 0)))),
    readout('code validity', computed(() => code.validity.value ?? 'valid')),
  ], 'u2demo-page');
}

function asyncPage(): HTMLElement {
  const local = new Combobox({items: CITIES, placeholder: 'Pick a city…'});

  const fetchCities = async (q: string, abort: AbortSignal): Promise<string[]> => {
    await new Promise((resolve) => setTimeout(resolve, 400));
    if (abort.aborted)
      return [];
    return CITIES.filter((c) => c.toLowerCase().includes(q.toLowerCase()));
  };

  const remote = new Combobox({items: fetchCities, placeholder: 'Type to search (simulated 400 ms latency)…'});

  const view = AsyncView.owned(fetchCities,
    (items) => divV(items.map((c) => span(c)), 'u2demo-card'),
    {empty: 'No cities match', skeleton: true});
  view.refresh();
  const filter = new TextInput({inline: true, search: true, placeholder: 'Filter the async view…',
    onChanged: (v) => view.refresh(v)});

  return divV([
    span('One async contract — debounce, cancellation, loading/empty/error — shared by every async control.',
      'u2demo-hint'),
    h3('Combobox, local items'), local,
    h3('Combobox, async items'), remote,
    readout('local', local.value), readout('remote', remote.value),
    h3('AsyncView'), filter, view,
  ], 'u2demo-page');
}

function listsPage(): HTMLElement {
  const list = new VirtualList<number>({
    keyOf: (i) => String(i),
    render: (item, index, row) => {
      row.title = `row ${index}`;
      return span(`item ${item} ${'▪'.repeat(1 + item % 7)}`);
    },
  });
  list.setItems(Array.from({length: 100000}, (_, i) => i));
  list.root.classList.add('u2demo-list');

  const projects: TreeNode<string>[] = [{
    id: 'demo', label: 'Demo project', children: [
      {id: 'demo/tables', label: 'Tables', children: [
        {id: 'demo/tables/demog', label: 'demog'},
        {id: 'demo/tables/cars', label: 'cars'},
      ]},
      {id: 'demo/queries', label: 'Queries', children: [
        {id: 'demo/queries/orders', label: 'orders by country'},
      ]},
    ],
  }, {
    id: 'remote', label: 'Server (lazy)', children: async () => {
      await new Promise((resolve) => setTimeout(resolve, 600));
      return Array.from({length: 5}, (_, i) => ({id: `remote/${i}`, label: `dataset-${i}.csv`}));
    },
  }];
  const tree = new VirtualTree<string>();
  tree.setRoots(projects);
  tree.root.classList.add('u2demo-list');

  return divV([
    span('Both containers render only the visible rows — the list below holds 100,000 items.', 'u2demo-hint'),
    h3('VirtualList'), list,
    readout('selectedIndex', computed(() => String(list.selectedIndex.value))),
    h3('VirtualTree'), tree,
    divH([
      button('Reveal demog', () => void tree.expandPath(['demo', 'demo/tables', 'demo/tables/demog'])),
      span('— expands ancestors (awaiting lazy loads), then selects and reveals.'),
    ], 'u2demo-row'),
    readout('selectedNode', computed(() => tree.selectedNode.value?.label ?? '(none)')),
  ], 'u2demo-page');
}

function layoutPage(): HTMLElement {
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

  return divV([
    h3('Splitter'), split,
    h3('Accordion'), acc,
    h3('Breadcrumbs'), crumbs,
    span('Click a segment to truncate the path.', 'u2demo-hint'),
    readout('path', computed(() => crumbs.items.value.join(' > '))),
  ], 'u2demo-page');
}

function popupsPage(status: Signal<string>): HTMLElement {
  const scope = Scope.ambient!;
  const log = (action: string) => status.value = action;

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
    divH([menuButton, button('Open dialog', () => dialog.show({modal: true}))], 'u2demo-row'),
    zone, hoverMe,
    readout('last action', computed(() => status.value)),
  ], 'u2demo-page');
}

function formPage(): HTMLElement {
  const first = new TextInput({label: 'First name', name: 'first', value: 'Ada'});
  first.addValidator((v) => v.trim() ? null : 'Required');
  const last = new TextInput({label: 'Last name', name: 'last'});
  last.addValidator((v) => v.trim() ? null : 'Required');
  const age = new NumberInput({label: 'Age', name: 'age', mode: 'int', min: 0, max: 120, value: 36});
  const role = new ChoiceInput({label: 'Role', name: 'role', items: ['Viewer', 'Editor', 'Admin'], value: 'Editor'});
  const subscribed = new BoolInput({label: 'Subscribe', name: 'subscribed', switch: true, value: true});

  const form = new Form();
  form.addAll([first, last, age, role, subscribed]);
  form.addButtons((row) => {
    row.append(bigButton('Submit', () => {
      if (form.validate())
        grok.shell.info(JSON.stringify(form.getValues()));
    }));
    row.append(button('Reset', () => form.setValues({first: 'Ada', last: '', age: 36, role: 'Editor',
      subscribed: true})));
  });

  const pg = new PropertyGrid();
  pg.setProperties([
    {name: 'title', type: 'string', description: 'Viewer title'},
    {name: 'opacity', type: 'float', min: 0, max: 1},
    {name: 'width', type: 'int', min: 100, max: 2000, category: 'Layout'},
    {name: 'height', type: 'int', min: 100, max: 2000, category: 'Layout'},
    {name: 'showLegend', type: 'bool', category: 'Legend'},
    {name: 'position', type: 'choice', choices: ['left', 'right', 'top', 'bottom'], category: 'Legend'},
  ], {title: 'Scatter plot', opacity: 0.8, width: 600, height: 400, showLegend: true, position: 'right'});

  return divV([
    h3('Form'),
    span('Aligned label column, aggregated validity, Enter submits.', 'u2demo-hint'),
    form,
    readout('form validity', computed(() => form.validity.value ?? 'valid')),
    h3('PropertyGrid'),
    span('The value record is replaced on every edit — bind to it directly.', 'u2demo-hint'),
    pg,
    readout('last change', computed(() => {
      const c = pg.onChanged.value;
      return c ? `${c.name} → ${String(c.value)}` : '(none)';
    })),
    readout('values', computed(() => JSON.stringify(pg.values.value))),
  ], 'u2demo-page');
}

function platformPage(): HTMLElement {
  const demog = grok.data.demo.demog(100);
  const table = tableInput('Open table');
  const column = columnInput('Demog column', demog);

  const bridged = new TextInput({label: 'Compound', value: 'Aspirin'});
  bridged.addValidator((v) => v.trim() ? null : 'Required');
  const dart = asDartInput(bridged);

  const leaks = signal('(press the button)');

  return divV([
    span('The dg entry point is the only layer that touches datagrok-api: ' +
      'host() joins the DG.Widget kill-walk, asDartInput() makes any u2 input a DG.InputBase.', 'u2demo-hint'),
    h3('Pickers'),
    table, column,
    readout('table', table.value), readout('column', column.value),
    divH([
      button('Add demog to workspace', () => grok.shell.addTableView(demog)),
      span('— then reopen this tab to see the pickers follow.'),
    ], 'u2demo-row'),
    h3('DartInput bridge'),
    divH([
      button('Open a DG.Dialog with a u2 input', () => ui.dialog('u2 input, platform dialog')
        .add(dart)
        .onOK(() => grok.shell.info(`Compound: ${bridged.value.value}`))
        .show()),
    ], 'u2demo-row'),
    readout('bridged value', bridged.value),
    h3('Leak detector'),
    divH([
      button('leakReport()', () => leaks.value = JSON.stringify(leakReport())),
      span(leaks),
    ], 'u2demo-row'),
    span('Close this view and call it again from another U2 Demo instance: both numbers return to baseline.',
      'u2demo-hint'),
  ], 'u2demo-page');
}

function objectFormPage(): HTMLElement {
  const group = DG.Group.create(`u2demo-${Math.random().toString(36).slice(2, 8)}`);
  const props: PropertyLike[] = [
    {name: 'friendlyName', caption: 'Name', type: 'string', nullable: false,
      get: (g) => (g as DG.Group).friendlyName, set: (g, v) => (g as DG.Group).friendlyName = v as string},
    {name: 'personal', caption: 'Personal', type: 'bool',
      get: (g) => (g as DG.Group).personal, set: (g, v) => (g as DG.Group).personal = v as boolean},
    {name: 'hidden', caption: 'Hidden', type: 'bool',
      get: (g) => (g as DG.Group).hidden, set: (g, v) => (g as DG.Group).hidden = v as boolean},
  ];
  const form = propertyForm(props, group);
  const result = signal('(not saved)');
  form.addButtons((row) => {
    row.append(
      bigButton('SAVE', async () => {
        if (!form.validate()) {
          result.value = 'Fix validation first';
          return;
        }
        const saved = await grok.dapi.groups.save(group);
        result.value = `Saved: ${saved.friendlyName} (${saved.id})`;
      }),
      button('Delete', async () => {
        await grok.dapi.groups.delete(group);
        result.value = 'Deleted';
      }));
  });
  return divV([
    span('propertyForm(): everything below is generated from PropertyLike metadata over a real DG.Group — ' +
      'editors chosen by type, required from nullable: false, edits write through to the entity, ' +
      'SAVE persists it via grok.dapi.groups.', 'u2demo-hint'),
    form.root,
    readout('result', result),
  ], 'u2demo-page');
}

function entitiesPage(): HTMLElement {
  const scope = Scope.ambient!;
  const user = userInput('Search users…');
  const group = entityInput('Search groups…', () => grok.dapi.groups,
    {filter: (q) => q ? `friendlyName like "${q}"` : undefined});

  const picked = divH([], 'u2demo-row');
  const chips: EntityChip[] = [];
  scope.effect(() => {
    const items = [user.selected.value, group.selected.value].filter((x) => x != null);
    for (const c of chips)
      c.dispose();
    chips.length = 0;
    picked.textContent = '';
    for (const x of items) {
      const c = chip(x);
      chips.push(c);
      picked.append(c.root);
    }
  });

  return divV([
    span('Entity pickers: TypeAhead + dapiSource(() => grok.dapi.*) + the object\'s own ObjectHandler ' +
      'rendering. The picked value is the entity itself; chips render handler markup, show the handler ' +
      'tooltip, and click makes the object current (see the context panel).', 'u2demo-hint'),
    h3('User picker (curated rows)'), user,
    h3('Group picker (generic handlerRenderer)'), group,
    h3('Picked, as chips'), picked,
    h3('Current user'), divH([chip(grok.shell.user).root], 'u2demo-row'),
    readout('user', computed(() => user.selected.value?.friendlyName ?? '(none)')),
    readout('group', computed(() => group.selected.value?.friendlyName ?? '(none)')),
  ], 'u2demo-page');
}

function spacesPage(): HTMLElement {
  const scope = Scope.ambient!;
  // the spaces endpoint does not honor smart filters (probed 2026-08-14); a deployment has a
  // handful of spaces, so fetch once and match client-side
  const space = new TypeAhead<DG.Project>({
    source: async (q) => (await grok.dapi.spaces.list({pageSize: 100}))
      .filter((s) => s.friendlyName.toLowerCase().includes(q.toLowerCase())),
    renderer: handlerRenderer<DG.Project>(),
    placeholder: 'Choose a space…',
  });
  const project = entityInput<DG.Project>('Choose a project or dashboard…', () => grok.dapi.projects,
    {filter: (q) => q ? `friendlyName like "${q}"` : undefined});

  const picked = signal<DG.Project | null>(null);
  scope.effect(() => { const s = space.selected.value; if (s) picked.value = s; });
  scope.effect(() => { const p = project.selected.value; if (p) picked.value = p; });

  const preview = divV([], 'u2demo-row');
  let shown: Scope | undefined;
  scope.own(() => shown?.dispose());
  scope.effect(() => {
    const x = picked.value;
    shown?.dispose();
    shown = new Scope();
    preview.replaceChildren(Scope.runWith(shown, () => x ?
      divV([
        entityCard(x).root,
        divH([chip(x).root, button('Open', () => void x.open())], 'u2demo-row'),
      ], 'u2demo-row') :
      span('Pick a space or a project to preview it.', 'u2demo-hint')));
  });

  return divV([
    span('Both pickers are entityInput over grok.dapi (spaces ARE projects server-side, so one ' +
      'ObjectHandler serves both). The preview is the handler\'s own renderCard — the platform\'s ' +
      'gallery card, framed by u2 — and clicking it or the chip makes the entity current.', 'u2demo-hint'),
    h3('Space'), space,
    h3('Project / dashboard'), project,
    h3('Preview'), preview,
  ], 'u2demo-page');
}

const DRUGS = [
  {name: 'Aspirin', smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O'},
  {name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'},
  {name: 'Ibuprofen', smiles: 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O'},
  {name: 'Paracetamol', smiles: 'CC(=O)NC1=CC=C(C=C1)O'},
  {name: 'Naproxen', smiles: 'CC(C1=CC2=CC=C(C=C2C=C1)OC)C(=O)O'},
];

function moleculesPage(): HTMLElement {
  const structure = moleculeInput('Structure');
  structure.value.value = DRUGS[0].smiles;

  const compound = {name: 'Aspirin', smiles: DRUGS[0].smiles, mw: 180.16};
  const props: PropertyLike[] = [
    {name: 'name', caption: 'Name', type: 'string', nullable: false,
      get: (o) => (o as typeof compound).name, set: (o, v) => (o as typeof compound).name = v as string},
    {name: 'smiles', caption: 'Structure', type: 'string', semType: 'Molecule',
      get: (o) => (o as typeof compound).smiles, set: (o, v) => (o as typeof compound).smiles = v as string},
    {name: 'mw', caption: 'MW, Da', type: 'float', min: 0,
      get: (o) => (o as typeof compound).mw, set: (o, v) => (o as typeof compound).mw = v as number},
  ];
  const edits = signal(0);
  const form = propertyForm(props, compound, {onChanged: () => edits.value++});

  const mols = moleculeRenderer({width: 110, height: 70});
  const chips = divH(DRUGS.map((d) =>
    chip(d.smiles, {renderer: mols, currentOnClick: false}).root), 'u2demo-row');

  const rowMol = moleculeRenderer({width: 90, height: 55});
  const picker = new TypeAhead<{name: string, smiles: string}>({
    source: DRUGS,
    renderer: {
      caption: (d) => d.name,
      listItem: (d) => divH([rowMol.listItem!(d.smiles), span(d.name)], 'u2demo-mol-row'),
    },
    placeholder: 'Find a compound…',
  });

  return divV([
    span('Everything below is bridged, not ported: the sketcher input comes through ' +
      'fromDartInput(ui.input.molecule), the form picks it by semType through the u2 editor ' +
      'registry, and depictions render through Chem\'s drawMolecule — u2 contains no chemistry.',
    'u2demo-hint'),
    h3('moleculeInput — the platform sketcher as a u2 Input<string>'),
    structure,
    readout('smiles', computed(() => structure.value.value || '(empty)')),
    h3('propertyForm — the Structure field is picked by semType: Molecule'),
    form.root,
    readout('compound', computed(() => {
      edits.value;
      return JSON.stringify(compound);
    })),
    h3('Structure chips — moleculeRenderer, hover for a larger depiction'),
    chips,
    h3('TypeAhead with structure rows'),
    picker,
    readout('picked', computed(() => picker.selected.value?.name ?? '(none)')),
  ], 'u2demo-page');
}

export function buildDemo(): Component {
  return Component.build(() => {
    const status = signal('Ready');
    const compact = signal(false);

    const bar = new Toolbar({ariaLabel: 'U2 demo toolbar'});
    bar.addButton('Run', () => status.value = 'Run clicked', {tooltip: 'Pretend to run'});
    bar.addButton('Save', () => status.value = 'Saved', {tooltip: 'Pretend to save'});
    bar.addSeparator();
    bar.addToggle('Compact', compact, {tooltip: 'Toolbar toggles bind to a signal'});
    bar.addMenu('Help', (m) => m
      .item('About u2', () => status.value = 'u2: next-gen Datagrok UI library')
      .item('Gallery', () => window.open('https://github.com/datagrok-ai/public/tree/master/libraries/u2')));

    const tabs = new TabStrip({tabs: [
      {id: 'inputs', label: 'Inputs', content: inputsPage()},
      {id: 'async', label: 'Combobox & async', content: asyncPage()},
      {id: 'lists', label: 'Lists & trees', content: listsPage()},
      {id: 'layout', label: 'Layout', content: layoutPage()},
      {id: 'popups', label: 'Popups', content: popupsPage(status)},
      {id: 'form', label: 'Form & properties', content: formPage()},
      {id: 'objectform', label: 'Object form', content: objectFormPage()},
      {id: 'entities', label: 'Entities', content: entitiesPage()},
      {id: 'spaces', label: 'Spaces & dashboards', content: spacesPage()},
      {id: 'molecules', label: 'Molecules', content: moleculesPage()},
      {id: 'dg', label: 'Platform bridge', content: platformPage()},
    ]});
    tabs.root.classList.add('u2demo-tabs');

    const statusBar = divH([
      span('Status: '), span(status),
      span(computed(() => compact.value ? ' · compact' : ''), 'u2demo-dim'),
      link('u2 sources', 'https://github.com/datagrok-ai/public/tree/master/libraries/u2'),
    ], 'u2demo-statusbar');

    return divV([bar, tabs, statusBar], 'u2demo-root');
  });
}
