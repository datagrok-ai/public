import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {
  signal, computed, Scope,
  divV, divH, span, h3, button,
  TextInput, TypeAhead,
} from '@datagrok-libraries/u2';
import {asDartInput, tableInput, columnInput, leakReport, propertyForm, IProperty,
  userInput, entityInput, chip, EntityChip, entityCard, handlerRenderer,
  moleculeInput, moleculeRenderer, fileInput, filesInput, rsaInput, tablesInput,
  columnsInput, columnsMapInput, aggregatedColumnsInput}
  from '@datagrok-libraries/u2/src/dg/index.js';
import {readout, DRUGS} from './common';

export function filesPage(): HTMLElement {
  const file = fileInput('File');
  const files = filesInput('Files');
  const key = rsaInput('Private key');

  return divV([
    span('Inputs that need the platform: file shares and keys through dapi — the pickers, not the ' +
      'file manager (that is Browse ▸ Files). On a stand a typed path resolves against the real ' +
      'shares; in local mode (login.html?mode=local) dapi answers every existence check with ' +
      'nothing, so a path settles on "does not exist" — the state the machine is supposed to reach.',
    'u2demo-hint'),
    h3('File and key pickers'),
    file, readout('file', computed(() => file.value.value?.name ?? '(none)')),
    files, readout('files', computed(() =>
      files.value.value.map((f) => f.name).join(', ') || '(none)')),
    key, readout('key', computed(() => {
      const v = key.value.value;
      return v ? `${v.length} characters` : '(none)';
    })),
  ], 'u2demo-page');
}

export function dataframesPage(): HTMLElement {
  const demog = grok.data.demo.demog(100);

  const table = tableInput('Open table');
  const demogColumn = columnInput('Demog column', demog);
  const column = columnInput('Column', demog, {rich: true});
  const tables = tablesInput('Tables');
  const columns = columnsInput('Columns', demog);
  const mapping = columnsMapInput('Mapping', ['x', 'y'], demog);
  const aggregations = aggregatedColumnsInput('Aggregations', demog);

  return divV([
    h3('Tables and columns'),
    table, readout('table', table.value),
    demogColumn, readout('demog column', demogColumn.value),
    column, readout('column', computed(() => column.value.value ?? '(none)')),
    tables, readout('tables', computed(() => tables.value.value.join(', ') || '(none)')),
    columns, readout('columns', computed(() => columns.value.value.join(', ') || '(none)')),
    mapping, readout('mapping', computed(() => JSON.stringify(mapping.value.value))),
    aggregations, readout('aggregations', computed(() => JSON.stringify(aggregations.value.value))),
  ], 'u2demo-page');
}

export function entitiesPage(): HTMLElement {
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

export function spacesPage(): HTMLElement {
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

export function moleculesPage(): HTMLElement {
  // everything here goes through Chem: the sketcher input, the semType editor and drawMolecule.
  // Without the package `grok.chem.drawMolecule` still returns a host and then fails to resolve
  // its asset — five console errors on app open in local mode, one per depiction
  if (DG.Func.find({package: 'Chem', name: 'drawMolecule'}).length === 0) {
    return divV([span('This sub-demo needs the Chem package: the sketcher input, the semType-driven ' +
      'editor and the structure depictions are all bridged to it. Open the app on a stand that ' +
      'has Chem published.', 'u2demo-hint')], 'u2demo-page');
  }
  const structure = moleculeInput('Structure');
  structure.value.value = DRUGS[0].smiles;

  const compound = {name: 'Aspirin', smiles: DRUGS[0].smiles, mw: 180.16};
  const props: IProperty[] = [
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

export function bridgePage(): HTMLElement {
  const bridged = new TextInput({label: 'Compound', value: 'Aspirin'});
  bridged.addValidator((v) => v.trim() ? null : 'Required');
  const dart = asDartInput(bridged);

  const leaks = signal('(press the button)');

  return divV([
    span('The u2 core lives in datagrok-api (DG.U2): every DG.Widget is a u2 Control, ' +
      'asDartInput() makes any u2 input a DG.InputBase.', 'u2demo-hint'),
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
