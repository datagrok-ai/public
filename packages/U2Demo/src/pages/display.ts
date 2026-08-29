/* The Display group: the u2 controls that show rather than collect — cards and KPIs, the feedback
   surfaces (progress, toasts, the guided tour), the two table controls, the collapsible section
   and the stepped wizard, and the compose box. Anything that edits a value is an input and lives
   under Inputs (pages/inputs.ts). */
import {
  signal, computed, Scope,
  divV, divH, span, h3, button,
  Card, StatCard, ProgressBar, Tour, BasicTable, VirtualGrid,
  MultiSelect, Section, Wizard, TextInput, notify,
} from '@datagrok-libraries/u2';
import type {NotifyHandle} from '@datagrok-libraries/u2';
import {messageInput} from '@datagrok-libraries/u2/src/dg/index.js';
import {readout, DRUGS} from './common';

const PICKED = DRUGS.slice(0, 3);

export function cardsPage(): HTMLElement {
  const clicks = signal(0);
  const picked = PICKED.map(() => signal(false));
  const setAll = (value: boolean): void => {
    for (const p of picked)
      p.value = value;
  };

  const revenue = signal(1240000);
  const delta = signal(0.12);
  const move = (by: number, d: number): void => {
    revenue.value += by;
    delta.value = d;
  };

  return divV([
    span('One bordered surface: header, media, body and footer render only where they are given. ' +
      'clickable makes the whole card a button; selectable adopts the signal you hand it.',
    'u2demo-hint'),
    h3('Card'),
    divH([
      new Card({title: 'Aspirin', subtitle: 'NSAID', icon: 'capsules',
        children: [span('Acetylsalicylic acid, 180.16 g/mol.')],
        footer: [span('ChEMBL25', 'u2demo-dim')]}),
      new Card({title: 'Open the compound', subtitle: 'the whole card is a button',
        icon: 'external-link', children: [span('Click it, or focus it and press Enter.')],
        onClick: () => clicks.value++}),
    ], 'u2demo-cards'),
    readout('clicks', computed(() => String(clicks.value))),
    h3('Selectable cards'),
    divH(PICKED.map((drug, i) => new Card({title: drug.name, icon: 'flask',
      selectable: true, selected: picked[i], children: [span(drug.smiles, 'u2demo-code')]})),
    'u2demo-cards'),
    divH([button('Select all', () => setAll(true)), button('Clear', () => setAll(false))],
      'u2demo-row'),
    readout('selected', computed(() =>
      PICKED.filter((_, i) => picked[i].value).map((d) => d.name).join(', ') || '(none)')),
    h3('Stat cards'),
    divH([
      new StatCard({label: 'Revenue', value: revenue, delta, icon: 'chart-line',
        format: (x) => `${(Number(x) / 1e6).toFixed(2)}M`}),
      new StatCard({label: 'Error rate', value: '0.8%', delta: -0.05, deltaInverted: true,
        icon: 'exclamation-triangle'}),
      new StatCard({label: 'Queries today'}),
    ], 'u2demo-cards'),
    divH([button('+120k', () => move(120000, 0.12)), button('-120k', () => move(-120000, -0.09))],
      'u2demo-row'),
    readout('delta', computed(() => delta.value.toFixed(2))),
  ], 'u2demo-page');
}

export function feedbackPage(): HTMLElement {
  const scope = Scope.ambient!;

  const bar = new ProgressBar({showPercent: true, description: 'Indexing compounds'});
  const busy = new ProgressBar({indeterminate: true, description: 'Waiting for the server'});
  let timer = 0;
  const stop = (): void => {
    clearInterval(timer);
    timer = 0;
  };
  scope.own(stop);
  const run = (): void => {
    stop();
    bar.value.value = 0;
    timer = window.setInterval(() => {
      bar.value.value = Math.min(1, bar.value.peek() + 0.05);
      if (bar.value.peek() >= 1)
        stop();
    }, 120);
  };

  // a warning or an error balloon stays until it is closed, and it lives in the global overlay
  // host: leaving this page has to take the balloons this page raised with it
  const balloons: NotifyHandle[] = [];
  const toast = (handle: NotifyHandle | null): void => {
    if (handle)
      balloons.push(handle);
  };
  scope.own(() => {
    for (const handle of balloons)
      handle.close();
    balloons.length = 0;
  });

  const infoButton = button('Info', () => toast(notify.info('Indexed 1,204 compounds.')));
  const warnButton = button('Warning',
    () => toast(notify.warning('Two structures could not be parsed.')));
  const errorButton = button('Error', () => toast(notify.error('Connection refused.')));

  const outcome = signal('(not started)');
  let tour: Tour | undefined;
  // the tour is a document-level overlay, not a control: nothing disposes it with the page unless
  // the page says so, and a stuck dim layer would swallow every click in the app
  scope.own(() => tour?.finish());

  return divV([
    span('Progress, toasts and the guided tour — the three ways the library reports on itself. ' +
      'Each one is opt-in: nothing here starts by itself.', 'u2demo-hint'),
    h3('Progress'),
    bar,
    divH([button('Run', run), button('Reset', () => {
      stop();
      bar.value.value = 0;
    })], 'u2demo-row'),
    busy,
    readout('progress', computed(() => bar.value.value.toFixed(2))),
    h3('Notifications'),
    span('Info auto-hides; warnings and errors stay until closed. Leaving this sub-demo closes ' +
      'the ones it raised.', 'u2demo-hint'),
    divH([infoButton, warnButton, errorButton,
      button('Close all', () => notify.closeAll())], 'u2demo-row'),
    h3('Guided tour'),
    span('A walkthrough over the UI above: the dim layer is click-through, so a step can ask you ' +
      'to use the control it points at. Esc skips, Enter advances.', 'u2demo-hint'),
    divH([
      button('Start tour', () => {
        outcome.value = 'running';
        tour = Tour.run({
          steps: [
            {target: bar.root, content: 'The progress bar follows a 0..1 signal.'},
            {target: infoButton, position: 'top',
              content: 'Toasts come from one stack; click Info while the tour is up.'},
            {target: errorButton, content: 'Errors stay until you close them.'},
            {target: 'noSuchControl',
              content: 'You never see this step: its target does not exist, so the tour skips it.'},
          ],
          onFinish: (result) => outcome.value = result,
        });
      }, {primary: true}),
      button('Finish now', () => tour?.finish()),
    ], 'u2demo-row'),
    readout('tour', outcome),
  ], 'u2demo-page');
}

export function tablesPage(): HTMLElement {
  const table = new BasicTable<{name: string, smiles: string}>({
    columns: [
      {header: 'Name', render: (drug) => drug.name},
      {header: 'SMILES', render: (drug) => span(drug.smiles, 'u2demo-code')},
      {header: 'Atoms', render: (drug) => String(drug.smiles.length), align: 'right', width: '70px'},
    ],
    items: DRUGS,
    selectable: true,
  });

  const grid = new VirtualGrid<number>({
    cellWidth: 96, cellHeight: 44,
    keyOf: (item) => String(item),
    render: (item, index, cell) => {
      cell.title = `cell ${index}`;
      return span(`#${item}`);
    },
  });
  grid.setItems(Array.from({length: 20000}, (_, i) => i));
  grid.root.classList.add('u2demo-list');

  return divV([
    span('Two ways to show rows: a real table for data you can read whole, and a virtualized grid ' +
      'for anything bigger than the screen.', 'u2demo-hint'),
    h3('BasicTable'), table,
    readout('selectedIndex', computed(() => String(table.selectedIndex.value))),
    h3('VirtualGrid'),
    span('20,000 cells; the column count follows the pane width and only the visible cells exist.',
      'u2demo-hint'),
    grid,
    readout('columns', computed(() => String(grid.columns.value))),
    readout('selected cell', computed(() => String(grid.selectedIndex.value))),
  ], 'u2demo-page');
}

export function sectionsPage(): HTMLElement {
  const advanced = new Section({title: 'Advanced', expanded: false})
    .add(span('Collapsed content is hidden, never detached, so its state survives.'),
      new TextInput({label: 'Threshold', value: '0.85'}));

  const project = new TextInput({label: 'Project', placeholder: 'Required to continue'});
  // a step builds its own controls: a DOM node lives in one place, so a step that reused the
  // page's own control would move it out of the page and BACK would not bring it back
  const targets = new MultiSelect({label: 'Targets', selectAll: true, value: ['Kinase'],
    items: ['GPCR', 'Kinase', 'Ion channel', 'Protease', 'Transporter', 'Nuclear receptor']});
  const finished = signal('(not finished)');
  const wizard = new Wizard({
    steps: [
      {id: 'name', title: 'Name', content: () => divV([project]),
        canProceed: computed(() => project.value.value.trim() ? null : 'Enter a project name')},
      {id: 'data', title: 'Data', content: () =>
        divV([span('Lazily built on first visit; its state survives BACK.'), targets])},
      {id: 'review', title: 'Review', content: () => divV([span(computed(() =>
        `Creating "${project.value.value}" for ${targets.value.value.join(', ') || 'no target'}.`))])},
    ],
    onFinish: () => finished.value = 'finished',
  });

  return divV([
    span('Two ways to give a long form a shape: one collapsible section you can nest anywhere, ' +
      'and a wizard that spreads the same content over gated steps.', 'u2demo-hint'),
    h3('Section'), advanced,
    readout('expanded', computed(() => String(advanced.expanded.value))),
    h3('Wizard'),
    span('Panels are hidden, never detached, so a step\'s state survives BACK; a step marked ✓ is ' +
      'one you have completed, not merely one you have seen.', 'u2demo-hint'),
    wizard,
    readout('step', wizard.currentStep), readout('wizard', finished),
  ], 'u2demo-page');
}

export function messagingPage(): HTMLElement {
  const sent = signal('(nothing sent)');
  const message = messageInput({placeholder: 'Type @ to mention a user, Ctrl+Enter to send',
    onSend: (markup) => {
      sent.value = markup;
    }});

  return divV([
    span('A contenteditable whose value is Datagrok markup: a mention becomes an atomic chip that ' +
      'serializes back to the platform token.', 'u2demo-hint'),
    h3('MessageInput'), message,
    readout('sent', sent),
  ], 'u2demo-page');
}
