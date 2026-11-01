import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {domains} from '@datagrok-libraries/ems/src/dg';
import type {DomainApp, DomainViewEntry} from '@datagrok-libraries/ems/src/dg';
import {moltrackDb} from '../generated/db';
import {Compounds} from './compounds';
import {MolTrackTables, SCHEMA} from './tables';
import {AppHeader} from '@datagrok-libraries/u2';
import type {RowView} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';
import {_package} from '../package';

import '@datagrok-libraries/ems/src/dg/styles.js';
import '../../css/moltrack.css';

export const APP_PATH = '/apps/Moltrack';

const FILTERS = {filters: {mode: 'builder' as const, placement: 'left' as const}};

const BATCH_GRID = {open: 'batch_regno', columns: ['batch_regno', 'source', 'source_batch_code', 'acquired_date',
  'synthesized_date', 'purity', 'project']};

const STRUCTURE = (row: RowView) => row.structure ? grok.chem.drawMolecule(String(row.structure), 120, 56) : '';
const COMPOUND_GRID = {rowHeight: 60, cells: {structure: STRUCTURE}};

const KEYS: Record<string, string> = {[`${SCHEMA}.compound`]: 'molregno', [`${SCHEMA}.compound_batch`]: 'batch_regno'};

const ENTRIES: DomainViewEntry[] = [
  {table: `${SCHEMA}.compound`, name: 'Compounds', options: {...FILTERS, sort: '!molregno',
    mode: 'cards', itemHeight: 88, exclude: Compounds.DERIVED, setup: Compounds.setup, menu: Compounds.menu,
    grid: COMPOUND_GRID,
    children: {tables: ['compound_batch'], grids: {compound_batch: BATCH_GRID}}}},
  {table: `${SCHEMA}.compound_batch`, name: 'Batches', options: {...FILTERS, sort: '!batch_regno',
    children: {tables: ['batch_addition', 'assay_result']}}},
  {table: `${SCHEMA}.assay`, name: 'Assays', options: {...FILTERS, sort: 'name', children: {tables: ['assay_run']}}},
  {table: `${SCHEMA}.assay_run`, name: 'Assay Runs', options: {...FILTERS, sort: '!run_date,name',
    children: {tables: ['assay_result']}}},
  {table: `${SCHEMA}.assay_result`, name: 'Assay Results',
    options: {...FILTERS, sort: '!created_on', children: false}},
  {table: `${SCHEMA}.addition`, name: 'Additions', options: {...FILTERS, sort: 'name'}},
];

const VIEWS = domains.views({
  path: APP_PATH,
  handle: (address) => MolTrackTables.table(address.slice(SCHEMA.length + 1)),
  ready: (app, table) => MolTrackEms.ready(app, table),
  views: ENTRIES,
});

export class MolTrackEms {
  static async open(address = ''): Promise<DG.ViewBase | null> {
    const target = await MolTrackEms._corporateLink(address) ?? address;
    return target.replace(/^\//, '') === '' ? MolTrackEms.home() : VIEWS.open(target);
  }

  static home(): DG.View {
    const rows = ENTRIES.map((entry) => {
      const count = ui.link('', () => void MolTrackEms.show(entry.name.replace(/ /g, '')));
      void grok.dapi.domains.table(entry.table).count()
        .then((n) => count.textContent = String(n), () => count.textContent = '?');
      return [entry.name, count];
    });
    const counts = ui.table(rows, (row) => row);
    counts.classList.add('moltrack-home-counts');
    const actions = MolTrackEms._quickActions([
      {label: 'Register compound', run: () => void VIEWS.create(`${SCHEMA}.compound`)},
      {label: 'Register batch', run: () => void VIEWS.create(`${SCHEMA}.compound_batch`)},
      {label: 'Register assay', run: () => void VIEWS.create(`${SCHEMA}.assay`)},
      {label: 'Register bulk', run: async () => void await domains.import(await MolTrackTables.compounds())},
      {label: 'Search', run: () => void MolTrackEms.show('Compounds')},
    ]);
    void domains.extensible(SCHEMA).then((tables) => {
      for (const [label, table] of [['Compound properties', 'compound'], ['Batch properties', 'compound_batch']]) {
        if (tables.has(table)) {
          actions.querySelector('.moltrack-quick-links')!.append(ui.link(label,
            () => void domains.editProperties(SCHEMA, table).catch((e) => grok.shell.error(e.message))));
        }
      }
    });
    const content = ui.divH([counts, actions]);
    content.classList.add('moltrack-home');
    const page = ui.divV([content]);
    const view = DG.View.fromRoot(page);
    page.prepend(ViewHost.of(view).build(() => MolTrackEms._header()));
    view.name = 'MolTrack';
    view.path = APP_PATH;
    return view;
  }

  static show(address: string): Promise<DG.ViewBase> {
    return VIEWS.show(address);
  }

  static treeBrowser(appNode: DG.TreeViewGroup): void {
    VIEWS.tree(appNode);
  }

  static ready(app: DomainApp, table: string): void {
    MolTrackTables.warmUp();
    const key = KEYS[table];
    if (key !== undefined)
      MolTrackEms._keySearch(app, key);
  }

  private static _header(): HTMLElement {
    return new AppHeader({
      icon: _package.webRoot + '/images/moltrack.png',
      learnMore: 'https://github.com/datagrok-ai/public/blob/master/packages/MolTrack/README.md',
      items: [
        'Chemical compound registration system',
        'Analyze assay data',
        'Find contextual information on molecules.',
      ],
      title: 'MolTrack',
      subtitle: 'Track, analyze, and manage chemical data',
      bottomLine: true,
    }).root;
  }

  private static _quickActions(links: {label: string; run: () => void}[]): HTMLElement {
    const boltIcon = ui.iconFA('bolt');
    boltIcon.classList.add('moltrack-bolt-icon');
    const quickActionsLabel = ui.label('Quick Actions');
    quickActionsLabel.classList.add('moltrack-quick-label');
    const quickActionsTitle = ui.divH([boltIcon, quickActionsLabel]);
    quickActionsTitle.classList.add('moltrack-quick-actions');
    const quickLinks = ui.divV(links.map((link) => ui.link(link.label, link.run)));
    quickLinks.classList.add('moltrack-quick-links');
    const container = ui.divV([quickActionsTitle, quickLinks]);
    container.classList.add('moltrack-quick-container');
    return container;
  }

  // the server cannot search a formatted key: a typed `DG-000060` becomes `molregno = 60`, only once the key is whole
  private static _keySearch(app: DomainApp, column: string): void {
    const {prefix, width} = MolTrackTables.keyFormat(app.table, column);
    app.effect(() => {
      const typed = app.listSource.search.value;
      const key = typed.trim().toUpperCase();
      const n = key.startsWith(prefix) && key.length >= prefix.length + width ?
        MolTrackTables.keyNumber(app.table, column, typed) : null;
      if (n === null)
        return;
      app.listSource.search.value = '';
      app.listSource.query.value = `${column} = ${n}`;
    });
  }

  private static async _corporateLink(address: string): Promise<string | null> {
    const params = new URLSearchParams(address.includes('?') ? address.slice(address.indexOf('?')) :
      window.location.search);
    const compound = params.get('corporate_compound_id');
    const batch = params.get('corporate_batch_id');
    if (compound !== null) {
      const n = MolTrackTables.keyNumber(await MolTrackTables.compounds(), 'molregno', compound);
      const row = n === null ? null :
        await moltrackDb.compounds.first({filter: DG.cond('molregno', '=', n), columns: ['molregno']});
      return row === null ? 'Compounds' : `Compounds?entity=${row.id}`;
    }
    if (batch !== null) {
      const n = MolTrackTables.keyNumber(await MolTrackTables.table('compound_batch'), 'batch_regno', batch);
      const row = n === null ? null :
        await moltrackDb.compoundBatches.first({filter: DG.cond('batch_regno', '=', n), columns: ['batch_regno']});
      return row === null ? 'Batches' : `Batches?entity=${row.id}`;
    }
    return null;
  }
}
