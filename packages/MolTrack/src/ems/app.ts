import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {domains, DomainApp} from '@datagrok-libraries/ems/src/dg';
import type {DomainAppOptions} from '@datagrok-libraries/ems/src/dg';
import {Filters, isEmpty, signal, span} from '@datagrok-libraries/u2';
import type {RowView} from '@datagrok-libraries/u2';
import type {CompoundBatchRow, CompoundRow} from '../generated/db';
import {moltrackEmsDb} from '../generated/db';
import {CompoundsApp} from './compounds';
import {StructureSearch} from './structure-search';
import {MolTrackTables, SCHEMA} from './tables';

import '@datagrok-libraries/ems/src/dg/styles.js';
import '../../css/moltrack-ems.css';

export const APP_PATH = '/apps/Moltrack/MolTrackEMS';

const FILTERS = {filters: {mode: 'builder' as const, placement: 'left' as const}};

const KEYS: Record<string, string> = {[`${SCHEMA}.compound`]: 'molregno', [`${SCHEMA}.compound_batch`]: 'batch_regno'};

// a batch page shows the compound it is a lot of, above its fields; the picker stays the way to set it
class BatchesApp extends DomainApp {
  private readonly _compound = ui.divV([], 'moltrack-ems-parent');
  /** The batch's compound, re-read on every edit; the card redraws only when it changes. */
  private readonly _compoundId = signal<string | null>(null);

  constructor(options: DomainAppOptions) {
    super(options);
    this.effect(() => {
      const source = this.entitySource.value;
      const edit = source?.edit.value;
      this._compoundId.value = BatchesApp._compoundOf(source?.currentRow.value ?? null);
      if (source === null || edit === undefined)
        return;
      const changes = edit.onChanged.subscribe(() =>
        this._compoundId.value = BatchesApp._compoundOf(source.currentRow.peek()));
      return () => changes.unsubscribe();
    });
    this.effect(() => {
      const id = this._compoundId.value;
      this._compound.remove();
      // a compound that does not load leaves the header without its card
      if (id !== null)
        this._showCompound(id).catch(() => {});
    });
  }

  private static _compoundOf(row: RowView | null): string | null {
    const id = (row as RowView<CompoundBatchRow> | null)?.compound_id;
    return isEmpty(id) ? null : id!;
  }

  private async _showCompound(id: string): Promise<void> {
    const compound = await moltrackEmsDb.compounds.first({filter: DG.cond('id', '=', id)});
    const card = compound === null ? null : await MolTrackEms.compoundCard(compound);
    if (this._compoundId.peek() !== id)
      return;
    this._compound.replaceChildren(ui.h3('Compound'), card ?? ui.divText('The compound was not found'));
    this.header.append(this._compound);
  }
}

const VIEWS = domains.views({
  path: APP_PATH,
  handle: (address) => MolTrackTables.table(address.slice(SCHEMA.length + 1)),
  ready: (app, table) => MolTrackEms.ready(app, table),
  views: [
    {table: `${SCHEMA}.compound`, name: 'Compounds', options: {...FILTERS, sort: '!molregno', app: CompoundsApp,
      mode: 'cards', itemHeight: 88, include: CompoundsApp.FIELDS, children: {tables: ['compound_batch'], mode: 'list'},
      // no Bulk edit: it would write past the identity a structure derives
      menu: {bulkEdit: false}}},
    {table: `${SCHEMA}.compound_batch`, name: 'Batches', options: {...FILTERS, sort: '!batch_regno', app: BatchesApp,
      children: {tables: ['batch_addition', 'assay_result']}}},
    {table: `${SCHEMA}.assay`, name: 'Assays', options: {...FILTERS, sort: 'name', children: {tables: ['assay_run']}}},
    {table: `${SCHEMA}.assay_run`, name: 'Assay Runs', options: {...FILTERS, sort: '!run_date,name',
      children: {tables: ['assay_result']}}},
    {table: `${SCHEMA}.assay_result`, name: 'Assay Results',
      options: {...FILTERS, sort: '!created_on', children: false}},
    {table: `${SCHEMA}.addition`, name: 'Additions', options: {...FILTERS, sort: 'name'}},
  ],
});

export class MolTrackEms {
  static async open(address = ''): Promise<DG.ViewBase> {
    return VIEWS.open(await MolTrackEms._corporateLink(address) ?? address);
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

  static async search(table: string, filter: string, structure: string, operator: string): Promise<DG.DataFrame> {
    const parsed = Filters.parse(filter);
    if (parsed.problems.length > 0)
      throw new Error(`${filter}: ${parsed.problems[0].message}`);
    if (structure !== '' && (table !== 'compound' || !StructureSearch.OPERATORS.includes(operator)))
      throw new Error(`"${operator}" does not search ${table}`);
    const root = structure === '' ? parsed.root :
      Filters.group('and', [parsed.root, Filters.cond('structure', operator, structure)]);
    return MolTrackTables.frame(table, root);
  }

  static async compoundByGrokId(corporateId: string): Promise<CompoundRow | null> {
    const molregno = MolTrackTables.keyNumber(await MolTrackTables.compounds(), 'molregno', corporateId);
    return molregno === null ? null : moltrackEmsDb.compounds.first({filter: DG.cond('molregno', '=', molregno)});
  }

  static async batchByGrokId(corporateId: string): Promise<CompoundBatchRow | null> {
    const regno = MolTrackTables.keyNumber(await MolTrackTables.table('compound_batch'), 'batch_regno', corporateId);
    return regno === null ? null : moltrackEmsDb.compoundBatches.first({filter: DG.cond('batch_regno', '=', regno)});
  }

  static async compoundPanel(corporateId: string): Promise<DG.Widget> {
    const compound = await MolTrackEms.compoundByGrokId(corporateId);
    return DG.Widget.fromRoot(compound === null ? ui.divText(`${corporateId} is not registered`) :
      await MolTrackEms.compoundCard(compound));
  }

  static async batchPanel(corporateId: string): Promise<DG.Widget> {
    const batch = await MolTrackEms.batchByGrokId(corporateId);
    if (batch === null)
      return DG.Widget.fromRoot(ui.divText(`${corporateId} is not registered`));
    const batches = await MolTrackTables.table('compound_batch');
    const id = MolTrackTables.corporateId(batches, 'batch_regno', batch.batch_regno);
    const compound = await moltrackEmsDb.compounds.first({filter: DG.cond('id', '=', batch.compound_id)});
    return DG.Widget.fromRoot(ui.divV([
      ui.tableFromMap({'Batch': ui.link(id, () => void MolTrackEms.show(`Batches?entity=${batch.id}`)),
        'Project': batch.project ?? '', 'Source batch code': batch.source_batch_code ?? '',
        'Responsible party': batch.responsible_party ?? ''}),
      ...(compound === null ? [] : [ui.h3('Compound'), await MolTrackEms.compoundCard(compound)]),
    ]));
  }

  static async compoundCard(compound: CompoundRow): Promise<HTMLElement> {
    const id = MolTrackTables.corporateId(await MolTrackTables.compounds(), 'molregno', compound.molregno);
    const batches = span('');
    void moltrackEmsDb.compoundBatches.count({filter: DG.cond('compound_id', '=', compound.id)})
      .then((count) => batches.textContent = String(count), () => batches.textContent = '?');
    return ui.divV([
      grok.chem.drawMolecule(compound.structure, 200, 150),
      ui.tableFromMap({'Grok ID': ui.link(id, () => void MolTrackEms.show(`Compounds?entity=${compound.id}`)),
        'Name': compound.common_name ?? '', 'CAS': compound.cas ?? '', 'Formula': compound.formula ?? '',
        'InChIKey': compound.inchikey ?? '', 'Batches': batches}),
    ]);
  }

  // the server cannot search a formatted key: a typed `DG-000060` becomes `molregno = 60`
  private static _keySearch(app: DomainApp, column: string): void {
    const {prefix} = MolTrackTables.keyFormat(app.table, column);
    app.effect(() => {
      const typed = app.listSource.search.value;
      const n = typed.trim().toUpperCase().startsWith(prefix) ?
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
      const row = await MolTrackEms.compoundByGrokId(compound);
      return row === null ? 'Compounds' : `Compounds?entity=${row.id}`;
    }
    if (batch !== null) {
      const row = await MolTrackEms.batchByGrokId(batch);
      return row === null ? 'Batches' : `Batches?entity=${row.id}`;
    }
    return null;
  }
}
