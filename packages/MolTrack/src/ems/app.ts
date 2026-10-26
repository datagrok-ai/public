import * as DG from 'datagrok-api/dg';
import {domains, DomainApp} from '@datagrok-libraries/ems/src/dg';
import type {CompoundBatchRow, CompoundRow} from '../generated/db';
import {moltrackEmsDb} from '../generated/db';
import {CompoundsApp} from './compounds';
import {MolTrackTables, SCHEMA} from './tables';

import '@datagrok-libraries/ems/src/dg/styles.js';
import '../../css/moltrack-ems.css';

export const APP_PATH = '/apps/Moltrack/MolTrackEMS';

const FILTERS = {filters: {mode: 'builder' as const, placement: 'left' as const}};

// a compound's batches are typed in its grid; the Grok Batch ID opens the batch's own page
const BATCH_GRID = {open: 'batch_regno', columns: ['batch_regno', 'source', 'source_batch_code', 'acquired_date',
  'synthesized_date', 'purity', 'project']};

const KEYS: Record<string, string> = {[`${SCHEMA}.compound`]: 'molregno', [`${SCHEMA}.compound_batch`]: 'batch_regno'};

const VIEWS = domains.views({
  path: APP_PATH,
  handle: (address) => MolTrackTables.table(address.slice(SCHEMA.length + 1)),
  ready: (app, table) => MolTrackEms.ready(app, table),
  views: [
    {table: `${SCHEMA}.compound`, name: 'Compounds', options: {...FILTERS, sort: '!molregno', app: CompoundsApp,
      mode: 'cards', itemHeight: 88, include: CompoundsApp.FIELDS,
      children: {tables: ['compound_batch'], grids: {compound_batch: BATCH_GRID}},
      // no Bulk edit: it would write past the identity a structure derives
      menu: {bulkEdit: false}}},
    {table: `${SCHEMA}.compound_batch`, name: 'Batches', options: {...FILTERS, sort: '!batch_regno',
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

  static async compoundByGrokId(corporateId: string): Promise<CompoundRow | null> {
    const molregno = MolTrackTables.keyNumber(await MolTrackTables.compounds(), 'molregno', corporateId);
    return molregno === null ? null : moltrackEmsDb.compounds.first({filter: DG.cond('molregno', '=', molregno)});
  }

  static async batchByGrokId(corporateId: string): Promise<CompoundBatchRow | null> {
    const regno = MolTrackTables.keyNumber(await MolTrackTables.table('compound_batch'), 'batch_regno', corporateId);
    return regno === null ? null : moltrackEmsDb.compoundBatches.first({filter: DG.cond('batch_regno', '=', regno)});
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
