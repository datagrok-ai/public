import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DomainApp} from '@datagrok-libraries/ems/src/dg';
import type {DomainAppViewOptions} from '@datagrok-libraries/ems/src/dg';
import {span} from '@datagrok-libraries/u2';
import type {CompoundRow} from '../generated/db';
import {moltrackEmsDb} from '../generated/db';
import {CompoundsApp} from './compounds';
import {MolTrackTables} from './tables';

import '@datagrok-libraries/ems/src/dg/styles.js';
import '../../css/moltrack-ems.css';

export const APP_PATH = '/apps/Moltrack/MolTrackEMS';

interface ViewEntry {
  table: string;
  name: string;
  sort: string;
  options: DomainAppViewOptions;
}

const KEYS: Record<string, string> = {compound: 'molregno', compound_batch: 'batch_regno'};

const VIEWS: ViewEntry[] = [
  {table: 'compound', name: 'Compounds', sort: '!molregno', options: {app: CompoundsApp, mode: 'cards',
    include: CompoundsApp.FIELDS, children: {tables: ['compound_batch'], mode: 'list'}}},
  {table: 'compound_batch', name: 'Batches', sort: '!batch_regno',
    options: {children: {tables: ['batch_addition', 'assay_result']}}},
  {table: 'assay', name: 'Assays', sort: 'name', options: {children: {tables: ['assay_run']}}},
  {table: 'assay_run', name: 'Assay Runs', sort: '!run_date,name', options: {children: {tables: ['assay_result']}}},
  {table: 'assay_result', name: 'Assay Results', sort: '!created_on', options: {children: false}},
  {table: 'addition', name: 'Additions', sort: 'name', options: {}},
];

export class MolTrackEms {
  // the router passes the path alone: a cold deep link's query is in the address bar, when it is this view's
  static async open(address = ''): Promise<DG.ViewBase> {
    const {entry, rest} = MolTrackEms._parse(await MolTrackEms._corporateLink(address) ?? address);
    const {view, app} = await MolTrackEms._create(entry);
    const own = window.location.pathname.toLowerCase().endsWith(MolTrackEms._base(entry).toLowerCase());
    void app.open(rest || (own ? undefined : ''));
    return view;
  }

  static async show(address: string): Promise<DG.ViewBase> {
    const {entry, rest} = MolTrackEms._parse(address);
    const open = DomainApp.openAt(MolTrackEms._base(entry), rest);
    if (open !== undefined)
      return open;
    const {view, app} = await MolTrackEms._create(entry);
    grok.shell.addView(view);
    void app.open(rest);
    return view;
  }

  static treeBrowser(appNode: DG.TreeViewGroup): void {
    for (const entry of VIEWS)
      appNode.item(entry.name).onSelected.subscribe(() => void MolTrackEms.show(MolTrackEms._segment(entry)));
  }

  static async compoundPanel(corporateId: string): Promise<DG.Widget> {
    const molregno = MolTrackTables.keyNumber(await MolTrackTables.compounds(), 'molregno', corporateId);
    const compound = molregno === null ? null :
      await moltrackEmsDb.compounds.first({filter: DG.cond('molregno', '=', molregno)});
    return DG.Widget.fromRoot(compound === null ? ui.divText(`${corporateId} is not registered`) :
      await MolTrackEms.compoundCard(compound));
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

  private static async _create(entry: ViewEntry): Promise<{view: DG.ViewBase, app: DomainApp}> {
    MolTrackTables.warmUp();
    const view = (await MolTrackTables.table(entry.table)).app({filters: {mode: 'builder', placement: 'left'},
      ...entry.options, name: entry.name, path: MolTrackEms._base(entry)});
    const app = DomainApp.of(view)!;
    app.listSource.sort.value = entry.sort;
    const key = KEYS[entry.table];
    if (key !== undefined)
      MolTrackEms._keySearch(app, key);
    return {view, app};
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
      const n = MolTrackTables.keyNumber(await MolTrackTables.compounds(), 'molregno', compound);
      const row = n === null ? null :
        await moltrackEmsDb.compounds.first({filter: DG.cond('molregno', '=', n), columns: ['molregno']});
      return row === null ? 'Compounds' : `Compounds?entity=${row.id}`;
    }
    if (batch !== null) {
      const n = MolTrackTables.keyNumber(await MolTrackTables.table('compound_batch'), 'batch_regno', batch);
      const row = n === null ? null :
        await moltrackEmsDb.compoundBatches.first({filter: DG.cond('batch_regno', '=', n), columns: ['batch_regno']});
      return row === null ? 'Batches' : `Batches?entity=${row.id}`;
    }
    return null;
  }

  private static _segment(entry: ViewEntry): string {
    return entry.name.replace(/\s/g, '');
  }

  private static _base(entry: ViewEntry): string {
    return `${APP_PATH}/${MolTrackEms._segment(entry)}`;
  }

  private static _parse(address: string): {entry: ViewEntry, rest: string} {
    const [, segment = '', rest = ''] = /^\/?([^/?]*)(.*)$/.exec(address) ?? [];
    const entry = VIEWS.find((e) => MolTrackEms._segment(e).toLowerCase() === segment.toLowerCase());
    return entry === undefined ? {entry: VIEWS[0], rest: address.startsWith('?') ? address : ''} : {entry, rest};
  }
}
