import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DomainApp, DomainPick, DomainTable} from '@datagrok-libraries/ems/src/dg';
import type {ChildEntry} from '@datagrok-libraries/ems/src/dg';
import {Rows, isEmpty, span} from '@datagrok-libraries/u2';
import type {RowView} from '@datagrok-libraries/u2';
import type {CompoundBatchRow, CompoundRow} from '../generated/db';
import {moltrackEmsDb} from '../generated/db';
import {CompoundsApp} from './compounds-app';
import {Registration} from './registration';
import {MolTrackTables, SCHEMA} from './tables';

import '@datagrok-libraries/ems/src/dg/styles.js';
import '../../css/moltrack-ems.css';

export const APP_PATH = '/apps/Moltrack/MolTrackEMS';

/** A table's view: its child tabs — whose New opens the child's own create page with `childPages` — and the orders
 * its list offers, the first the one it opens in. */
type ViewEntry = {table: string, name: string, children?: string[], childPages?: boolean,
  sorts: [label: string, sort: string][]};

const NEWEST: [string, string][] = [['Newest first', '!created_on'], ['Oldest first', 'created_on']];

export const VIEWS: ViewEntry[] = [
  {table: 'compound', name: 'Compounds', children: ['compound_batch'], childPages: true,
    sorts: [['Newest first', '!molregno'], ['Oldest first', 'molregno'], ['By name', 'common_name,molregno']]},
  {table: 'compound_batch', name: 'Batches', children: ['batch_addition', 'assay_result'],
    sorts: [['Newest first', '!batch_regno'], ['Oldest first', 'batch_regno']]},
  {table: 'addition', name: 'Additions', sorts: [['By name', 'name'], ...NEWEST]},
  {table: 'assay', name: 'Assays', children: ['assay_run'], childPages: true, sorts: [['By name', 'name'], ...NEWEST]},
  {table: 'assay_run', name: 'Assay Runs', children: ['assay_result'],
    sorts: [['Latest run first', '!run_date,name'], ['Earliest run first', 'run_date,name'], ['By name', 'name']]},
  {table: 'assay_result', name: 'Assay Results', sorts: NEWEST},
];

export class MolTrackEms {
  static async view(table: string): Promise<DG.ViewBase> {
    MolTrackTables.warmUp();
    DomainPick.opener = MolTrackEms.openRow;
    const entry = VIEWS.find((v) => v.table === table)!;
    const handle = await MolTrackTables.table(table);
    return handle.app({
      name: entry.name, path: MolTrackEms.baseOf(entry),
      mode: table === 'compound' ? 'cards' : 'brief',
      filters: {mode: 'builder', placement: 'left'}, sorts: entry.sorts,
      ...(entry.children === undefined ? {children: false} :
        {children: {tables: await Promise.all(entry.children.map((child) => MolTrackTables.table(child))),
          ...(entry.childPages ? {onNew: (child: ChildEntry, values: Record<string, unknown>) =>
            void MolTrackEms.createIn(child.table.address, values)} : {})}}),
      // the server keeps both write-once: a compound's structure, and the compound a batch is of
      ...(table === 'compound' ? {app: CompoundsApp, include: [...CompoundsApp.FIELDS, ...CompoundsApp.DETAILS],
        immutable: ['original_structure'],
        sections: [{title: 'Registration details', include: CompoundsApp.DETAILS, expanded: false,
          hint: 'Registration hashes, matching rule, InChI'}],
        context: (row: RowView) => Rows.isDraft(row) ? null :
          MolTrackEms.compoundCard(handle, row as RowView<CompoundRow>)} : {}),
      ...(table === 'compound_batch' ? {immutable: ['compound_id'], context: (row: RowView) => {
        const id = (row as RowView<CompoundBatchRow>).compound_id;
        return isEmpty(id) ? null : ui.wait(async () =>
          MolTrackEms.compoundCard(await MolTrackTables.compounds(), await moltrackEmsDb.compounds.get(id)));
      }} : {}),
    });
  }

  static async register(df: DG.DataFrame, rejectAll: boolean): Promise<DG.DataFrame> {
    return await Registration.register(df, Registration.autoMap(df, await MolTrackTables.compounds()), {rejectAll});
  }

  /** A new view at an address (`Compounds?entity=new`, `Batches`, `Register` — the bulk registration over the
   * compounds) — what the app func answers the router, which docks it. The router hands over the path alone: a cold
   * deep link's query is replayed from the address bar by the view itself. */
  static async open(address = ''): Promise<DG.ViewBase> {
    const {entry, rest, bulk} = MolTrackEms._route(address);
    const view = await MolTrackEms.view(entry.table);
    void DomainApp.of(view)!.open(rest);
    if (bulk)
      MolTrackEms._import(view);
    return view;
  }

  /** An address from inside the app — the tree, a link: the open view comes forward at it, else a new one docks. */
  static async show(address: string): Promise<DG.ViewBase> {
    const {entry, rest, bulk} = MolTrackEms._route(address);
    const open = DomainApp.openAt(MolTrackEms.baseOf(entry), rest);
    if (open !== undefined && bulk)
      MolTrackEms._import(open);
    return open ?? grok.shell.addView(await MolTrackEms.open(address));
  }

  /** A row of a table with a view, opened there: what a reference's open control does when no view shows it. */
  static openRow(address: string, id: string): boolean {
    const entry = VIEWS.find((v) => `${SCHEMA}.${v.table}` === address);
    if (entry === undefined)
      return false;
    void MolTrackEms.show(`${MolTrackEms.segmentOf(entry)}?entity=${id}`);
    return true;
  }

  /** The create page of a table's view with `values` filled in — the open view's, else a new one's. */
  static async createIn(address: string, values: Record<string, unknown>): Promise<void> {
    const entry = VIEWS.find((v) => `${SCHEMA}.${v.table}` === address)!;
    if (!DomainApp.createAt(MolTrackEms.baseOf(entry), values))
      await DomainApp.of(await MolTrackEms.show(MolTrackEms.segmentOf(entry)))!.create(values);
  }

  static segmentOf(entry: ViewEntry): string {
    return entry.name.replace(/ /g, '');
  }

  static baseOf(entry: ViewEntry): string {
    return `${APP_PATH}/${MolTrackEms.segmentOf(entry)}`;
  }

  static async compoundPanel(corporateId: string): Promise<DG.Widget> {
    const digits = DomainPick.digits(corporateId);
    const compound = digits === '' ? null :
      await moltrackEmsDb.compounds.first({filter: `molregno = ${Number(digits)}`});
    return DG.Widget.fromRoot(compound === null ? ui.divText(`${corporateId} is not registered in MolTrack EMS`) :
      MolTrackEms.compoundCard(await MolTrackTables.compounds(), compound));
  }

  /** A compound at a glance, its Grok ID the way to its page — the Grok ID panel, and the context panel of a compound
   * and of a batch's compound; the batch count follows. */
  static compoundCard(compounds: DomainTable, compound: CompoundRow): HTMLElement {
    const corporateId = MolTrackTables.corporateId(compounds, 'molregno', compound.molregno);
    const open = ui.link(corporateId, () => void MolTrackEms.show(`Compounds?entity=${compound.id}`),
      'Open in MolTrack EMS');
    const batches = span('');
    void moltrackEmsDb.compoundBatches.count({filter: DG.cond('compound_id', '=', compound.id)})
      .then((count) => batches.textContent = String(count), () => batches.textContent = '?');
    return ui.divV([
      grok.chem.drawMolecule(compound.structure, 200, 150),
      ui.tableFromMap({'Grok ID': open, 'Name': compound.common_name ?? '', 'CAS': compound.cas ?? '',
        'Formula': compound.formula ?? '', 'InChIKey': compound.inchikey ?? '', 'Batches': batches}),
    ]);
  }

  static treeBrowser(appNode: DG.TreeViewGroup): void {
    const register = appNode.getOrCreateGroup('Register');
    // an assay is defined by the results it measures, which the schema fixes: Assays › New covers the rest
    for (const [name, address] of [['Compound', 'Compounds?entity=new'], ['Batch', 'Batches?entity=new'],
      ['Bulk', 'Register']])
      register.item(name).onSelected.subscribe(() => void MolTrackEms.show(address));
    const search = appNode.getOrCreateGroup('Search');
    for (const entry of VIEWS)
      search.item(entry.name).onSelected.subscribe(() => void MolTrackEms.show(MolTrackEms.segmentOf(entry)));
  }

  private static _route(address: string): {entry: ViewEntry, rest: string, bulk: boolean} {
    const [, segment, rest] = /^\/?([^/?]*)(.*)$/.exec(address)!;
    const bulk = segment === 'Register';
    return {entry: bulk ? VIEWS[0] : VIEWS.find((v) => MolTrackEms.segmentOf(v) === segment) ?? VIEWS[0], rest, bulk};
  }

  private static _import(view: DG.ViewBase): void {
    const app = DomainApp.of(view)!;
    void app.importRows().then((landed) => landed ? app.listSource.refresh() : null);
  }
}
