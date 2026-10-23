import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {domains, DomainTable} from '@datagrok-libraries/ems/src/dg';
import {Rows, divH, divV, isEmpty, span} from '@datagrok-libraries/u2';
import type {RowView} from '@datagrok-libraries/u2';
import type {AssayResultRow, CompoundBatchRow, CompoundRow} from '../generated/db';
import {moltrackEmsDb} from '../generated/db';
import {Chemistry, HASH_COLUMNS, MATCH_RULES} from './chemistry';
import {CompoundsApp} from './compounds';
import {BulkRegistration} from './registration';
import type {MatchRule} from './chemistry';

export const SCHEMA = 'moltrack_ems';

const CAS = /^[1-9]\d{1,6}-\d{2}-\d$/;

// one handle per table, acquired once; what MolTrack declares on a table is configured on every handle over it
export class MolTrackTables {
  private static readonly _handles = new Map<string, Promise<DomainTable<any>>>();

  static matchRule(): MatchRule {
    const rule = Chemistry.package().settings?.['CompoundMatchingRule'];
    return MATCH_RULES.includes(rule) ? rule : 'ALL_LAYERS';
  }

  static keyFormat(table: DomainTable, column: string): {prefix: string, width: number} {
    const format = table.properties.find((p) => p.name === column)?.format ?? '';
    const width = /0+$/.exec(format)?.[0].length ?? 0;
    return {prefix: format.slice(0, format.length - width), width};
  }

  static keyNumber(table: DomainTable, column: string, text: string): number | null {
    const {prefix} = MolTrackTables.keyFormat(table, column);
    const typed = text.trim().toUpperCase();
    const digits = typed.startsWith(prefix) ? typed.slice(prefix.length) : typed;
    return /^\d+$/.test(digits) ? Number(digits) : null;
  }

  // the server captions a reference by the bare number
  static corporateId(table: DomainTable, column: string, value: unknown): string {
    const {prefix, width} = MolTrackTables.keyFormat(table, column);
    return isEmpty(value) ? '' : `${prefix}${String(value).padStart(width, '0')}`;
  }

  static registeredAs(hash: string, rule: MatchRule): Promise<CompoundRow | null> {
    return moltrackEmsDb.compounds.first({filter: DG.cond(HASH_COLUMNS[rule], '=', hash),
      columns: ['molregno', 'common_name'], sort: 'molregno'});
  }

  static table(name: string): Promise<DomainTable<any>> {
    let handle = MolTrackTables._handles.get(name);
    if (handle === undefined) {
      handle = domains.table<any>(`${SCHEMA}.${name}`);
      handle.catch(() => MolTrackTables._handles.delete(name));
      MolTrackTables._handles.set(name, handle);
    }
    return handle;
  }

  static compounds(): Promise<DomainTable<CompoundRow>> {
    return MolTrackTables.table('compound');
  }

  static warmUp(): void {
    Chemistry.load().catch((e) => grok.log.warning(`MolTrack: Crux did not load: ${e}`));
  }

  // a number as the column's own format shows it, the platform's default where it has none
  private static _number(table: DomainTable, column: string, value: number): string {
    return DG.format(value, table.properties.find((p) => p.name === column)?.format || undefined);
  }

  /** Registered at load, so a view's handle, a child tab's and the `/domains` route's are all dressed alike. */
  static configure(): void {
    const id = {id: 'moltrack'};
    domains.configure<CompoundRow>(`${SCHEMA}.compound`, (table) => MolTrackTables._configureCompounds(table), id);
    domains.configure<CompoundBatchRow>(`${SCHEMA}.compound_batch`, async (table) =>
      MolTrackTables._configureBatches(table, await MolTrackTables.compounds()), id);
    domains.configure<AssayResultRow>(`${SCHEMA}.assay_result`, async (table) =>
      MolTrackTables._configureAssayResults(table, await MolTrackTables.table('compound_batch')), id);
  }

  private static _configureCompounds(table: DomainTable<CompoundRow>): void {
    CompoundsApp.configure(table);
    table.importer = BulkRegistration.importer(table);
    table.validators.add('cas', (value) =>
      typeof value === 'string' && value !== '' && !CAS.test(value) ? 'CAS must look like 50-78-2' : null);
    const caption = (row: RowView<CompoundRow>) => Rows.isDraft(row) ? 'New compound' :
      MolTrackTables.corporateId(table, 'molregno', row.molregno);
    const card = (row: RowView<CompoundRow>) => divH([
      row.structure ? grok.chem.drawMolecule(row.structure, 84, 56) : span(''),
      divV([
        span(caption(row), 'u2-domain-card-title'),
        span(row.common_name ?? '', 'u2-domain-card-description'),
        span(row.formula ?? '', 'u2-domain-card-description'),
      ], 'moltrack-ems-card-text'),
    ], 'moltrack-ems-card');
    table.renderer = {...table.renderer, caption, card};
  }

  private static _configureBatches(table: DomainTable<CompoundBatchRow>, compounds: DomainTable<CompoundRow>): void {
    const caption = (row: RowView<CompoundBatchRow>) => Rows.isDraft(row) ? 'New batch' :
      MolTrackTables.corporateId(table, 'batch_regno', row.batch_regno);
    const compound = (row: RowView<CompoundBatchRow>) =>
      MolTrackTables.corporateId(compounds, 'molregno', row[Rows.caption('compound_id')]);
    table.renderer = {...table.renderer, caption, listItem: (row) => divH([span(caption(row), 'u2-domain-list-name'),
      span([compound(row), row.project, row.source_batch_code].filter((d) => !isEmpty(d)).join(' · '),
        'u2-domain-card-description')], 'moltrack-ems-line')};
  }

  private static _configureAssayResults(table: DomainTable<AssayResultRow>, batches: DomainTable<CompoundBatchRow>):
    void {
    const caption = (row: RowView<AssayResultRow>) => Rows.isDraft(row) ? 'New assay result' :
      [MolTrackTables.corporateId(batches, 'batch_regno', row[Rows.caption('batch_id')]),
        row[Rows.caption('assay_run_id')]].filter((d) => !isEmpty(d)).join(' · ');
    const clint = (row: RowView<AssayResultRow>) => isEmpty(row.reported_clint) ? '' :
      `CLint ${row.reported_clint_qualifier === '=' ? '' : row.reported_clint_qualifier ?? ''}` +
      MolTrackTables._number(table, 'reported_clint', row.reported_clint!);
    const dosed = (row: RowView<AssayResultRow>) => isEmpty(row.dosed_concentration) ? '' :
      `${MolTrackTables._number(table, 'dosed_concentration', row.dosed_concentration!)} uM`;
    table.renderer = {...table.renderer, caption, listItem: (row) => divH([span(caption(row), 'u2-domain-list-name'),
      span([dosed(row), clint(row)].filter((d) => d !== '').join(' · '), 'u2-domain-card-description'),
    ], 'moltrack-ems-line')};
  }
}

MolTrackTables.configure();
