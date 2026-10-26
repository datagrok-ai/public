import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {domains, DomainTable} from '@datagrok-libraries/ems/src/dg';
import {Filters, Rows, divH, divV, isEmpty, span} from '@datagrok-libraries/u2';
import type {FilterGroup, RowView} from '@datagrok-libraries/u2';
import type {AssayResultRow, CompoundBatchRow, CompoundRow} from '../generated/db';
import {moltrackEmsDb} from '../generated/db';
import {Chemistry, HASH_COLUMNS, MATCH_RULES} from './chemistry';
import {BulkRegistration} from './registration';
import {StructureSearch} from './structure-search';
import type {MatchRule} from './chemistry';

export const SCHEMA = 'moltrack_ems';

const CAS = /^[1-9]\d{1,6}-\d{2}-\d$/;

type Resolve<T> = (value: T) => void;

// one handle per table, acquired once
export class MolTrackTables {
  private static readonly _handles = new Map<string, Promise<DomainTable<any>>>();
  private static readonly _names = new Map<string, Promise<string>>();
  private static _asked: Map<string, Resolve<string>> | null = null;

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

  /** A table's rows as a data frame, structure conditions answered in the browser first. */
  static async frame(name: string, root: FilterGroup, scope: Omit<DG.DomainQuerySpec, 'filter'> = {}):
    Promise<DG.DataFrame> {
    const resolved = await StructureSearch.resolve(root);
    const df = await grok.dapi.domains.table(`${SCHEMA}.${name}`).queryDf({...scope,
      ...(resolved.nodes.length === 0 ? {} : {filter: Filters.toDomainTree(resolved) as DG.DomainFilter})});
    df.name = (await MolTrackTables.table(name)).info.pluralName || name;
    return df;
  }

  static warmUp(): void {
    Chemistry.load().catch((e) => grok.log.warning(`MolTrack: Crux did not load: ${e}`));
  }

  private static _number(table: DomainTable, column: string, value: number): string {
    return DG.format(value, table.properties.find((p) => p.name === column)?.format || undefined);
  }

  /** Registered at load, so every handle over a table is configured alike. */
  static configure(): void {
    const id = {id: 'moltrack'};
    domains.configure<CompoundRow>(`${SCHEMA}.compound`, (table) => MolTrackTables._configureCompounds(table), id);
    domains.configure<CompoundBatchRow>(`${SCHEMA}.compound_batch`, async (table) =>
      MolTrackTables._configureBatches(table, await MolTrackTables.compounds()), id);
    domains.configure<AssayResultRow>(`${SCHEMA}.assay_result`, async (table) =>
      MolTrackTables._configureAssayResults(table, await MolTrackTables.table('compound_batch')), id);
  }

  private static _configureCompounds(table: DomainTable<CompoundRow>): void {
    table.importer = BulkRegistration.importer(table);
    table.validators.add('cas', (value) =>
      typeof value === 'string' && value !== '' && !CAS.test(value) ? 'CAS must look like 50-78-2' : null);
    const caption = (row: RowView<CompoundRow>) => Rows.isDraft(row) ? 'New compound' :
      MolTrackTables.corporateId(table, 'molregno', row.molregno);
    const card = (row: RowView<CompoundRow>) => divH([
      row.structure ? grok.chem.drawMolecule(row.structure, 120, 80) : span(''),
      divV([
        span(caption(row), 'u2-domain-card-title'),
        MolTrackTables._text(row.common_name ?? '', 'u2-domain-card-description'),
        span(row.formula ?? '', 'u2-domain-card-description'),
      ], 'moltrack-ems-card-text'),
    ], 'moltrack-ems-card');
    table.renderer = {...table.renderer, caption, card};
  }

  private static _configureBatches(table: DomainTable<CompoundBatchRow>, compounds: DomainTable<CompoundRow>): void {
    const caption = (row: RowView<CompoundBatchRow>) => Rows.isDraft(row) ? 'New batch' :
      MolTrackTables.corporateId(table, 'batch_regno', row.batch_regno);
    const listItem = (row: RowView<CompoundBatchRow>) => {
      const compound = MolTrackTables.corporateId(compounds, 'molregno', row[Rows.caption('compound_id')]);
      const rest = [row.project, row.source_batch_code].filter((d) => !isEmpty(d));
      const details = (name: string) => [[compound, name].filter((d) => d !== '').join(' '), ...rest]
        .filter((d) => d !== '').join(' · ');
      const line = MolTrackTables._text(details(''), 'u2-domain-list-details');
      if (!isEmpty(row.compound_id)) {
        void MolTrackTables.compoundName(row.compound_id).then((name) => {
          if (name !== '')
            line.textContent = line.title = details(name);
        });
      }
      return divH([span(caption(row), 'u2-domain-list-name'), line], 'u2-domain-list-line');
    };
    table.renderer = {...table.renderer, caption, listItem};
  }

  private static _configureAssayResults(table: DomainTable<AssayResultRow>, batches: DomainTable<CompoundBatchRow>):
    void {
    const batch = (row: RowView<AssayResultRow>) =>
      MolTrackTables.corporateId(batches, 'batch_regno', row[Rows.caption('batch_id')]);
    // the server captions a run by its key columns joined with ", "
    const run = (row: RowView<AssayResultRow>) =>
      String(row[Rows.caption('assay_run_id')] ?? '').replace(/, /g, ' · ');
    const caption = (row: RowView<AssayResultRow>) => Rows.isDraft(row) ? 'New assay result' :
      [batch(row), run(row)].filter((d) => d !== '').join(' · ');
    const clint = (row: RowView<AssayResultRow>) => isEmpty(row.reported_clint) ? '' :
      `CLint ${row.reported_clint_qualifier === '=' ? '' : row.reported_clint_qualifier ?? ''}` +
      MolTrackTables._number(table, 'reported_clint', row.reported_clint!);
    const dosed = (row: RowView<AssayResultRow>) => isEmpty(row.dosed_concentration) ? '' :
      `${MolTrackTables._number(table, 'dosed_concentration', row.dosed_concentration!)} uM`;
    table.renderer = {...table.renderer, caption, listItem: (row) => Rows.isDraft(row) ?
      span(caption(row), 'u2-domain-list-name') : divH([
        span(batch(row), 'u2-domain-list-name'),
        MolTrackTables._text(run(row), 'moltrack-ems-run'),
        MolTrackTables._text([dosed(row), clint(row)].filter((d) => d !== '').join(' · '), 'u2-domain-list-details'),
      ], 'u2-domain-list-line moltrack-ems-result')};
  }

  /** Text that may be cut, whole on hover. */
  private static _text(text: string, cls: string): HTMLElement {
    const el = span(text, cls);
    el.title = text;
    return el;
  }

  /** A compound's common name; names asked for together are read in one query. */
  static compoundName(id: string): Promise<string> {
    let name = MolTrackTables._names.get(id);
    if (name === undefined) {
      const asked = MolTrackTables._asked ??= MolTrackTables._readNames();
      name = new Promise((resolve) => asked.set(id, resolve));
      MolTrackTables._names.set(id, name);
    }
    return name;
  }

  private static _readNames(): Map<string, Resolve<string>> {
    const asked = new Map<string, Resolve<string>>();
    setTimeout(async () => {
      MolTrackTables._asked = null;
      const ids = [...asked.keys()];
      const rows = await moltrackEmsDb.compounds.query({filter: DG.cond('id', '=', ids), columns: ['common_name'],
        limit: ids.length}).catch(() => null);
      for (const [id, resolve] of asked) {
        if (rows === null)
          MolTrackTables._names.delete(id);
        resolve(rows?.find((r) => r.id === id)?.common_name ?? '');
      }
    });
    return asked;
  }
}

MolTrackTables.configure();
