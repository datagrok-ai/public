import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {DomainQueries} from '@datagrok-libraries/ems';
import {domains, DomainSkip, DomainTable} from '@datagrok-libraries/ems/src/dg';
import {moleculeRenderer} from '@datagrok-libraries/u2/src/dg/index.js';
import {Rows, divH, divV, isEmpty, span, text} from '@datagrok-libraries/u2';
import type {RowView} from '@datagrok-libraries/u2';
import type {AdditionRow, AssayResultRow, CompoundRow, CompoundBatchRow, MoltrackEmsTransactionOp}
  from '../generated/db';
import {moltrackEmsDb} from '../generated/db';
import {Chemistry, HASH_COLUMNS, HashColumn, MatchRule, MATCH_RULES} from './chemistry';
import {StructureSearch} from './structure-search';

export const SCHEMA = 'moltrack_ems';

const CAS = /^[1-9]\d{1,6}-\d{2}-\d$/;

export interface Registered {
  id: string;
  molregno: number;
}

export class MolTrackTables {
  private static readonly _handles = new Map<string, Promise<DomainTable<any>>>();
  private static readonly _batches = new Map<string, {keys: Set<string>, done: Promise<Map<string, unknown>>}>();

  static matchRule(): MatchRule {
    const rule = Chemistry.package().settings?.['CompoundMatchingRule'];
    return MATCH_RULES.includes(rule) ? rule : 'ALL_LAYERS';
  }

  static corporateId(table: DomainTable, column: string, value: unknown): string {
    return isEmpty(value) ? '' : DomainTable.captionText(table.properties, column, Number(value));
  }

  static alreadyRegistered(corporateId: string): string {
    return `This structure is already registered as ${corporateId}`;
  }

  /** The earliest compound registered under each of `hashes` in `column`. */
  static async registered(column: HashColumn | 'hash_mol', hashes: string[]): Promise<Map<string, Registered>> {
    const rows = await DomainQueries.chunked(hashes, (chunk) => DomainQueries.all((page) =>
      moltrackEmsDb.compounds.query({filter: DG.cond(column, '=', chunk), columns: [column, 'molregno'],
        sort: 'molregno', ...page}), DomainQueries.MAX_ROWS));
    const found = new Map<string, Registered>();
    for (const row of rows) {
      if (!found.has(row[column]!))
        found.set(row[column]!, {id: row.id, molregno: row.molregno});
    }
    return found;
  }

  // batched: a bulk run derives every row at once, and its lookups go out as one query per chunk
  static lookup(column: HashColumn, hash: string): Promise<Registered | undefined> {
    return MolTrackTables._batched(column, hash, (hashes) => MolTrackTables.registered(column, hashes));
  }

  /** One configured handle per table for the page: every caller shares its derivation, validators and renderer. */
  static table(name: string): Promise<DomainTable<any>> {
    let handle = MolTrackTables._handles.get(name);
    if (handle === undefined) {
      handle = MolTrackTables._open(name).catch((e) => {
        MolTrackTables._handles.delete(name);
        throw e;
      });
      MolTrackTables._handles.set(name, handle);
    }
    return handle;
  }

  static compounds(): Promise<DomainTable<CompoundRow>> {
    return MolTrackTables.table('compound');
  }

  static batches(): Promise<DomainTable<CompoundBatchRow>> {
    return MolTrackTables.table('compound_batch');
  }

  static line(caption: string, ...details: unknown[]): HTMLElement {
    return divH([span(caption, 'u2-domain-list-name'),
      span(details.filter((d) => !isEmpty(d)).join(' · '), 'u2-domain-card-description')], 'moltrack-ems-line');
  }

  /** `key` looked up together with every key asked for under `group` in the same tick: one query per tick. */
  private static _batched<T>(group: string, key: string, run: (keys: string[]) => Promise<Map<string, T>>):
    Promise<T | undefined> {
    let batch = MolTrackTables._batches.get(group);
    if (batch === undefined) {
      const keys = new Set<string>();
      const done = new Promise((resolve) => setTimeout(resolve)).then(() => {
        MolTrackTables._batches.delete(group);
        return run([...keys]);
      });
      MolTrackTables._batches.set(group, batch = {keys, done});
    }
    batch.keys.add(key);
    return batch.done.then((found) => found.get(key) as T | undefined);
  }

  // `hash_mol` is the identity a compound was registered under: counted where it would differ, never rewritten
  static async reidentify(): Promise<string> {
    const rows = await DomainQueries.all((page) => moltrackEmsDb.compounds.query({
      columns: ['original_structure', 'structure', 'match_rule', 'hash_mol'], sort: 'molregno', ...page}));
    let moved = 0;
    await DomainQueries.chunked(rows, async (chunk) => {
      const ops: MoltrackEmsTransactionOp[] = [];
      for (const row of chunk) {
        const identity =
          await Chemistry.identify(row.original_structure ?? row.structure, row.match_rule ?? 'ALL_LAYERS');
        ops.push({op: 'update', table: 'compound', id: row.id, values: {hash_all_layers: identity.hash_all_layers,
          hash_stereo_insensitive: identity.hash_stereo_insensitive,
          hash_tautomer_insensitive: identity.hash_tautomer_insensitive}});
        if (identity.hash_mol !== row.hash_mol)
          moved++;
      }
      return moltrackEmsDb.transaction(ops);
    });
    return `${rows.length} compounds re-identified; ${moved} would now have a different registration hash`;
  }

  static warmUp(): void {
    Chemistry.load().catch((e) => grok.log.warning(`MolTrack: Crux did not load: ${e}`));
  }

  private static async _open(name: string): Promise<DomainTable<any>> {
    const table = await domains.table<any>(`${SCHEMA}.${name}`);
    if (name === 'compound')
      MolTrackTables._configureCompounds(table);
    else if (name === 'compound_batch')
      MolTrackTables._configureBatches(table);
    else if (name === 'addition')
      MolTrackTables._configureAdditions(table);
    else if (name === 'assay_result')
      MolTrackTables._configureAssayResults(table);
    return table;
  }

  private static _configureCompounds(table: DomainTable<CompoundRow>): void {
    table.derived.add({
      from: ['original_structure'],
      to: ['structure', 'hash_mol', 'match_rule', 'hash_all_layers', 'hash_stereo_insensitive',
        'hash_tautomer_insensitive', 'inchi', 'inchikey', 'formula'],
      unique: 'hash_mol',
      derive: async (row) => {
        const rule = MolTrackTables.matchRule();
        const identity = await Chemistry.identify(row.original_structure!, rule);
        // a saved compound is the one registered under its hash; an imported row has no id yet
        const existing = row.id !== undefined && !Rows.isDraft(row) ? undefined :
          await MolTrackTables.lookup(HASH_COLUMNS[rule], identity.hash_mol);
        if (existing !== undefined) {
          throw new DomainSkip(
            MolTrackTables.alreadyRegistered(MolTrackTables.corporateId(table, 'molregno', existing.molregno)));
        }
        return {...identity};
      },
    });
    table.validators.add('cas', (value) =>
      typeof value === 'string' && value !== '' && !CAS.test(value) ? 'CAS must look like 50-78-2' : null);
    StructureSearch.register(table);
    const thumbnail = moleculeRenderer({width: 84, height: 56});
    const caption = table.renderer.caption;
    const card = (row: RowView<CompoundRow>) => divH([
      row.structure ? thumbnail.listItem!(row.structure) : span(''),
      divV([
        span(caption(row), 'u2-domain-card-title'),
        span(row.common_name ?? '', 'u2-domain-card-description'),
        span(row.formula ?? '', 'u2-domain-card-description'),
      ], 'moltrack-ems-card-text'),
    ], 'moltrack-ems-card');
    table.renderer = {...table.renderer, card,
      tooltip: (row) => row.structure ? thumbnail.tooltip!(row.structure) : caption(row)};
  }

  // a reference's projected caption is in its target's display format already: DG-000123
  private static _configureBatches(table: DomainTable<CompoundBatchRow>): void {
    const caption = table.renderer.caption;
    table.renderer = {...table.renderer, listItem: (row) => MolTrackTables.line(caption(row),
      row[Rows.caption('compound_id')], row.epa_batch_id, row.project)};
  }

  private static _configureAdditions(table: DomainTable<AdditionRow>): void {
    const caption = table.renderer.caption;
    const role = (row: RowView<AdditionRow>) => isEmpty(row.role) ? '' :
      `${row.role.charAt(0)}${row.role.slice(1).toLowerCase()}`;
    const weight = (row: RowView<AdditionRow>) => isEmpty(row.molecular_weight) || row.molecular_weight! <= 0 ? '' :
      `MW ${DomainTable.captionText(table.properties, 'molecular_weight', row.molecular_weight)}`;
    table.renderer = {...table.renderer,
      listItem: (row) => MolTrackTables.line(caption(row), role(row), row.formula, weight(row))};
  }

  private static _configureAssayResults(table: DomainTable<AssayResultRow>): void {
    const caption = (row: RowView<AssayResultRow>) => Rows.isDraft(row) ? 'New assay result' :
      `${text(row[Rows.caption('batch_id')])} · ${text(row[Rows.caption('assay_run_id')])}`;
    const clint = (row: RowView<AssayResultRow>) => {
      if (isEmpty(row.reported_clint))
        return '';
      const qualifier = row.reported_clint_qualifier === '=' ? '' : row.reported_clint_qualifier ?? '';
      return `CLint ${qualifier}${row.reported_clint}`;
    };
    table.renderer = {...table.renderer, caption, listItem: (row) => MolTrackTables.line(caption(row),
      isEmpty(row.dosed_concentration) ? '' : `${row.dosed_concentration} uM`, clint(row))};
  }
}
