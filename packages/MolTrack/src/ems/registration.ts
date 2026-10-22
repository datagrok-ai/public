import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {DomainErrors, DomainTable, SYSTEM_COLUMNS} from '@datagrok-libraries/ems/src/dg';
import {isEmpty} from '@datagrok-libraries/u2';
import type {IProperty} from '@datagrok-libraries/u2';
import type {CompoundRow, CompoundBatchRow} from '../generated/db';
import {moltrackEmsDb} from '../generated/db';
import {MolTrackTables, Registered, SCHEMA} from './tables';

export type RegistrationStatus = 'success' | 'failed' | 'not_processed';

export interface RegistrationMapping {
  structure: string;
  compound: Record<string, string>;
  batch?: Record<string, string>;
}

export interface RegistrationOptions {
  rejectAll?: boolean;
}

export interface RegistrationPreview {
  verdict: 'new' | 'existing' | 'error';
  message: string;
  /** The earlier row of the file whose new compound this row is too. */
  sameAs?: number;
}

interface RowOutcome {
  status: RegistrationStatus;
  message: string;
  compoundId?: string;
  batchId?: string;
}

type Values = Record<string, unknown>;

interface Plan {
  rows: Values[];
  outcomes: RowOutcome[];
  /** each distinct structure's hash → the row whose values the compound gets */
  firstOf: Map<string, number>;
  registered: Map<string, Registered>;
}

const UNMAPPED = ['original_structure', 'molregno', 'batch_regno', 'compound_id'];

export class Registration {
  static readonly STATUS = 'Registration Status';
  static readonly MESSAGE = 'Registration Message';
  static readonly COMPOUND_ID = 'Grok ID';
  static readonly BATCH_ID = 'Grok Batch ID';
  static readonly OUTCOME =
    [Registration.STATUS, Registration.MESSAGE, Registration.COMPOUND_ID, Registration.BATCH_ID];
  static readonly MAX_OPS = 1000;

  static targets(table: DomainTable): IProperty[] {
    const writable = table.access.columnPolicy().forDraft();
    return table.properties.filter((p) => writable.field(p.name!) === 'editable' &&
      !SYSTEM_COLUMNS.some(([name]) => name === p.name) && !table.derived.targets.includes(p.name!) &&
      !UNMAPPED.includes(p.name!));
  }

  // batch columns win a shared name: a batch file's "Source" is the batch's
  static autoMap(df: DG.DataFrame, compounds: DomainTable<CompoundRow>,
    batches?: DomainTable<CompoundBatchRow>): RegistrationMapping {
    const targets = (table: DomainTable) => new Map(Registration.targets(table)
      .flatMap((p) => [p.name!, p.friendlyName ?? p.name!].map((n) => [DomainTable.matchKey(n), p.name!] as const)));
    const compound = targets(compounds);
    const batch = batches === undefined ? undefined : targets(batches);
    const mapping: RegistrationMapping = {structure: '', compound: {}, ...(batch === undefined ? {} : {batch: {}})};
    for (const name of df.columns.names()) {
      const k = DomainTable.matchKey(name);
      if (Registration.OUTCOME.includes(name))
        continue;
      if (mapping.structure === '' && ['smiles', 'structure', 'molecule', 'molfile', 'molblock'].includes(k))
        mapping.structure = name;
      else if (batch?.has(k))
        mapping.batch![name] = batch.get(k)!;
      else if (compound.has(k))
        mapping.compound[name] = compound.get(k)!;
    }
    return mapping;
  }

  static async register(df: DG.DataFrame, mapping: RegistrationMapping,
    options: RegistrationOptions = {}): Promise<DG.DataFrame> {
    const compounds = await MolTrackTables.compounds();
    const plan = await Registration._plan(df, mapping, compounds);
    if (options.rejectAll) {
      if (!plan.outcomes.some((o) => o.status === 'failed'))
        await Registration._atomic(df, mapping, plan, compounds);
      const failed = plan.outcomes.findIndex((o) => o.status === 'failed');
      if (failed >= 0) {
        for (const outcome of plan.outcomes.filter((o) => o.status === 'not_processed' && o.message === ''))
          outcome.message = `Not registered: row ${failed + 1} has a problem`;
      }
    } else
      await Registration._rowByRow(df, mapping, plan, compounds);
    return Registration._report(df, plan.outcomes, compounds);
  }

  /** What {@link register} would do with each row, without writing. */
  static async preview(df: DG.DataFrame, mapping: RegistrationMapping): Promise<RegistrationPreview[]> {
    const compounds = await MolTrackTables.compounds();
    const {rows, outcomes, firstOf, registered} = await Registration._plan(df, mapping, compounds);
    const batch = mapping.batch !== undefined;
    return outcomes.map((outcome, i): RegistrationPreview => {
      if (outcome.status === 'failed')
        return {verdict: 'error', message: outcome.message};
      const hash = rows[i].hash_mol as string;
      const existing = registered.get(hash);
      if (existing !== undefined) {
        const id = MolTrackTables.corporateId(compounds, 'molregno', existing.molregno);
        return {verdict: 'existing', message: batch ? `Registered as ${id}, adds a batch` : `Registered as ${id}`};
      }
      const first = firstOf.get(hash)!;
      if (first !== i)
        return {verdict: 'new', message: `Same compound as row ${first + 1}, adds a batch`, sameAs: first};
      return {verdict: 'new', message: batch ? 'New compound and batch' : 'New compound'};
    });
  }

  private static async _plan(df: DG.DataFrame, mapping: RegistrationMapping,
    compounds: DomainTable<CompoundRow>): Promise<Plan> {
    const structure = df.col(mapping.structure);
    const outcomes: RowOutcome[] = Array.from({length: df.rowCount}, () => ({status: 'not_processed', message: ''}));
    const rows = Array.from({length: df.rowCount}, (_, i) => Registration._values(df, i, mapping.compound,
      {original_structure: structure === null ? null : DomainTable.cellValue(structure, i)}));
    for (const [i, row] of rows.entries()) {
      if (isEmpty(row.original_structure))
        Object.assign(outcomes[i], {status: 'failed', message: 'Structure is required'});
    }
    for (const failure of await compounds.deriveRows(rows))
      Object.assign(outcomes[failure.index], {status: 'failed', message: failure.problem});
    const firstOf = new Map<string, number>();
    for (const [i, row] of rows.entries()) {
      if (outcomes[i].status !== 'failed' && !firstOf.has(row.hash_mol as string))
        firstOf.set(row.hash_mol as string, i);
    }
    // the derivation has just looked every one of these hashes up
    const registered = MolTrackTables.lookedUp(firstOf.keys());
    if (mapping.batch === undefined) {
      for (const [i, row] of rows.entries()) {
        if (outcomes[i].status === 'failed')
          continue;
        const hash = row.hash_mol as string;
        const existing = registered.get(hash);
        if (existing !== undefined) {
          Object.assign(outcomes[i], {status: 'failed', message: MolTrackTables.alreadyRegistered(
            MolTrackTables.corporateId(compounds, 'molregno', existing.molregno))});
        } else if (firstOf.get(hash) !== i)
          Object.assign(outcomes[i], {status: 'failed', message: `Same structure as row ${firstOf.get(hash)! + 1}`});
      }
    }
    return {rows, outcomes, firstOf, registered};
  }

  private static async _rowByRow(df: DG.DataFrame, mapping: RegistrationMapping, plan: Plan,
    compounds: DomainTable<CompoundRow>): Promise<void> {
    const {rows, outcomes, firstOf, registered} = plan;
    const idOf = new Map<string, string>([...registered].map(([hash, r]) => [hash, r.id]));
    const problemOf = new Map<string, string>();
    const written = [...firstOf.entries()]
      .filter(([hash, i]) => !registered.has(hash) && outcomes[i].status !== 'failed')
      .map(([, i]) => i);
    const report = await moltrackEmsDb.compounds.batch(written.map((i) => rows[i]), {mode: 'insert'});
    for (const r of report.rows) {
      const hash = rows[written[r.index]].hash_mol as string;
      if (r.status === 'inserted')
        idOf.set(hash, r.id!);
      else
        problemOf.set(hash, Registration._problem(r, compounds));
    }
    for (const [i, row] of rows.entries()) {
      if (outcomes[i].status === 'failed')
        continue;
      const id = idOf.get(row.hash_mol as string);
      Object.assign(outcomes[i], id === undefined ? {status: 'failed', message: problemOf.get(row.hash_mol as string)} :
        {status: 'success', compoundId: id});
    }
    if (mapping.batch === undefined)
      return;
    const batched = outcomes.map((o, i) => o.status === 'success' ? i : -1).filter((i) => i >= 0);
    const batchReport = await moltrackEmsDb.compoundBatches.batch(
      batched.map((i) => Registration._values(df, i, mapping.batch!, {compound_id: outcomes[i].compoundId})),
      {mode: 'insert'});
    const batches = await MolTrackTables.batches();
    for (const r of batchReport.rows) {
      const outcome = outcomes[batched[r.index]];
      if (r.status === 'inserted')
        outcome.batchId = r.id ?? undefined;
      else
        Object.assign(outcome, {status: 'failed', message: Registration._problem(r, batches)});
    }
  }

  // one transaction: a batch refers to a compound inserted beside it by its `$ref`
  private static async _atomic(df: DG.DataFrame, mapping: RegistrationMapping, plan: Plan,
    compounds: DomainTable<CompoundRow>): Promise<void> {
    const {rows, outcomes, firstOf, registered} = plan;
    const ops: DG.DomainTransactionOp[] = [];
    const rowOf: number[] = [];
    const refOf = new Map<string, string>();
    for (const [hash, i] of firstOf) {
      if (registered.has(hash))
        continue;
      refOf.set(hash, `c${refOf.size}`);
      ops.push({op: 'insert', table: 'compound', ref: refOf.get(hash), values: rows[i]});
      rowOf.push(i);
    }
    if (mapping.batch !== undefined) {
      for (const [i, row] of rows.entries()) {
        const hash = row.hash_mol as string;
        const compound = registered.get(hash)?.id ?? `$${refOf.get(hash)}`;
        ops.push({op: 'insert', table: 'compound_batch',
          values: Registration._values(df, i, mapping.batch, {compound_id: compound})});
        rowOf.push(i);
      }
    }
    if (ops.length > Registration.MAX_OPS) {
      for (const outcome of outcomes) {
        outcome.message = `Reject all writes at most ${Registration.MAX_OPS} rows at once (this file needs ` +
          `${ops.length}): split the file`;
      }
      return;
    }
    try {
      const results = await grok.dapi.domains.transaction(SCHEMA, ops);
      const idOf = new Map<string, string>([...registered].map(([hash, r]) => [hash, r.id]));
      for (const [k, result] of results.entries()) {
        if (ops[k].table === 'compound')
          idOf.set(rows[rowOf[k]].hash_mol as string, result.id);
        else
          outcomes[rowOf[k]].batchId = result.id;
      }
      for (const [i, row] of rows.entries())
        Object.assign(outcomes[i], {status: 'success', compoundId: idOf.get(row.hash_mol as string)});
    } catch (e) {
      const at = (e as DG.DomainError).opIndex;
      const row = at === undefined ? undefined : rowOf[at];
      const table = at !== undefined && ops[at].table === 'compound_batch' ? await MolTrackTables.batches() : compounds;
      const errors = (e as DG.DomainValidationError).rows?.[0]?.errors;
      const message = errors?.length ? Registration._problem({status: 'error', errors}, table) :
        DomainErrors.message(e);
      if (row === undefined) {
        for (const outcome of outcomes)
          outcome.message = message;
      } else
        Object.assign(outcomes[row], {status: 'failed', message});
    }
  }

  private static _values(df: DG.DataFrame, i: number, mapping: Record<string, string>, extra: Values): Values {
    const values: Values = {...extra};
    for (const [source, target] of Object.entries(mapping)) {
      const raw = DomainTable.cellValue(df.col(source)!, i);
      values[target] = typeof raw === 'string' ? (raw.trim() === '' ? null : raw.trim()) : raw;
    }
    return values;
  }

  private static _problem(row: {status: string, errors?: {column?: string, message: string}[]},
    table: DomainTable): string {
    if (row.errors?.length) {
      return row.errors.map((e) => {
        const prop = table.properties.find((p) => p.name === e.column);
        return prop === undefined ? e.message : `${prop.friendlyName ?? prop.name}: ${e.message}`;
      }).join('; ');
    }
    return row.status === 'duplicate' ? MolTrackTables.alreadyRegistered() : row.status;
  }

  private static async _report(df: DG.DataFrame, outcomes: RowOutcome[],
    compounds: DomainTable<CompoundRow>): Promise<DG.DataFrame> {
    const report = df.clone();
    for (const name of Registration.OUTCOME) {
      if (report.col(name) !== null)
        report.columns.remove(name);
    }
    const numbers = async (table: 'compounds' | 'compoundBatches', column: 'molregno' | 'batch_regno',
      of: (o: RowOutcome) => string | undefined) => {
      const wanted = [...new Set(outcomes.map(of).filter((x): x is string => x !== undefined))];
      const rows: Values[] = await MolTrackTables.inChunks(wanted, (chunk) =>
        (moltrackEmsDb[table] as DG.DomainTableClient).query(
          {filter: MolTrackTables.within('id', chunk), columns: [column], limit: chunk.length}));
      return new Map(rows.map((row) => [row.id as string, row[column]]));
    };
    const molregnos = await numbers('compounds', 'molregno', (o) => o.compoundId);
    const batches = await MolTrackTables.batches();
    const regnos = await numbers('compoundBatches', 'batch_regno', (o) => o.batchId);
    report.columns.addNewString(Registration.STATUS).init((i) => outcomes[i].status);
    report.columns.addNewString(Registration.MESSAGE).init((i) => outcomes[i].message);
    report.columns.addNewString(Registration.COMPOUND_ID).init((i) => outcomes[i].compoundId === undefined ? '' :
      MolTrackTables.corporateId(compounds, 'molregno', molregnos.get(outcomes[i].compoundId!)));
    if (outcomes.some((o) => o.batchId !== undefined)) {
      report.columns.addNewString(Registration.BATCH_ID).init((i) => outcomes[i].batchId === undefined ? '' :
        MolTrackTables.corporateId(batches, 'batch_regno', regnos.get(outcomes[i].batchId!)));
    }
    const outcome = Registration.OUTCOME.filter((name) => report.col(name) !== null);
    report.columns.setOrder([...outcome, ...report.columns.names().filter((name) => !outcome.includes(name))]);
    return report;
  }
}
