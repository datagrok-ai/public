import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {DomainErrors, DomainTable, SYSTEM_COLUMNS} from '@datagrok-libraries/ems/src/dg';
import {isEmpty} from '@datagrok-libraries/u2';
import type {IProperty} from '@datagrok-libraries/u2';
import type {CompoundRow} from '../generated/db';
import {moltrackEmsDb} from '../generated/db';
import {MolTrackTables, SCHEMA} from './tables';

/** `existing`: the row's compound was registered already, before the run or by an earlier row of the file. */
export type RegistrationStatus = 'success' | 'existing' | 'failed' | 'not_processed';

/** The file's structure column, and the file's other columns by the table column each fills. */
export interface RegistrationMapping {
  structure: string;
  columns: Record<string, string>;
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
}

type Values = Record<string, unknown>;

interface Plan {
  rows: Values[];
  outcomes: RowOutcome[];
  /** each distinct structure's hash → the row whose values the compound gets */
  firstOf: Map<string, number>;
}

export class Registration {
  static readonly STATUS = 'Registration Status';
  static readonly MESSAGE = 'Registration Message';
  static readonly COMPOUND_ID = 'Grok ID';
  static readonly OUTCOME = [Registration.STATUS, Registration.MESSAGE, Registration.COMPOUND_ID];
  static readonly MAX_OPS = 1000;

  /** The columns of `table` a file may fill. */
  static targets(table: DomainTable): IProperty[] {
    const writable = table.access.columnPolicy().forDraft();
    return table.properties.filter((p) => writable.field(p.name!) === 'editable' &&
      !SYSTEM_COLUMNS.some(([name]) => name === p.name) && !table.derived.targets.includes(p.name!) &&
      p.name !== 'original_structure');
  }

  /** The file's columns matched to the columns of `table` by name or caption. */
  static autoMap(df: DG.DataFrame, table: DomainTable): RegistrationMapping {
    const targets = new Map(Registration.targets(table)
      .flatMap((p) => [p.name!, p.friendlyName ?? p.name!].map((n) => [DomainTable.matchKey(n), p.name!] as const)));
    const mapping: RegistrationMapping = {structure: '', columns: {}};
    for (const name of df.columns.names()) {
      const k = DomainTable.matchKey(name);
      if (Registration.OUTCOME.includes(name))
        continue;
      if (mapping.structure === '' && ['smiles', 'structure', 'molecule', 'molfile', 'molblock'].includes(k))
        mapping.structure = name;
      else if (targets.has(k))
        mapping.columns[name] = targets.get(k)!;
    }
    return mapping;
  }

  /** Registers each new structure of the file once; a row whose compound is registered already is reported with
   * its Grok ID and not written. */
  static async register(df: DG.DataFrame, mapping: RegistrationMapping,
    options: RegistrationOptions = {}): Promise<DG.DataFrame> {
    const compounds = await MolTrackTables.compounds();
    const plan = await Registration._plan(df, mapping, compounds);
    // the first row of each new structure registers it
    const written = [...plan.firstOf.values()].filter((i) => plan.outcomes[i].status === 'not_processed');
    if (written.length > 0 && !options.rejectAll)
      await Registration._rowByRow(plan, written, compounds);
    else if (written.length > 0 && !plan.outcomes.some((o) => o.status === 'failed'))
      await Registration._atomic(plan, written, compounds);
    const failed = plan.outcomes.findIndex((o) => o.status === 'failed');
    if (options.rejectAll && failed >= 0) {
      for (const outcome of plan.outcomes.filter((o) => o.status === 'not_processed' && o.message === ''))
        outcome.message = `Not registered: row ${failed + 1} has a problem`;
    }
    return Registration._report(df, plan.outcomes, compounds);
  }

  /** What {@link register} would do with each row, without writing. */
  static async preview(df: DG.DataFrame, mapping: RegistrationMapping): Promise<RegistrationPreview[]> {
    const compounds = await MolTrackTables.compounds();
    const {rows, outcomes, firstOf} = await Registration._plan(df, mapping, compounds);
    return outcomes.map((outcome, i): RegistrationPreview => {
      if (outcome.status === 'failed')
        return {verdict: 'error', message: outcome.message};
      if (outcome.status === 'existing')
        return {verdict: 'existing', message: outcome.message};
      const first = firstOf.get(rows[i].hash_mol as string)!;
      return first === i ? {verdict: 'new', message: 'New compound'} :
        {verdict: 'new', message: `Same compound as row ${first + 1}`, sameAs: first};
    });
  }

  private static async _plan(df: DG.DataFrame, mapping: RegistrationMapping,
    compounds: DomainTable<CompoundRow>): Promise<Plan> {
    const structure = df.col(mapping.structure);
    const outcomes: RowOutcome[] = Array.from({length: df.rowCount}, () => ({status: 'not_processed', message: ''}));
    const rows = Array.from({length: df.rowCount}, (_, i) => Registration.values(df, i, mapping.columns,
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
    for (const [i, row] of rows.entries()) {
      const existing = registered.get(row.hash_mol as string);
      if (outcomes[i].status !== 'failed' && existing !== undefined) {
        Object.assign(outcomes[i], {status: 'existing', compoundId: existing.id, message:
          MolTrackTables.alreadyRegistered(MolTrackTables.corporateId(compounds, 'molregno', existing.molregno))});
      }
    }
    return {rows, outcomes, firstOf};
  }

  /** Every row still pending, by the compound its structure got: the first row registered it, a later one repeats
   * it. */
  private static _settle({rows, outcomes, firstOf}: Plan, idOf: Map<string, string>,
    problemOf: Map<string, string> = new Map()): void {
    for (const [i, row] of rows.entries()) {
      if (outcomes[i].status !== 'not_processed')
        continue;
      const hash = row.hash_mol as string;
      const id = idOf.get(hash);
      const first = firstOf.get(hash)!;
      Object.assign(outcomes[i], id === undefined ? {status: 'failed', message: problemOf.get(hash)} :
        first === i ? {status: 'success', compoundId: id} :
          {status: 'existing', compoundId: id, message: `Same compound as row ${first + 1}`});
    }
  }

  private static async _rowByRow(plan: Plan, written: number[], compounds: DomainTable<CompoundRow>):
    Promise<void> {
    const report = await moltrackEmsDb.compounds.batch(written.map((i) => plan.rows[i]), {mode: 'insert'});
    const idOf = new Map<string, string>();
    const problemOf = new Map<string, string>();
    for (const r of report.rows) {
      const hash = plan.rows[written[r.index]].hash_mol as string;
      if (r.status === 'inserted')
        idOf.set(hash, r.id!);
      else
        problemOf.set(hash, Registration._problem(r, compounds));
    }
    Registration._settle(plan, idOf, problemOf);
  }

  // one transaction: a refusal of any compound writes none
  private static async _atomic(plan: Plan, written: number[], compounds: DomainTable<CompoundRow>):
    Promise<void> {
    const pending = plan.outcomes.filter((o) => o.status === 'not_processed');
    if (written.length > Registration.MAX_OPS) {
      for (const outcome of pending) {
        outcome.message = `Reject all writes at most ${Registration.MAX_OPS} rows at once (this file needs ` +
          `${written.length}): split the file`;
      }
      return;
    }
    try {
      const results = await grok.dapi.domains.transaction(SCHEMA,
        written.map((i): DG.DomainTransactionOp => ({op: 'insert', table: 'compound', values: plan.rows[i]})));
      Registration._settle(plan, new Map(results.map((result, k) =>
        [plan.rows[written[k]].hash_mol as string, result.id])));
    } catch (e) {
      const at = (e as DG.DomainError).opIndex;
      const errors = (e as DG.DomainValidationError).rows?.[0]?.errors;
      const message = errors?.length ? Registration._problem({status: 'error', errors}, compounds) :
        DomainErrors.message(e);
      if (at === undefined) {
        for (const outcome of pending)
          outcome.message = message;
      } else
        Object.assign(plan.outcomes[written[at]], {status: 'failed', message});
    }
  }

  /** Row `i` of the file under the table's column names: text trimmed, an empty cell null. */
  static values(df: DG.DataFrame, i: number, mapping: Record<string, string>, extra: Values): Values {
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
    const wanted = [...new Set(outcomes.map((o) => o.compoundId).filter((x): x is string => x !== undefined))];
    const rows = await MolTrackTables.inChunks(wanted, (chunk) => moltrackEmsDb.compounds.query(
      {filter: MolTrackTables.within('id', chunk), columns: ['molregno'], limit: chunk.length}));
    const molregnos = new Map(rows.map((row) => [row.id, row.molregno]));
    report.columns.addNewString(Registration.STATUS).init((i) => outcomes[i].status);
    report.columns.addNewString(Registration.MESSAGE).init((i) => outcomes[i].message);
    report.columns.addNewString(Registration.COMPOUND_ID).init((i) => outcomes[i].compoundId === undefined ? '' :
      MolTrackTables.corporateId(compounds, 'molregno', molregnos.get(outcomes[i].compoundId!)));
    report.columns.setOrder([...Registration.OUTCOME,
      ...report.columns.names().filter((name) => !Registration.OUTCOME.includes(name))]);
    return report;
  }
}
