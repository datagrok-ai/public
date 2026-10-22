import type * as DG from 'datagrok-api/dg';
import {DomainTable, SYSTEM_COLUMNS} from '@datagrok-libraries/ems/src/dg';
import {isEmpty} from '@datagrok-libraries/u2';
import type {IProperty} from '@datagrok-libraries/u2';
import type {CompoundRow} from '../generated/db';
import {moltrackEmsDb} from '../generated/db';
import {MolTrackTables} from './tables';
import type {Registered} from './tables';

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
  compound?: Registered;
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
    const rejectAll = options.rejectAll === true;
    // the first row of each new structure registers it
    const written = [...plan.firstOf.values()].filter((i) => plan.outcomes[i].status === 'not_processed');
    if (written.length > 0 && !(rejectAll && plan.outcomes.some((o) => o.status === 'failed')))
      await Registration._write(plan, written, compounds, rejectAll);
    const failed = plan.outcomes.findIndex((o) => o.status === 'failed');
    if (rejectAll && failed >= 0) {
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
        Object.assign(outcomes[i], {status: 'existing', compound: existing, message:
          MolTrackTables.alreadyRegistered(MolTrackTables.corporateId(compounds, 'molregno', existing.molregno))});
      }
    }
    return {rows, outcomes, firstOf};
  }

  /** Every row still pending, by the compound its structure got: the first row registered it, a later one repeats
   * it. */
  private static _settle({rows, outcomes, firstOf}: Plan, landed: Map<string, Registered>,
    problemOf: Map<string, string>): void {
    for (const [i, row] of rows.entries()) {
      if (outcomes[i].status !== 'not_processed')
        continue;
      const hash = row.hash_mol as string;
      const compound = landed.get(hash);
      const first = firstOf.get(hash)!;
      Object.assign(outcomes[i], compound === undefined ? {status: 'failed', message: problemOf.get(hash)} :
        first === i ? {status: 'success', compound} :
          {status: 'existing', compound, message: `Same compound as row ${first + 1}`});
    }
  }

  /** One batch of the rows `written`: with `allOrNothing` a row the server refuses leaves every row unwritten. */
  private static async _write(plan: Plan, written: number[], compounds: DomainTable<CompoundRow>,
    allOrNothing: boolean): Promise<void> {
    const hashOf = (k: number) => plan.rows[written[k]].hash_mol as string;
    const report = await moltrackEmsDb.compounds.batch(written.map((i) => plan.rows[i]),
      {mode: 'insert', allOrNothing});
    const problemOf = new Map(report.rows.filter((r) => r.status !== 'inserted')
      .map((r): [string, string] => [hashOf(r.index), Registration._problem(r, compounds)]));
    if (report.error !== undefined) {
      for (const r of report.rows)
        Object.assign(plan.outcomes[written[r.index]], {status: 'failed', message: problemOf.get(hashOf(r.index))});
      return;
    }
    // the report names at most 1000 rows: the written ones are read back by their hash
    Registration._settle(plan, await MolTrackTables.registered('hash_mol',
      written.map((_, k) => hashOf(k)).filter((hash) => !problemOf.has(hash))), problemOf);
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

  private static _report(df: DG.DataFrame, outcomes: RowOutcome[], compounds: DomainTable<CompoundRow>):
    DG.DataFrame {
    const report = df.clone();
    for (const name of Registration.OUTCOME) {
      if (report.col(name) !== null)
        report.columns.remove(name);
    }
    report.columns.addNewString(Registration.STATUS).init((i) => outcomes[i].status);
    report.columns.addNewString(Registration.MESSAGE).init((i) => outcomes[i].message);
    report.columns.addNewString(Registration.COMPOUND_ID).init((i) =>
      MolTrackTables.corporateId(compounds, 'molregno', outcomes[i].compound?.molregno));
    report.columns.setOrder([...Registration.OUTCOME,
      ...report.columns.names().filter((name) => !Registration.OUTCOME.includes(name))]);
    return report;
  }
}
