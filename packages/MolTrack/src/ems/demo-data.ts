import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {DomainPick, DomainTable} from '@datagrok-libraries/ems/src/dg';
import {moltrackEmsDb} from '../generated/db';
import type {AssayResultInsert, AssayRunInsert} from '../generated/db';
import {Chemistry} from './chemistry';
import {Registration} from './registration';
import {MolTrackTables, SCHEMA} from './tables';

const ASSAY = 'Hepatocyte Stability';

export class DemoData {
  static async table(name: string): Promise<DG.DataFrame> {
    return DG.DataFrame.fromCsv(await Chemistry.package().files.readAsText(`ems-demo/${name}.csv`));
  }

  static records(df: DG.DataFrame): Record<string, any>[] {
    return Array.from({length: df.rowCount}, (_, i) =>
      Object.fromEntries(df.columns.toList().map((c) => [c.name, DomainTable.cellValue(c, i)])));
  }

  // additions and the assay are upserted; every other part is skipped where its table has rows already
  static async load(): Promise<string> {
    const additions = await DemoData.table('additions');
    await moltrackEmsDb.additions.batch(DemoData.records(additions).map((r) => ({...r,
      molecular_weight: Number(r.molecular_weight) > 0 ? Number(r.molecular_weight) : null})), {mode: 'upsert'});

    const batches = await DemoData.table('batches');
    let failed = 0;
    if (await moltrackEmsDb.compounds.count() === 0) {
      // a column the batch has too (Source) is the batch's
      const ofBatch = Registration.autoMap(batches, await MolTrackTables.batches()).columns;
      const mapping = Registration.autoMap(batches, await MolTrackTables.compounds());
      for (const name of Object.keys(ofBatch))
        delete mapping.columns[name];
      const registered = await Registration.register(batches, mapping);
      const grokIds = registered.col(Registration.COMPOUND_ID)!.toList() as (string | null)[];
      const idOf = await DomainPick.lookUp(`${SCHEMA}.compound`, grokIds.filter((id): id is string => !!id));
      const rows = grokIds.flatMap((id, i) => !id ? [] :
        [Registration.values(batches, i, ofBatch, {compound_id: idOf.get(id)})]);
      await moltrackEmsDb.compoundBatches.batch(rows, {mode: 'insert'});
      failed = batches.rowCount - rows.length;
    }

    await moltrackEmsDb.assays.upsert({name: ASSAY, assay_format: 'in cellulo', biological_system: 'hepatocyte'});
    const assay = (await moltrackEmsDb.assays.getByKey({name: ASSAY}))!;
    const runs = DemoData.records(await DemoData.table('assay_runs')).map((r): AssayRunInsert => ({
      name: `${ASSAY}, ${String(r['cell lot']).trim()}, ${DemoData.day(r['assay date'])}`, assay_id: assay.id,
      cell_species: String(r['cell species']).trim(), cell_lot: String(r['cell lot']).trim(),
      run_date: `${DemoData.day(r['assay date'])}T00:00:00Z`, assayer: String(r['assayer']).trim()}));
    if (await moltrackEmsDb.assayRuns.count() === 0)
      await moltrackEmsDb.assayRuns.batch(runs);
    const runOf = new Map((await MolTrackTables.allPages((page) =>
      moltrackEmsDb.assayRuns.query({columns: ['cell_lot', 'run_date'], sort: 'id', ...page})))
      .map((r) => [`${r.cell_lot}|${DemoData.day(r.run_date)}`, r.id]));

    const batchOf = new Map((await MolTrackTables.allPages((page) =>
      moltrackEmsDb.compoundBatches.query({columns: ['epa_batch_id'], sort: 'batch_regno', ...page})))
      .map((b) => [b.epa_batch_id, b.id]));
    const results = DemoData.records(await DemoData.table('assay_results')).map((r): AssayResultInsert => {
      const date = DemoData.day(r['assay_date']);
      const clint = DemoData.qualified(r['Reported CLint (µL/min/ 106 cells)']);
      return {
        assay_run_id: runOf.get(`${String(r['Cell Lot']).trim()}|${date}`)!, batch_id: batchOf.get(r['epa_batch_id'])!,
        dosed_concentration: DemoData.number(r['Concentration (µM)']),
        reported_clint: clint.value, reported_clint_qualifier: clint.qualifier,
        mean_htc_recovery: DemoData.number(r['Mean HTC recovery (%; n=2)']),
        sd_htc_recovery: DemoData.number(r['SD HTC recovery (%; n=2)']),
      };
    });
    if (await moltrackEmsDb.assayResults.count() === 0)
      await moltrackEmsDb.assayResults.batch(results);
    const [additionCount, compoundCount, batchCount, runCount, resultCount] = await Promise.all([
      moltrackEmsDb.additions.count(), moltrackEmsDb.compounds.count(), moltrackEmsDb.compoundBatches.count(),
      moltrackEmsDb.assayRuns.count(), moltrackEmsDb.assayResults.count()]);
    const summary = `MolTrack EMS demo: ${additionCount} additions, ${compoundCount} compounds, ` +
      `${batchCount} batches (${failed} rows failed this time), ${runCount} assay runs, ${resultCount} assay results`;
    grok.log.info(summary);
    return summary;
  }

  // dates are UTC midnights: text is read as written and a stored date in UTC, else either shifts a day
  static day(x: unknown): string {
    const text = typeof x === 'string' ? x.trim() : '';
    const iso = /^\d{4}-\d{2}-\d{2}/.exec(text);
    if (iso !== null)
      return iso[0];
    const us = /^(\d{1,2})\/(\d{1,2})\/(\d{4})$/.exec(text);
    if (us !== null)
      return `${us[3]}-${us[1].padStart(2, '0')}-${us[2].padStart(2, '0')}`;
    return dayjs(x as dayjs.ConfigType).toISOString().slice(0, 10);
  }

  static number(x: unknown): number | undefined {
    const n = typeof x === 'number' ? x : Number(String(x ?? '').trim());
    return x === null || x === undefined || String(x).trim() === '' || Number.isNaN(n) ? undefined : n;
  }

  static qualified(x: unknown): {value?: number, qualifier?: '=' | '<' | '>'} {
    const text = String(x ?? '').trim();
    const qualifier = text.startsWith('<') ? '<' : text.startsWith('>') ? '>' : '=';
    const value = DemoData.number(qualifier === '=' ? text : text.slice(1));
    return value === undefined ? {} : {value, qualifier};
  }
}
