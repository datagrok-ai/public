import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {moltrackEmsDb} from '../generated/db';
import type {AssayResultInsert, AssayRunInsert, CompoundBatchInsert, CompoundInsert} from '../generated/db';
import {Chemistry, IDENTITY_COLUMNS} from './chemistry';
import {BulkRegistration} from './registration';
import {MolTrackTables} from './tables';

const ASSAY = 'Hepatocyte Stability';
const ALL = 10000;

type Values = Record<string, any>;

// additions and the assay are upserted; every other part loads only into an empty table
export class DemoData {
  static async load(): Promise<string> {
    const additions = DemoData.records(await DemoData.table('additions'));
    await moltrackEmsDb.additions.batch(additions.map((r) => ({...r,
      molecular_weight: Number(r.molecular_weight) > 0 ? Number(r.molecular_weight) : null})), {mode: 'upsert'});
    if (await moltrackEmsDb.compounds.count() === 0)
      await DemoData._compoundsAndBatches();
    await DemoData._assays();
    const [additionCount, compoundCount, batchCount, runCount, resultCount] = await Promise.all([
      moltrackEmsDb.additions.count(), moltrackEmsDb.compounds.count(), moltrackEmsDb.compoundBatches.count(),
      moltrackEmsDb.assayRuns.count(), moltrackEmsDb.assayResults.count()]);
    const summary = `MolTrack EMS demo: ${additionCount} additions, ${compoundCount} compounds, ` +
      `${batchCount} batches, ${runCount} assay runs, ${resultCount} assay results`;
    grok.log.info(summary);
    return summary;
  }

  static async table(name: string): Promise<DG.DataFrame> {
    return DG.DataFrame.fromCsv(await Chemistry.package().files.readAsText(`ems-demo/${name}.csv`));
  }

  static records(df: DG.DataFrame): Values[] {
    return Array.from({length: df.rowCount}, (_, i) =>
      Object.fromEntries(df.columns.toList().map((c) => [c.name, c.isNone(i) ? null : c.get(i)])));
  }

  private static async _compoundsAndBatches(): Promise<void> {
    const df = await DemoData.table('batches');
    const compounds = await MolTrackTables.compounds();
    const {ready, hashes} = await BulkRegistration.prepare(df, compounds);
    const label = (column: string) => BulkRegistration.label(compounds, column);
    const identity = (r: Values) =>
      Object.fromEntries(IDENTITY_COLUMNS.map((c) => [c, r[label(c)]])) as Omit<CompoundInsert, 'original_structure'>;
    await moltrackEmsDb.compounds.batch(DemoData.records(ready).map((r): CompoundInsert => ({...identity(r),
      original_structure: r[label('original_structure')], common_name: r['Common Name']?.trim(), cas: r['CAS'],
      epa_compound_id: r['EPA Compound ID'], source: r['Source']})));
    const compoundOf = new Map((await moltrackEmsDb.compounds.query({columns: ['hash_mol'], limit: ALL}))
      .map((c) => [c.hash_mol, c.id]));
    await moltrackEmsDb.compoundBatches.batch(DemoData.records(df).flatMap((r, i): CompoundBatchInsert[] =>
      compoundOf.has(hashes[i]) ? [{compound_id: compoundOf.get(hashes[i])!, epa_batch_id: r['EPA Batch ID'],
        source: r['Source'], source_batch_code: r['Source Batch Code'], responsible_party: r['Responsible Party'],
        eln_reference: r['ELN Reference'], purity: DemoData.number(r['Purity']), project: r['Project'],
        acquired_date: DemoData.date(r['Acquired Date']),
        synthesized_date: DemoData.date(r['Synthesized Date'])}] : []));
  }

  private static async _assays(): Promise<void> {
    await moltrackEmsDb.assays.upsert({name: ASSAY, assay_format: 'in cellulo', biological_system: 'hepatocyte'});
    const assay = (await moltrackEmsDb.assays.getByKey({name: ASSAY}))!;
    if (await moltrackEmsDb.assayRuns.count() === 0) {
      const runs = DemoData.records(await DemoData.table('assay_runs'));
      await moltrackEmsDb.assayRuns.batch(runs.map((r): AssayRunInsert => ({
        name: `${ASSAY}, ${String(r['cell lot']).trim()}, ${DemoData.day(r['assay date'])}`, assay_id: assay.id,
        cell_species: String(r['cell species']).trim(), cell_lot: String(r['cell lot']).trim(),
        run_date: DemoData.date(r['assay date']), assayer: String(r['assayer']).trim()})));
    }
    if (await moltrackEmsDb.assayResults.count() > 0)
      return;
    const runOf = new Map((await moltrackEmsDb.assayRuns.query({columns: ['cell_lot', 'run_date'], limit: ALL}))
      .map((r) => [`${r.cell_lot}|${DemoData.day(r.run_date)}`, r.id]));
    const batchOf = new Map((await moltrackEmsDb.compoundBatches.query({columns: ['epa_batch_id'], limit: ALL}))
      .map((b) => [b.epa_batch_id, b.id]));
    await moltrackEmsDb.assayResults.batch(DemoData.records(await DemoData.table('assay_results'))
      .map((r): AssayResultInsert => {
        const clint = DemoData.qualified(r['Reported CLint (µL/min/ 106 cells)']);
        return {
          assay_run_id: runOf.get(`${String(r['Cell Lot']).trim()}|${DemoData.day(r['assay_date'])}`)!,
          batch_id: batchOf.get(r['epa_batch_id'])!,
          dosed_concentration: DemoData.number(r['Concentration (µM)']),
          reported_clint: clint.value, reported_clint_qualifier: clint.qualifier,
          mean_htc_recovery: DemoData.number(r['Mean HTC recovery (%; n=2)']),
          sd_htc_recovery: DemoData.number(r['SD HTC recovery (%; n=2)']),
        };
      }));
  }

  // dates are UTC midnights read as written, or a day shifts with the time zone
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

  static date(x: unknown): string | undefined {
    return x === null || x === undefined || String(x).trim() === '' ? undefined : `${DemoData.day(x)}T00:00:00Z`;
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
