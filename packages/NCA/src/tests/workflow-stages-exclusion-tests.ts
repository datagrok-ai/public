/**
 * nca-workflow-stages S5 — Exclude λz ("not reportable") reaches the statistics exactly as PKNCA's
 * `exclude()` does (rule 18, DATA-GAP-47; U1 / AD-6).
 *
 * The oracle is the REGENERATED `theoph-exclusion-reference.json` (PKNCA 0.12.1), read at test
 * time — never a transcribed literal. Tolerances are the NFR-04/05 bands: AUC / Cmax 0.1 %,
 * λz 0.5 %, derived 1 %; N is exact.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';

import {DEFAULT_RULES} from '../adapter/defaults';
import {reportableResults} from '../adapter/reportable';
import {aggregateNcaResults} from '../adapter/results-aggregator';
import {loadFixture, runPipelineOrThrow, withDoseMg, theophContract, theophScalars} from './test-utils';

const TIMEOUT = 60000;

interface ExclusionOracle {
  pknca_version: string;
  excluded: {subject: string; reason: string}[];
  summary: Record<string, {mean: number | null; sd: number | null; n: number; n_excluded: number}>;
}

/** Relative tolerance per parameter (NFR-04 / NFR-05). */
const REL_TOL: Readonly<Record<string, number>> = {
  CMAX: 0.001, AUC_LAST: 0.001, AUC_INF: 0.001, AUMC_INF: 0.001,
  LAMBDA_Z: 0.005, HALF_LIFE: 0.01, CL: 0.01, VZ: 0.01, MRT: 0.01, PCT_EXTRAP: 0.01, TMAX: 0.001,
};

category('Workflow stages: exclusion oracle (S5.3, rule 18)', () => {
  test('AC-S5.3.1 aggregate(reportable(theoph, 2 excluded)) = PKNCA exclude() MEAN/SD/N', async () => {
    const oracle = JSON.parse(await grok.dapi.files.readAsText(
      'System:AppData/Nca/theoph-exclusion-reference.json')) as ExclusionOracle;
    const theoph = await loadFixture('theoph.csv');
    const {results} = runPipelineOrThrow(
      withDoseMg(theoph), theophContract(theoph), theophScalars(), DEFAULT_RULES);
    const excluded = new Set(oracle.excluded.map((e) => e.subject));
    expect(excluded.size, 2, 'the oracle excludes two subjects');
    const rep = reportableResults(results, excluded);
    expect(rep.excludedCount, 2, 'both subjects found in the results frame');

    const grouped = rep.df.clone();
    grouped.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ALL', new Array(grouped.rowCount).fill('all')));
    const agg = aggregateNcaResults(grouped, 'ALL');

    // The oracle's CL / VZ / VSS are PKNCA's ABSOLUTE values (theoph doses mg) — the app reports them
    // on the absolute companions since nca-cl-bw-normalised-units (CL / VZ / VSS are per kg).
    const COMPANION: Readonly<Record<string, string>> = {CL: 'CL_ABS', VZ: 'VZ_ABS', VSS: 'VSS_ABS'};
    let checked = 0;
    for (const [oracleParam, want] of Object.entries(oracle.summary)) {
      const param = COMPANION[oracleParam] ?? oracleParam;
      const nCol = agg.col(`N_${param}`);
      if (nCol === null) continue; // not an aggregated parameter here
      expect(Number(nCol.get(0)), want.n, `N_${param} (PKNCA n, ${want.n_excluded} excluded)`);
      const tol = REL_TOL[param] ?? 0.01;
      if (want.mean !== null) {
        const mean = Number(agg.col(`MEAN_${param}`)!.get(0));
        expectFloat(mean, want.mean, Math.abs(want.mean) * tol, `MEAN_${param}`);
      }
      if (want.sd !== null && agg.col(`SD_${param}`) !== null) {
        const sd = Number(agg.col(`SD_${param}`)!.get(0));
        expectFloat(sd, want.sd, Math.abs(want.sd) * tol, `SD_${param}`); // same NFR band as the mean (rule 18)
      }
      checked++;
    }
    expect(checked >= 10, true, `enough parameters checked (${checked})`);
  }, {timeout: TIMEOUT});

  test('AC-S5.1.2 excluding one subject: λz-dependent N −1, N_CMAX unchanged', async () => {
    const theoph = await loadFixture('theoph.csv');
    const {results} = runPipelineOrThrow(
      withDoseMg(theoph), theophContract(theoph), theophScalars(), DEFAULT_RULES);
    const addAll = (df: DG.DataFrame): DG.DataFrame => {
      const c = df.clone();
      c.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ALL', new Array(c.rowCount).fill('all')));
      return c;
    };
    const raw = aggregateNcaResults(addAll(results), 'ALL');
    const ex = aggregateNcaResults(addAll(reportableResults(results, new Set(['3'])).df), 'ALL');
    const n = (a: DG.DataFrame, p: string): number => Number(a.col(`N_${p}`)!.get(0));
    for (const p of ['HALF_LIFE', 'AUC_INF', 'CL', 'VZ', 'LAMBDA_Z', 'MRT'])
      expect(n(ex, p), n(raw, p) - 1, `N_${p} drops by one`);
    for (const p of ['CMAX', 'TMAX', 'AUC_LAST'])
      expect(n(ex, p), n(raw, p), `N_${p} unchanged`);
  }, {timeout: TIMEOUT});
});
