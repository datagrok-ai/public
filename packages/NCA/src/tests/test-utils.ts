/**
 * Test utilities — fixture loaders, column-contract presets, and pipeline
 * driver for the UC-01 Datagrok test suite (`grok test`).
 *
 * **File loading.** Fixtures ship under `files/` and are uploaded to
 * the Datagrok server on every `grok publish`. Inside the test runtime they
 * are reachable via `grok.dapi.files.readAsText('System:AppData/Nca/<name>')`
 * (same path the Phase 6 `ncaApp` launcher uses for the "Use sample data"
 * dropdown).
 *
 * **Column-contract presets.** `theoph` and `indometh` use different CSV
 * column names (`Subject/Wt/Dose/Time/conc` vs `Subject/time/conc`), so the
 * presets are CSV-specific. The Subject column's runtime type depends on
 * how `DG.DataFrame.fromCsv` autodetects quoted-integer ids — we probe the
 * loaded DataFrame and build the matching `rawType` so the validator passes
 * regardless of platform-side autodetection.
 *
 * **Reference-fixture shape.** Reference JSONs ship snake_case from PKNCA
 * (v0.12.1): `cmax / tmax / auclast / aucinf / pct_aucextrap / lambda_z /
 * half_life / cl / vz`, plus `lambda_z_n_points / lambda_z_r_squared /
 * lambda_z_adj_r_squared / lambda_z_time_first / lambda_z_time_last /
 * clast_obs / c0_extrapolated` under `provenance`. No `lambda_z_points_used`
 * array — only the count. The reference-suite test asserts per-subject
 * params and the points-count regression marker (spec §8.1 row 11 is
 * narrowed to count-only because the fixture only carries the count).
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {nca} from '@datagrok-libraries/sci-comp';

import {runNcaPipeline} from '../controller/nca-controller';
import type {NcaRunOutcome, RunNcaPipelineOptions} from '../controller/nca-controller';
import {normalizeNumericNaMarkers} from '../adapter/na-normalize';
import {reparseFloatColumnsToFloat64} from '../adapter/csv-float64';
import {scanNumericSentinels, describeSentinelEntries} from '../adapter/adapter-utils';
import {expect} from '@datagrok-libraries/test/src/test';
import {normalizeRouteColumn, DERIVED_ROUTE_CODE_COL} from '../adapter/route-normalize';
import type {ColumnContract} from '../adapter/column-contract';
import type {ProfileScalars} from '../adapter/profile-inputs';

const FIXTURE_PATH_PREFIX = 'System:AppData/Nca';

/** Subset of the PKNCA reference JSON the reference-suite test asserts on. */
export interface FixtureParameters {
  readonly cmax: number | null;
  readonly tmax: number | null;
  readonly auclast: number | null;
  readonly aucinf: number | null;
  readonly pct_aucextrap: number | null;
  readonly lambda_z: number | null;
  readonly half_life: number | null;
  readonly cl: number | null;
  readonly vz: number | null;
  // FR-200 moment / lag params (sci-comp 0.7.0). `null` = not applicable for
  // the route: `vss` null for extravascular, `tlag` null for IV.
  readonly aumclast?: number | null;
  readonly aumcinf_obs?: number | null;
  readonly mrt?: number | null;
  readonly vss?: number | null;
  readonly tlag?: number | null;
  readonly pct_aumcextrap?: number | null;
}

export interface FixtureProvenance {
  readonly lambda_z_n_points: number | null;
  readonly lambda_z_r_squared: number | null;
  readonly lambda_z_adj_r_squared: number | null;
  readonly lambda_z_time_first: number | null;
  readonly lambda_z_time_last: number | null;
  readonly clast_obs: number | null;
  readonly c0_extrapolated: number | null;
  /** PKNCA `span.ratio` — the λz window's span in half-lives. OPTIONAL because
   *  the oracle exists on only part of the corpus: the theoph / indometh JSON
   *  fixtures carry it (regenerated from sci-comp's own 0.10.0 fixtures, which
   *  request it from PKNCA in `regen-fixtures.R`) and the rat / multicompound
   *  reference CSVs carry a `span.ratio` column, but the NONMEM and IV-infusion
   *  JSONs predate it. Absent ⇒ `undefined`, which every assertion must treat as
   *  "no oracle here", never as "expected null". */
  readonly span_ratio?: number | null;
}

export interface FixtureProfile {
  readonly profile_key: {readonly subject: string; readonly route: string};
  readonly parameters: FixtureParameters;
  readonly provenance: FixtureProvenance;
}

export interface FixtureFile {
  readonly dataset: string;
  readonly pknca_version: string;
  readonly profiles: ReadonlyArray<FixtureProfile>;
}

/**
 * Load a CSV fixture bundled with the package and parse it via the platform
 * CSV reader. Returns the raw DataFrame as the platform parses it; the
 * caller layers a column-contract preset over it (see {@link theophContract},
 * {@link indomethContract}).
 */
export async function loadFixture(csvName: string): Promise<DG.DataFrame> {
  const csv = await grok.dapi.files.readAsText(`${FIXTURE_PATH_PREFIX}/${csvName}`);
  const df = DG.DataFrame.fromCsv(csv);
  // BUG-01: DG parses numeric CSV through Float32 (precision lost at parse).
  // Re-parse FLOAT columns from the source text into Float64 so Cmax/Tmax hit
  // PKNCA bit-identity (NFR-04 / §8.1 exact-match). Disjoint from #15 below.
  reparseFloatColumnsToFloat64(df, csv);
  // Mirror the production launch path: NONMEM-style `.`/empty NA markers are
  // normalized STRING→FLOAT so DV/AMT/RATE are numeric for the reference suite
  // (Ed #15, VAL-01 wiring of 03-06). No-op for theoph/indometh (no markers).
  normalizeNumericNaMarkers(df);
  return df;
}

/** Load and parse a reference JSON. */
export async function loadReference(jsonName: string): Promise<FixtureFile> {
  const text = await grok.dapi.files.readAsText(`${FIXTURE_PATH_PREFIX}/${jsonName}`);
  return JSON.parse(text) as FixtureFile;
}

/**
 * Build a quick lookup `Map<subject, FixtureProfile>` so reference-suite
 * assertions can pull the expected row for a given subject in O(1).
 */
export function indexFixtureProfiles(fx: FixtureFile): Map<string, FixtureProfile> {
  const out = new Map<string, FixtureProfile>();
  for (const p of fx.profiles) out.set(p.profile_key.subject, p);
  return out;
}

/**
 * Inspect a column and return the `ColumnRawType` the contract should
 * declare. The Datagrok `DataFrame.fromCsv` autodetector resolves
 * quoted-integer columns to `int` and quoted non-integer strings to `string`;
 * pinning the contract to the inferred type keeps the validator happy on
 * either path.
 */
function subjectRawType(df: DG.DataFrame, name: string): 'string' | 'int32' {
  const col = df.col(name);
  if (col === null) throw new Error(`subjectRawType: column '${name}' missing`);
  return col.type === DG.COLUMN_TYPE.INT ? 'int32' : 'string';
}

/**
 * Add a `Dose_mg` per-row column derived as `Dose * Wt` for the theoph
 * fixture (Dose is per-kg in the CSV; sci-comp consumes absolute mg). The
 * resulting column is Float64-backed via `Column.fromFloat64Array`, matching
 * the adapter's Float64 invariant (spec §7.2 / §10).
 *
 * Returns the same DataFrame for chaining.
 */
export function withDoseMg(df: DG.DataFrame): DG.DataFrame {
  const wt = df.col('Wt');
  const dosePerKg = df.col('Dose');
  if (wt === null || dosePerKg === null)
    throw new Error(`withDoseMg: theoph CSV missing 'Wt' or 'Dose'`);
  const n = df.rowCount;
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) out[i] = dosePerKg.get(i) * wt.get(i);
  df.columns.add(DG.Column.fromFloat64Array('Dose_mg', out));
  return df;
}

/**
 * Build the column contract for the theoph fixture (oral, per-subject dose
 * `Dose * Wt` mg, body weight in kg). Subject rawType is probed against the
 * loaded DataFrame so int / string autodetection both validate.
 */
export function theophContract(df: DG.DataFrame): ColumnContract {
  return {
    observations: {
      subjectId: {name: 'Subject', required: true, rawType: subjectRawType(df, 'Subject')},
      time: {name: 'Time', required: true, rawType: 'float64', units: 'h'},
      concentration: {name: 'conc', required: true, rawType: 'float64', units: 'mg/L'},
      dose: {name: 'Dose_mg', required: false, rawType: 'float64', units: 'mg'},
      bodyWeight: {name: 'Wt', required: false, rawType: 'float64', units: 'kg'},
    },
  };
}

/**
 * Scalars for the theoph fixture. Dose / route are read per-row from the
 * `Dose_mg` / no-route-column paths above; the scalars only matter as
 * fallbacks (and as the metadata `concentrationUnits / timeUnits / doseUnits`
 * the engine carries through `ProfileInputs`).
 */
export function theophScalars(): ProfileScalars {
  return {
    dose: 0,
    route: nca.ROUTE_PO,
    doseUnits: 'mg',
    concentrationUnits: 'mg/L',
    timeUnits: 'h',
    bodyWeight: null,
    bodyWeightUnit: 'kg',
    infusionDuration: null,
    lloq: 0,
  };
}

/**
 * Build the column contract for the indometh fixture (IV bolus, fixed
 * 25 mg dose carried as a scalar, no per-row dose / body-weight columns).
 */
export function indomethContract(df: DG.DataFrame): ColumnContract {
  return {
    observations: {
      subjectId: {name: 'Subject', required: true, rawType: subjectRawType(df, 'Subject')},
      time: {name: 'time', required: true, rawType: 'float64', units: 'h'},
      concentration: {name: 'conc', required: true, rawType: 'float64', units: 'mg/L'},
    },
  };
}

/** Scalars for the indometh fixture — 25 mg fixed dose, IV-bolus route. */
export function indomethScalars(): ProfileScalars {
  return {
    dose: 25,
    route: nca.ROUTE_IV_BOLUS,
    doseUnits: 'mg',
    concentrationUnits: 'mg/L',
    timeUnits: 'h',
    bodyWeight: null,
    bodyWeightUnit: 'kg',
    infusionDuration: null,
    lloq: 0,
  };
}

/**
 * Build the column contract for the IV-infusion fixture (FR-200 AC-D3). Same
 * Subject/time/conc shape as indometh; the route + infusion duration are
 * scalars (see {@link infusionScalars}).
 */
export function infusionContract(df: DG.DataFrame): ColumnContract {
  return {
    observations: {
      subjectId: {name: 'Subject', required: true, rawType: subjectRawType(df, 'Subject')},
      time: {name: 'time', required: true, rawType: 'float64', units: 'h'},
      concentration: {name: 'conc', required: true, rawType: 'float64', units: 'mg/L'},
    },
  };
}

/**
 * Scalars for the IV-infusion fixture — 100 mg dose over a 1 h infusion (AC-D3).
 * The non-null `infusionDuration` is what threads through `scaleProfileInputs`
 * into `computeNca`'s IV-infusion branch, triggering the −T_inf/2 MRT/Vss
 * correction. This is the only fixture that exercises that path.
 */
export function infusionScalars(): ProfileScalars {
  return {
    dose: 100,
    route: nca.ROUTE_IV_INFUSION,
    doseUnits: 'mg',
    concentrationUnits: 'mg/L',
    timeUnits: 'h',
    bodyWeight: null,
    bodyWeightUnit: 'kg',
    infusionDuration: 1,
    lloq: 0,
  };
}

/**
 * Build the column contract for the rat multi-route fixtures (03_rat_iv / po /
 * sc). Preclinical schema: per-row `Dose_mg_kg` (per-kg) + `BodyWeight_g`
 * (grams) → the extractor canonicalizes dose to absolute mg via BW; conc in
 * ng/mL → mg/L; time in h. Route is supplied as a scalar (each file is
 * single-route), so no per-row Route mapping. SubjectID is the R-coded id.
 */
export function ratContract(df: DG.DataFrame): ColumnContract {
  return {
    observations: {
      subjectId: {name: 'SubjectID', required: true, rawType: subjectRawType(df, 'SubjectID')},
      time: {name: 'Time_hr', required: true, rawType: 'float64', units: 'h'},
      concentration: {name: 'Conc_ng_mL', required: true, rawType: 'float64', units: 'ng/mL'},
      dose: {name: 'Dose_mg_kg', required: false, rawType: 'float64', units: 'mg/kg'},
      bodyWeight: {name: 'BodyWeight_g', required: false, rawType: 'float64', units: 'g'},
    },
  };
}

/**
 * Scalars for the rat fixtures. Dose / body weight are read per-row from the
 * `Dose_mg_kg` / `BodyWeight_g` columns; the scalars carry the source units
 * (`mg/kg`, `ng/mL`, `g`) that drive canonicalization. `route` is per-file.
 */
export function ratScalars(route: nca.RouteCode): ProfileScalars {
  return {
    dose: 0,
    route,
    doseUnits: 'mg/kg',
    concentrationUnits: 'ng/mL',
    timeUnits: 'h',
    bodyWeight: null,
    bodyWeightUnit: 'g',
    infusionDuration: null,
    lloq: 0,
  };
}

/**
 * Build the column contract for the multicompound batch-screen fixture
 * (`05_multicompound_batch_screen.csv`, VAL-01). Discovery-PK schema: 20
 * compounds × 3 rats, all IV bolus, per-row `Dose_mg_kg` (per-kg) and **no
 * body-weight column** → the per-kg CL/Vz path (nca-cl-vz-convention AD-2,
 * Ed obs #10). `SubjectID` is globally unique across compounds (R001..R060,
 * one compound each), so grouping by `SubjectID` yields the 60 per-subject
 * profiles the reference is keyed on — no compound dimension needed here.
 *
 * BLQ rows already carry `Conc_ng_mL = 0.000`, so (mirroring {@link ratContract})
 * neither `blqFlag` nor `lloq` is mapped — the trailing-zero trim (BUG-05 fix,
 * sci-comp 0.7.1 `collectMeasurable`) handles them exactly as in dataset 03.
 */
export function multicompoundContract(df: DG.DataFrame): ColumnContract {
  return {
    observations: {
      subjectId: {name: 'SubjectID', required: true, rawType: subjectRawType(df, 'SubjectID')},
      time: {name: 'Time_hr', required: true, rawType: 'float64', units: 'h'},
      concentration: {name: 'Conc_ng_mL', required: true, rawType: 'float64', units: 'ng/mL'},
      dose: {name: 'Dose_mg_kg', required: false, rawType: 'float64', units: 'mg/kg'},
    },
  };
}

/**
 * Scalars for the multicompound fixture — IV-bolus, mg/kg dose read per-row,
 * conc in ng/mL, time in h, **no body weight** (per-kg path). Mirrors
 * {@link ratScalars} minus the body-weight unit (BW is unmapped, not gram-scaled).
 */
export function multicompoundScalars(): ProfileScalars {
  return {
    dose: 0,
    route: nca.ROUTE_IV_BOLUS,
    doseUnits: 'mg/kg',
    concentrationUnits: 'ng/mL',
    timeUnits: 'h',
    bodyWeight: null,
    bodyWeightUnit: 'kg',
    infusionDuration: null,
    lloq: 0,
  };
}

/**
 * Read a numeric cell from a reference DataFrame column, mapping platform null
 * (`isNone`) to `null` so absent reference parameters surface as `null` rather
 * than a sentinel number.
 */
function refCell(df: DG.DataFrame, name: string, i: number): number | null {
  const c = df.col(name);
  if (c === null || c.isNone(i))
    return null;
  return c.get(i) as number;
}

/**
 * Load a `*_nca_datagrok_match.csv` PKNCA reference (the rat multi-route study
 * shipped CSV, not JSON) into the {@link FixtureFile} shape, keyed by
 * `SubjectID`. Re-parses FLOAT columns to Float64 (BUG-01) so the reference
 * values themselves are full precision. Maps the PKNCA dot-named columns to the
 * fixture parameter/provenance fields; `c0_extrapolated` is absent from these
 * references (→ null).
 */
export async function loadReferenceCsv(csvName: string): Promise<FixtureFile> {
  const csv = await grok.dapi.files.readAsText(`${FIXTURE_PATH_PREFIX}/${csvName}`);
  const df = DG.DataFrame.fromCsv(csv);
  reparseFloatColumnsToFloat64(df, csv);

  const sid = df.col('SubjectID');
  if (sid === null)
    throw new Error(`loadReferenceCsv: '${csvName}' missing SubjectID column`);
  const route = df.col('Route');
  const profiles: FixtureProfile[] = [];
  for (let i = 0; i < df.rowCount; i++) {
    profiles.push({
      profile_key: {
        subject: String(sid.get(i)),
        route: route !== null ? String(route.get(i)) : '',
      },
      parameters: {
        cmax: refCell(df, 'cmax', i),
        tmax: refCell(df, 'tmax', i),
        auclast: refCell(df, 'auclast', i),
        aucinf: refCell(df, 'aucinf.obs', i),
        pct_aucextrap: refCell(df, 'aucpext.obs', i),
        lambda_z: refCell(df, 'lambda.z', i),
        half_life: refCell(df, 'half.life', i),
        cl: refCell(df, 'cl.obs', i),
        vz: refCell(df, 'vz.obs', i),
      },
      provenance: {
        lambda_z_n_points: refCell(df, 'lambda.z.n.points', i),
        lambda_z_r_squared: refCell(df, 'r.squared', i),
        lambda_z_adj_r_squared: refCell(df, 'adj.r.squared', i),
        lambda_z_time_first: refCell(df, 'lambda.z.time.first', i),
        lambda_z_time_last: refCell(df, 'lambda.z.time.last', i),
        clast_obs: refCell(df, 'clast.obs', i),
        c0_extrapolated: null,
        // Resolved BY HEADER NAME, never by column index — these references are
        // not position-compatible with each other: `span.ratio` is field 21 in
        // the three rat files and field 18 in `multicompound_nca_per_subject.csv`
        // (measured 2026-08-12). A positional read would compare the wrong column
        // and then pass or fail for the wrong reason. `refCell` returns null when
        // the column is absent, so a reference without the oracle degrades to
        // "not asserted" rather than to a false comparison.
        span_ratio: refCell(df, 'span.ratio', i),
      },
    });
  }
  return {dataset: csvName, pknca_version: '0.12.1', profiles};
}

/**
 * Convenience wrapper that runs the full UC-01 pipeline on a DataFrame and
 * narrows the outcome to the `completed` variant. Throws if the validators
 * rejected the DataFrame — the reference-suite tests want the failure to
 * surface as a stack frame rather than a fall-through to "rowCount = 0".
 */
export function runPipelineOrThrow(
  df: DG.DataFrame,
  contract: ColumnContract,
  scalars: ProfileScalars,
  rules: nca.NcaRules,
  options: RunNcaPipelineOptions = {},
): {results: DG.DataFrame; outcome: Extract<NcaRunOutcome, {status: 'completed'}>} {
  const outcome = runNcaPipeline(df, contract, scalars, rules, options);
  if (outcome.status !== 'completed') {
    const detail = outcome.errors
      .map((e) => `${e.code}:${e.column ?? '-'}${e.row !== undefined ? `[${e.row}]` : ''}`)
      .join(', ');
    throw new Error(`runPipelineOrThrow: validation-failed — ${detail || '(no errors reported)'}`);
  }
  return {results: outcome.results, outcome};
}

/**
 * Build the column contract for the combined rat multi-route fixture
 * (`03_rat_multiroute.csv`: bare `IV` / `PO` / `SC` in one `Route` column).
 * Mirrors {@link ratContract} but adds the encoded RouteCode slot
 * ({@link DERIVED_ROUTE_CODE_COL}) — call {@link setupMultirouteRoute} first to
 * materialise that column (Option A: pane-side encoding, replicated here).
 */
export function ratMultirouteContract(df: DG.DataFrame): ColumnContract {
  return {
    observations: {
      subjectId: {name: 'SubjectID', required: true, rawType: subjectRawType(df, 'SubjectID')},
      time: {name: 'Time_hr', required: true, rawType: 'float64', units: 'h'},
      concentration: {name: 'Conc_ng_mL', required: true, rawType: 'float64', units: 'ng/mL'},
      dose: {name: 'Dose_mg_kg', required: false, rawType: 'float64', units: 'mg/kg'},
      route: {name: DERIVED_ROUTE_CODE_COL, required: false, rawType: 'int32'},
      bodyWeight: {name: 'BodyWeight_g', required: false, rawType: 'float64', units: 'g'},
    },
  };
}

/**
 * Replicate the mapping pane's Option-A route encoding for the multi-route
 * fixture: tolerantly normalize the `Route` string column to a derived int32
 * RouteCode column ({@link DERIVED_ROUTE_CODE_COL}) added to `df`, returning the
 * per-subject route warnings (pass into `runPipelineOrThrow` options).
 */
export function setupMultirouteRoute(df: DG.DataFrame): RunNcaPipelineOptions {
  const norm = normalizeRouteColumn(df.col('Route')!, df.col('SubjectID')!, {});
  if (norm.errors.length > 0)
    throw new Error(`setupMultirouteRoute: unexpected blocking errors: ${norm.errors.map((e) => e.code).join(', ')}`);
  if (df.col(DERIVED_ROUTE_CODE_COL) !== null) df.columns.remove(DERIVED_ROUTE_CODE_COL);
  df.columns.add(DG.Column.fromInt32Array(DERIVED_ROUTE_CODE_COL, norm.codes));
  return {routeWarningsBySubject: norm.warningsBySubject};
}

/**
 * Build a small synthetic observations DataFrame for adapter unit tests.
 * Constructs the columns from typed-array fragments so the storage layout
 * is predictable (avoids autodetection surprises).
 *
 * Caller supplies parallel arrays for subject/time/conc plus optional
 * dose / route / bodyWeight / blq / lloq.
 */
export interface SyntheticObsSpec {
  readonly subjects: ReadonlyArray<string>;
  readonly times: Float64Array;
  readonly concs: Float64Array;
  readonly doses?: Float64Array;
  readonly routes?: Int32Array;
  readonly bodyWeights?: Float64Array;
  readonly blq?: Uint8Array;
  readonly lloqs?: Float64Array;
}

export function makeSyntheticObs(spec: SyntheticObsSpec): DG.DataFrame {
  const n = spec.subjects.length;
  if (spec.times.length !== n || spec.concs.length !== n)
    throw new Error('makeSyntheticObs: subjects/times/concs length mismatch');

  const cols: DG.Column[] = [
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', Array.from(spec.subjects)),
    DG.Column.fromFloat64Array('Time', spec.times),
    DG.Column.fromFloat64Array('conc', spec.concs),
  ];
  if (spec.doses !== undefined)
    cols.push(DG.Column.fromFloat64Array('Dose_mg', spec.doses));
  if (spec.routes !== undefined)
    cols.push(DG.Column.fromInt32Array('Route', spec.routes));
  if (spec.bodyWeights !== undefined)
    cols.push(DG.Column.fromFloat64Array('Wt', spec.bodyWeights));
  if (spec.blq !== undefined) {
    const flag = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'BLQ',
      Array.from(spec.blq).map((b) => b !== 0));
    cols.push(flag);
  }
  if (spec.lloqs !== undefined)
    cols.push(DG.Column.fromFloat64Array('LLOQ', spec.lloqs));
  return DG.DataFrame.fromColumns(cols);
}

/**
 * Default 5-point synthetic profile (one subject, IV bolus, clean log-linear
 * decay) used by unit tests when the input data itself is not the
 * variable under test.
 */
export function singleSubjectSyntheticDf(): DG.DataFrame {
  return makeSyntheticObs({
    subjects: ['S1', 'S1', 'S1', 'S1', 'S1'],
    times: Float64Array.from([0, 1, 2, 4, 8]),
    concs: Float64Array.from([10, 5, 2.5, 0.625, 0.039]),
  });
}

export function singleSubjectContract(df: DG.DataFrame): ColumnContract {
  return {
    observations: {
      subjectId: {name: 'Subject', required: true, rawType: subjectRawType(df, 'Subject')},
      time: {name: 'Time', required: true, rawType: 'float64', units: 'h'},
      concentration: {name: 'conc', required: true, rawType: 'float64', units: 'mg/L'},
    },
  };
}

export function singleSubjectScalars(): ProfileScalars {
  return {
    dose: 100,
    route: nca.ROUTE_IV_BOLUS,
    doseUnits: 'mg',
    concentrationUnits: 'mg/L',
    timeUnits: 'h',
    bodyWeight: null,
    bodyWeightUnit: 'kg',
    infusionDuration: null,
    lloq: 0,
  };
}

/**
 * Fail if any numeric column of `df` carries a raw IEEE NaN — the permanent
 * write-side drift guard behind the DG null-sentinel migration
 * (`nca-dg-native-missing-values` F4 Arm 1 / AC-12 / AC-13).
 *
 * Lives here rather than in the detector test file so that all NINE result-frame
 * producers can be covered from wherever their output already exists, instead of
 * duplicating four fixtures into one file. A raw NaN in a DG column is invisible
 * to `isNone`, is counted as PRESENT by group-by / filters / export, and poisons
 * every statistic on its column.
 *
 * `±Infinity` is REPORTED, not failed — `f64Column` deliberately passes it
 * through (AC-3) and the `ivRef === 0` defect that produces it is a separate,
 * user-deferred gap (AC-12b).
 */
export function expectNoRawNan(df: DG.DataFrame, frameName: string): void {
  const entries = scanNumericSentinels(df);
  const nans = entries.filter((e) => e.kind === 'raw-nan');
  expect(nans.length, 0,
    `${frameName}: raw IEEE NaN in a DG column — invisible to isNone, counted as PRESENT by ` +
    `group-by/filter/export, and it poisons the column's statistics. ` +
    `Offenders: ${describeSentinelEntries(frameName, nans)}`);
  const infs = entries.filter((e) => e.kind === 'non-finite');
  if (infs.length > 0)
    console.log(`[AC-12b report-only] ${describeSentinelEntries(frameName, infs)}`);
}
