/**
 * NCA EMS migration — pre-migration goldens (AC-1.4, REGRESSION_PIN).
 *
 * The EMS migration rewrites the open path (S2/S3) and moves the pane's pre-engine assembly
 * into `assembleEngineInputs` / `replayAnalysis` (S8a). Both are acceptable only with
 * byte-identical outputs, so the outputs of the CURRENT code are pinned here BEFORE any of
 * those slices land. Each case drives the real pane (`mountAnalysisIntoView` → `triggerCompute`,
 * the same path as the play icon) and captures, from the compute callback:
 *  - the reportable results frame (every column, every row, full float64 precision);
 *  - the sparse composite table when the analysis is sparse (not in the callback context —
 *    read from the docked grid);
 *  - the summary cache the persistence layer stores (`buildSummaryCache`, fixed label + time);
 *  - the reconciled workflow state (overrides, ledgers, stage) and the emitted recipe.
 * Who / when fields (`user`, `at`, …) and the `_lz_*` service columns (analyst names, local
 * times) are masked; everything else is compared as a string.
 *
 * Generation: a case with no golden returns its capture as the test result and passes; read
 * it with `grok test --csv`. The separate "every case has a golden" test fails until every
 * capture is committed into `goldens/ems-pre-migration.json`, so a missing golden cannot pass
 * silently. The file records the `assemblyVersion` it was generated under (AC-8a.7).
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {loadSampleByName, mountAnalysisIntoView} from '../views/nca-launch';
import {normalizeNumericNaMarkers} from '../adapter/na-normalize';
import {liveWorkflowState, setPendingWorkflowState, SPARSE_PARAM_GRID_TAG} from '../views/nca-mapping-pane';
import type {ComputeCallbackCtx, PendingWorkflowState} from '../views/nca-mapping-pane';
import {buildSummaryCache} from '../adapter/summary-cache';
import {getNcaRepository} from '../adapter/nca-repository';
import type {AnalysisRecipe} from '../adapter';
import {ROUTE_STRING_TO_CODE} from '../adapter';
import GOLDENS from './goldens/ems-pre-migration.json';

const TIMEOUT = 180000;
// eslint-disable-next-line @typescript-eslint/no-explicit-any
type Bag = Record<string, any>;

/** One golden case: the fixture, an optional recipe / workflow / data edit, the compute mode. */
interface GoldenCase {
  readonly file: string;
  readonly recipe?: () => Promise<AnalysisRecipe>;
  readonly workflow?: () => Promise<PendingWorkflowState>;
  /** Edits the freshly loaded frame before the pane binds (a changed point, D13). */
  readonly mutate?: () => Promise<(df: DG.DataFrame) => void>;
  readonly silent?: boolean;
}

// ── canonical serialization ──────────────────────────────────────────────────

const NUMERIC = new Set<string>([DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.BIG_INT,
  DG.COLUMN_TYPE.QNUM]);
const MASKED_KEYS = new Set(['at', 'user', 'supersededAt']);

/** A cell as JSON: null for a DG-missing cell (isNone first — guardrails), the raw float64
 *  (JSON keeps the shortest round-trip form, so equal doubles serialize equally), 'NaN' /
 *  'Infinity' as strings, and getString for everything else. `_lz_*` cells hold analyst names
 *  and local times, so only their presence is kept. */
function cellJson(col: DG.Column, i: number): unknown {
  if (col.name.startsWith('_lz_'))
    return col.isNone(i) || col.getString(i) === '' ? '' : '<set>';
  if (col.isNone(i)) return null;
  if (NUMERIC.has(col.type)) {
    const v = Number(col.get(i));
    return Number.isFinite(v) ? v : String(v);
  }
  return col.getString(i);
}

function frameJson(df: DG.DataFrame | null | undefined): unknown {
  if (df == null) return null;
  return {rows: df.rowCount, cols: df.columns.toList().map((c) => {
    const values: unknown[] = [];
    for (let i = 0; i < df.rowCount; i++) values.push(cellJson(c, i));
    return {name: c.name, type: c.type, values};
  })};
}

/** Deep copy with the who / when fields replaced by a marker (null stays null). A string that
 *  holds a JSON object / array (a decision serialized into a summary cell) is masked inside and
 *  re-serialized — its embedded `at` otherwise made the workflow goldens differ on every run. */
function masked(x: unknown): unknown {
  if (typeof x === 'string' && (x.startsWith('{') || x.startsWith('['))) {
    try {
      return JSON.stringify(masked(JSON.parse(x)));
    } catch {
      return x; // not JSON after all — a plain value
    }
  }
  if (Array.isArray(x)) return x.map(masked);
  if (x !== null && typeof x === 'object') {
    // A summary key inventory of an `_lz_*` service column (analyst names, override times):
    // only its presence is kept, as `cellJson` does for those cells in the results frame.
    const rec = x as Record<string, unknown>;
    if (typeof rec.key === 'string' && rec.key.startsWith('_lz_') && Array.isArray(rec.values))
      return {...rec, values: '<lz>'};
    const out: Record<string, unknown> = {};
    for (const [k, v] of Object.entries(x as Record<string, unknown>)) {
      // The analyst's user id on an override / decision (ems-modification, IR-12) postdates these
      // goldens and is environment-specific: dropped here, pinned by the controller and storage tests.
      if (k === 'userId') continue;
      out[k] = MASKED_KEYS.has(k) && v != null ? '<masked>' : masked(v);
    }
    return out;
  }
  return x;
}

// ── the runner ───────────────────────────────────────────────────────────────

async function waitFor(cond: () => boolean, ms = 8000): Promise<boolean> {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) {
    if (cond()) return true;
    await DG.delay(100);
  }
  return cond();
}

/** Mount the pane on the fixture as the open path does, compute once, return the capture.
 *  `persisted`: the frame first goes through the store — written as a snapshot, read back, the
 *  snapshot then removed — so the capture is what a reopened stored analysis computes (AC-7.2). */
async function runCase(c: GoldenCase, persisted = false): Promise<string> {
  let df = await loadSampleByName(c.file);
  if (c.mutate) (await c.mutate())(df);
  if (persisted) {
    const repo = getNcaRepository();
    const written = await repo.writeSnapshot(df, new Date().toISOString());
    try {
      df = await repo.readSource(written.source);
    } finally {
      await repo.deleteSnapshot(written.source);
    }
  }
  normalizeNumericNaMarkers(df); // mirror the production open path (recipe-round-trip roundTrip)
  const recipe = c.recipe ? await c.recipe() : undefined;
  const workflow = c.workflow ? await c.workflow() : undefined;
  const view = grok.shell.addTableView(df);
  let ctx: ComputeCallbackCtx | null = null;
  try {
    const handle = mountAnalysisIntoView(view, df, {recipe, rules: recipe?.ncaRules,
      onComputed: (x) => {ctx = x;}});
    if (workflow) setPendingWorkflowState(view, workflow);
    await DG.delay(100); // deferred autodetect / sex auto-bind settle
    handle.triggerCompute(c.silent ? {silent: true} : undefined);
    const got = ctx as ComputeCallbackCtx | null;
    if (got === null)
      return JSON.stringify({status: 'blocked', recipe: masked(handle.getCurrentRecipe())});
    const bag = view.temp as Bag;
    if (got.recipe.modes.sparse.on) await waitFor(() => bag[SPARSE_PARAM_GRID_TAG] != null);
    const summary = buildSummaryCache(got.resultsDf, got.perProfile, got.recipe, 'golden', got.observations);
    return JSON.stringify({
      status: 'computed',
      recipe: masked(got.recipe),
      results: frameJson(got.resultsDf),
      sparse: frameJson(bag[SPARSE_PARAM_GRID_TAG]?.dataFrame),
      summary: masked(summary),
      workflow: masked(got.workflow),
    });
  } finally {
    for (const d of DG.Dialog.getOpenDialogs()) d.close();
    view.close();
  }
}

/** The recipe the pane autodetects on a fixture (mounted, never computed). */
async function autodetected(file: string): Promise<AnalysisRecipe> {
  const df = await loadSampleByName(file);
  normalizeNumericNaMarkers(df);
  const view = grok.shell.addTableView(df);
  try {
    const h = mountAnalysisIntoView(view, df, {});
    await DG.delay(100);
    return h.getCurrentRecipe();
  } finally {
    view.close();
  }
}

/** NONMEM files whose ID column autodetect does not bind: map the standard roles explicitly
 *  (as the R2-F1 route test in recipe-round-trip-tests does). */
function withNonmemIds(r: AnalysisRecipe): AnalysisRecipe {
  return {...r, columns: {...r.columns, subject: r.columns.subject ?? 'ID', time: r.columns.time ?? 'TIME',
    concentration: r.columns.concentration ?? 'DV'}};
}

/** UTF-8 → base64: the `grok test --csv` report rewrites `"` as `'`, so a JSON capture cannot
 *  travel through it verbatim. */
function toBase64(s: string): string {
  const bytes = new TextEncoder().encode(s);
  let bin = '';
  for (const b of bytes) bin += String.fromCharCode(b);
  return btoa(bin);
}

// ── the override / Exclude / vanished-point / D13 fixture (theoph) ───────────

/** Recorded once through the LIVE controller of a computed theoph view, so every key
 *  (`timeKey`, `concKey`, fingerprint) is the one the app itself writes. Subject A carries a
 *  force-out of its earliest λz point; subject B (Flagged) an Exclude λz decision. */
interface TheophWorkflow {
  readonly state: PendingWorkflowState;
  readonly overrideSubject: string;
  readonly overrideTimeKey: string;
  readonly excludedSubject: string;
  readonly subjectCol: string;
  readonly timeCol: string;
  readonly concCol: string;
}
let theophWf: Promise<TheophWorkflow> | null = null;

function recordTheophWorkflow(): Promise<TheophWorkflow> {
  theophWf ??= (async () => {
    const df = await loadSampleByName('theoph.csv');
    normalizeNumericNaMarkers(df);
    const view = grok.shell.addTableView(df);
    try {
      const handle = mountAnalysisIntoView(view, df, {});
      await DG.delay(100);
      handle.triggerCompute();
      const bag = view.temp as Bag;
      if (!await waitFor(() => bag['ncaController']?.results != null))
        throw new Error('theoph did not compute');
      const ctrl = bag['ncaController'];
      const recipe = handle.getCurrentRecipe();
      const {units, fingerprints} = liveWorkflowState(view);
      const excluded = units.find((u) => u.allowedActions.includes('exclude'));
      if (excluded === undefined) throw new Error('no Flagged theoph profile to exclude');
      // The subject with the most λz points, other than the excluded one, gets the override.
      const fits = ctrl.fitRowsBySubject() as Map<string, Set<number>>;
      let target = '';
      for (const [s, rows] of fits)
        if (s !== excluded.key && rows.size > (fits.get(target)?.size ?? 0)) target = s;

      const earliest = Math.min(...fits.get(target)!);
      ctrl.overrideLambdaZPoints(target, [{rowIdx: earliest, action: 'force-out'}]);
      ctrl.decide(excluded.key, 'exclude', [{code: 'EXC_ERRATIC', text: 'erratic'}],
        fingerprints.get(excluded.key) ?? '');
      const overrides = ctrl.getLzOverrides();
      return {
        state: {lzOverrides: overrides, reviewDecisions: ctrl.getReviewLedgers()},
        overrideSubject: target, overrideTimeKey: overrides[0].timeKey, excludedSubject: excluded.key,
        subjectCol: recipe.columns.subject!, timeCol: recipe.columns.time!, concCol: recipe.columns.concentration!,
      };
    } finally {
      view.close();
    }
  })();
  return theophWf;
}

/** Scales the concentration of the matching rows (canonical getString keys). */
function scaleConc(df: DG.DataFrame, w: TheophWorkflow, subject: string,
  pick: (timeKeys: string[]) => string, factor: number): void {
  const s = df.col(w.subjectCol)!; const t = df.col(w.timeCol)!; const c = df.col(w.concCol)!;
  const rows: number[] = [];
  for (let i = 0; i < df.rowCount; i++) if (s.getString(i) === subject) rows.push(i);
  const key = pick(rows.map((i) => t.getString(i)));
  for (const i of rows) if (t.getString(i) === key) c.set(i, Number(c.get(i)) * factor);
}

// ── the cases ────────────────────────────────────────────────────────────────

const CASES: Record<string, GoldenCase> = {
  // (a) round-trip modes
  'dense theoph': {file: 'theoph.csv'},
  'dense indometh': {file: 'indometh.csv'},
  'nonmem multiroute (explicit ROUTE)': {file: '03_rat_multiroute_nonmem.csv'},
  'sparse destructive': {file: '04_mouse_sparse_destructive.csv'},
  'multicompound batch (silent through per-kg)': {file: '05_multicompound_batch_screen.csv', silent: true},
  // (c) assembler heuristic branches
  'nonmem day column + dose-group pooling': {file: '06_rat_tk28d_nonmem.csv'},
  'nonmem per-day reset': {
    file: '08_perday_reset_nonmem.csv',
    recipe: async () => withNonmemIds(await autodetected('08_perday_reset_nonmem.csv')),
  },
  'sampling-type pattern (two compounds)': {file: '04c_mouse_sparse_two_compounds.csv'},
  'nominal-time + animal fallback': {
    file: '04_mouse_sparse_destructive.csv',
    recipe: async () => {
      const r = await autodetected('04_mouse_sparse_destructive.csv');
      return {...r, modes: {...r.modes,
        sparse: {...r.modes.sparse, on: true, nominalTimeCol: null, animalIdCol: null}}};
    },
  },
  'NA markers: blank body weight': {file: '03_rat_iv_partial_bw.csv'},
  // (c) interactive preconditions
  // Subject mapped explicitly: autodetect leaves it null on this file, and an unmapped subject
  // blocks before the route gate is reached — the golden would then pin the wrong guard.
  'route confirm: blocks a user compute': {
    file: '04_iv_infusion_nonmem.csv',
    recipe: async () => withNonmemIds(await autodetected('04_iv_infusion_nonmem.csv')),
  },
  'route confirm: confirmed route computes': {
    file: '04_iv_infusion_nonmem.csv',
    recipe: async () => {
      const r = withNonmemIds(await autodetected('04_iv_infusion_nonmem.csv'));
      return {...r, modes: {...r.modes, confirmedRoute: ROUTE_STRING_TO_CODE['IV-infusion']}};
    },
  },
  // The pre-Run per-kg confirm was retired: a user compute of a /kg dose with no body weight
  // now computes, per kg as everywhere.
  'per-kg dose, no BW: a user compute computes (confirm retired)': {file: '05_multicompound_batch_screen.csv'},
  'sex unassigned: blocks': {file: '06_rat_tk28d_nonmem_sexcoded.csv'},
  // Positive control: the same file with its codes assigned computes, so the block above is
  // attributable to the sex gate and nothing else.
  'sex codes assigned: computes': {
    file: '06_rat_tk28d_nonmem_sexcoded.csv',
    recipe: async () => {
      const r = await autodetected('06_rat_tk28d_nonmem_sexcoded.csv');
      return {...r, columns: {...r.columns, sexCol: 'Sex', sexValueMap: {'1': 'M', '2': 'F'}}};
    },
  },
  'duplicate mapping: blocks': {
    file: 'theoph.csv',
    recipe: async () => {
      const r = await autodetected('theoph.csv');
      return {...r, columns: {...r.columns, bodyWeightSource: 'column', bodyWeightCol: r.columns.concentration}};
    },
  },
  'R_acc below minR: blocks': {
    file: 'theoph.csv',
    recipe: async () => ({...await autodetected('theoph.csv'), acceptanceRSquared: 0.5}),
  },
  // (b) λz overrides and review decisions
  'override + Exclude replay': {
    file: 'theoph.csv', workflow: async () => (await recordTheophWorkflow()).state,
  },
  'override whose point changed (reverts to auto)': {
    file: 'theoph.csv',
    workflow: async () => (await recordTheophWorkflow()).state,
    mutate: async () => {
      const w = await recordTheophWorkflow();
      return (df) => scaleConc(df, w, w.overrideSubject, () => w.overrideTimeKey, 1.5);
    },
  },
  'decision on a changed fit (D13 auto-reset)': {
    file: 'theoph.csv',
    workflow: async () => (await recordTheophWorkflow()).state,
    mutate: async () => {
      const w = await recordTheophWorkflow();
      // The last sample sits in every terminal window, so changing it changes the fit.
      return (df) => scaleConc(df, w, w.excludedSubject, (keys) => keys[keys.length - 1], 0.6);
    },
  },
};

const STORED = (GOLDENS as {cases: Record<string, string>}).cases;

/** The deepest first differing path and both values there, to make a failure readable (a
 *  top-level section name alone cannot tell a changed number from changed metadata). */
function firstDiff(actual: string, golden: string): string {
  const walk = (a: unknown, g: unknown, path: string): string | null => {
    if (JSON.stringify(a) === JSON.stringify(g)) return null;
    if (a !== null && g !== null && typeof a === 'object' && typeof g === 'object') {
      const keys = Array.isArray(a) && Array.isArray(g) ?
        Array.from({length: Math.max(a.length, g.length)}, (_, i) => String(i)) :
        [...new Set([...Object.keys(g as object), ...Object.keys(a as object)])];
      for (const k of keys) {
        const d = walk((a as Record<string, unknown>)[k], (g as Record<string, unknown>)[k], `${path}.${k}`);
        if (d !== null) return d;
      }
    }
    const cut = (x: unknown): string => (JSON.stringify(x) ?? 'undefined').slice(0, 200);
    return `${path}: actual ${cut(a)} vs golden ${cut(g)}`;
  };
  return walk(JSON.parse(actual), JSON.parse(golden), '$') ?? '(none)';
}

category('NCA EMS: pre-migration goldens', () => {
  for (const [name, c] of Object.entries(CASES)) {
    test(name, async () => {
      const actual = await runCase(c);
      const golden = STORED[name];
      if (golden === undefined) return `B64:${toBase64(actual)}`; // generation run: read it from --csv
      expect(actual === golden, true,
        `${name}: output differs from the pre-migration golden at ${firstDiff(actual, golden)}`);
      return 'matches golden';
    }, {timeout: TIMEOUT});
  }

  test('every case has a golden', async () => {
    const missing = Object.keys(CASES).filter((n) => STORED[n] === undefined);
    expect(missing.length, 0, `cases without a committed golden: ${missing.join(', ')}`);
    const stale = Object.keys(STORED).filter((n) => !(n in CASES));
    expect(stale.length, 0, `goldens without a case: ${stale.join(', ')}`);
  });
});

// AC-7.2: the same goldens after a trip through the store — a stored analysis, reopened, must
// compute exactly what the pre-migration code computed on the original file.
category('NCA EMS: persisted goldens', () => {
  for (const [name, c] of Object.entries(CASES)) {
    test(name, async () => {
      const golden = STORED[name];
      expect(golden !== undefined, true, `${name}: no committed golden`);
      const actual = await runCase(c, true);
      expect(actual === golden, true,
        `${name}: the stored round trip differs from the golden at ${firstDiff(actual, golden)}`);
      return 'matches golden after the store';
    }, {timeout: TIMEOUT});
  }
});
