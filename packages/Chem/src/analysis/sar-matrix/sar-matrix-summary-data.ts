/* The Summary tab's computation: one walk over every matrix that produces everything the tab
   shows. Kept apart from the rendering because it touches no DOM and the panel's own helpers never
   reach into it: a SummaryCollector is handed the viewer and the fold tier, and returns SummaryData. */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../../package';
import {AdditiveFit, fitAdditiveFromTriples} from './sar-matrix-assemble';
import {median} from './sar-matrix-decompose';
import {fitRoleEffects, RoleFit} from './sar-matrix-role-fit';
import {logSarTime, SarMatrix, SarMatrixCell} from './sar-matrix-types';
import {MatrixCellRef} from './sar-matrix-ui-common';

export const SUM_ROWS = 3;

export const SUM_POOL = 24;

/** Below this the leave-one-out fit does not support acting on a prediction. */
export const TRUST_R2 = 0.5;

/** A prediction resting on fewer measured compounds than this on either axis is extrapolation. */
export const MIN_SUPPORT = 3;

/** Measured pairs a pooled swap needs, and lineages it must span, before it is worth a row. */
export const SWAP_MIN_PAIRS = 3;

export const SWAP_MIN_SERIES = 2;

/** Measured cells kept per row when enumerating swaps; both extremes survive the trim. */
export const SWAP_ROW_CAP = 32;

/** Lineages an R-group must win in before it is ranked; two is the mean of two numbers. */
export const RGROUP_MIN_SERIES = 3;

/** Cross-validatable cells a fit needs before its R² is quotable as "best validated". */
export const BEST_FIT_MIN_N = 8;

/** Matrices a per-core outcome strip can carry before a row of squares stops being a mark. */
export const STRIP_SLOTS = 12;

/** Reliable losers shown under the winners; a loss over two lineages is not a rule, so the tier the
 *  winners get is not offered here. */
export const LOSER_ROWS = 2;

/** Rows the analog list holds. A chemist does not browse ten thousand; two hundred is more than a
 *  quarter's synthesis and fits one frame of rendered structures. */
export const ANALOG_LIST_MAX = 200;

export interface SummaryRow {
  /** What the pool deduplicates on; held rather than re-derived, since only the pool knows it. */
  key: string;
  matrix: SarMatrix;
  ri: number;
  ci: number;
  /** Always higher-is-better, so one comparator orders every pool. */
  score: number;
}

/** What the summary reads back from the viewer, and the things it asks the viewer to do. */
export interface SummaryHost {
  readonly matrices: SarMatrix[];
  /** Lineage root id per matrix, index-aligned with `matrices` — the one hierarchy the navigator uses. */
  readonly matrixRoots: string[];
  /** Displayed tier per matrix (the L1/L2/L3 on the navigator cards), index-aligned with `matrices`. */
  readonly matrixTiers: number[];
  /** Rows carrying an activity value at all. A compound the assay never reached and one no core could
   *  group are both outside every matrix, and only this separates them. */
  readonly assayedCount: number;
  /** Whether SAR transfer detection has run, and what it found. Detection is lazy and lives on its own
   *  tab, so the landing screen reports its state rather than triggering it. */
  readonly transferSummary: {scanned: boolean, count: number};
  readonly higherIsBetter: boolean;
  readonly scalingLabel: string;
  /** Whether a difference in activity units is already a log ratio. A raw column declared
   *  higher-is-better is a precomputed pIC50: no transform applied, but the numbers are logs. */
  readonly activityIsLog: boolean;
  readonly activityColumnName: string;
  /** Rows of the host table — the denominator the coverage line is a fraction of. */
  readonly hostRowCount: number;
  /** Assayed values the chosen scaling cannot represent; they are in no matrix and are not untested. */
  readonly unscalableCount: number;
  /** The fragment column every matrix varies, or null when these came from fragmentation. Only when
   *  it is set does one substituent label mean the same thing in two series. */
  readonly axisRole: string | null;
  /** The column the cores came from, or null when they came from fragmentation. A degrader set names
   *  it Linker, and "Best core" over an unnamed scaffold is the same finding nobody can act on. */
  readonly coreRole: string | null;
  /** Whether a matrix IS one core, so its fitted mean compares cores rather than groupings. */
  readonly coresAreSeries: boolean;
  /** Fragment columns that could run across the top, the current axis included. */
  readonly roleColumns: string[];
  setColumnAxis(name: string): void;
  /** Select the compounds carrying one value of one component column, in the table the analysis ran on. */
  selectRoleValue(role: string, value: string): void;
  readonly predictVirtual: boolean;
  readonly predictUnmeasured: boolean;
  readonly computing: boolean;
  noMatricesMessage(): string;
  formatActivity(value: number): string;
  cellIdText(cell: SarMatrixCell): string | null;
  observedNeighbours(matrix: SarMatrix, ri: number, ci: number): number;
  revealCell(matrix: SarMatrix, ri: number, ci: number, position?: string): void;
  revealMatrix(matrix: SarMatrix): void;
  showTab(name: string): void;
  addCellsToMakeList(cells: MatrixCellRef[], emptyMessage: string): void;
}

/** Both extremes by value once a group offers more pairs than the cap allows. A cap taken in the
 *  group's own order would keep the commonest substituents and truncate away the rare one that jumped
 *  two logs; this loses only mid-range pairs. */
export function keepExtremes<T>(items: T[], valueOf: (item: T) => number): T[] {
  if (items.length <= SWAP_ROW_CAP)
    return items;
  const half = SWAP_ROW_CAP >> 1;
  const sorted = [...items].sort((a, b) => valueOf(a) - valueOf(b));
  return [...sorted.slice(0, half), ...sorted.slice(sorted.length - half)];
}

export function supportOf(row: SummaryRow): number {
  return row.matrix.cells[row.ri][row.ci].support ?? 0;
}

/** Ties end at `matrix.id`: fragments arrive in worker-completion order, so anything resolved by
 *  index or Map order would make the cards depend on scheduling. */
export function finerSeries(a: SummaryRow, b: SummaryRow): boolean {
  return a.matrix.level !== b.matrix.level ? a.matrix.level < b.matrix.level : a.matrix.id < b.matrix.id;
}

export function betterSupported(a: SummaryRow, b: SummaryRow): boolean {
  const sa = supportOf(a);
  const sb = supportOf(b);
  return sa !== sb ? sa > sb : finerSeries(a, b);
}

export function betterEvidenced(a: SummaryRow, b: SummaryRow): boolean {
  const sa = supportOf(a);
  const sb = supportOf(b);
  if (sa !== sb)
    return sa > sb;
  const ra = a.matrix.confidence?.r2 ?? Number.NEGATIVE_INFINITY;
  const rb = b.matrix.confidence?.r2 ?? Number.NEGATIVE_INFINITY;
  return ra !== rb ? ra > rb : a.matrix.id < b.matrix.id;
}

/** A bounded best-of pool that deduplicates on insert. Deduplicating a pool of cells afterwards can
 *  leave it shorter than the card has room for, and the same compound occupies a cell in every tier
 *  that folded it. */
export class TopList {
  private readonly rows: SummaryRow[] = [];
  private readonly byKey = new Map<string, SummaryRow>();
  /** Keys that cleared every gate and are not in the pool for want of room. A set, not a counter: one
   *  structure is offered by every tier that folded its core, and a count would report it once per
   *  tier under a label the reader hears as structures. */
  private readonly missed = new Set<string>();

  /** `better` picks which occurrence of one key to keep; the score alone cannot, since a compound
   *  carries the same activity wherever it appears. */
  constructor(private readonly better: (a: SummaryRow, b: SummaryRow) => boolean,
    private readonly cap = SUM_POOL) {}

  /** The whole pool, best first — what a browsable list shows, as against the few rows a card has
   *  room for. */
  get all(): readonly SummaryRow[] {
    return this.rows;
  }

  /** How many further structures qualified. A capped list presented as a total states the cap as the
   *  finding, so every screen printing `all.length` has to be able to name this too. */
  get overflow(): number {
    return this.missed.size;
  }

  /** Drop the keys a late-arriving fact disqualifies. The structures the dataset already holds are
   *  only fully known once every matrix has been walked, and a pool offered to during that walk cannot
   *  have tested against the complete set. */
  prune(drop: (key: string) => boolean): void {
    for (let i = this.rows.length - 1; i >= 0; i--) {
      if (!drop(this.rows[i].key))
        continue;
      this.byKey.delete(this.rows[i].key);
      this.rows.splice(i, 1);
    }
    // The squeezed-out keys are disqualified by the same fact, and counting one of them as a structure
    // the list has no room for would offer to make something the dataset already holds.
    for (const key of this.missed) {
      if (drop(key))
        this.missed.delete(key);
    }
  }

  /** Primitive parameters: a full pool rejects a weak cell before anything is allocated. */
  offer(key: string, score: number, matrix: SarMatrix, ri: number, ci: number): void {
    const held = this.byKey.get(key);
    if (held === undefined && this.rows.length === this.cap && score <= this.rows[this.cap - 1].score) {
      this.missed.add(key);
      return;
    }
    const row: SummaryRow = {key, matrix, ri, ci, score};
    if (held !== undefined) {
      if (!this.better(row, held))
        return;
      this.rows.splice(this.rows.indexOf(held), 1);
    }
    let at = 0;
    while (at < this.rows.length && this.rows[at].score >= score)
      at++;
    this.rows.splice(at, 0, row);
    this.byKey.set(key, row);
    this.missed.delete(key);
    if (this.rows.length > this.cap) {
      const out = this.rows.pop()!;
      this.byKey.delete(out.key);
      this.missed.add(out.key);
    }
  }

  /** One row per series first — three analogs of one chemotype is one finding, not three — then top
   *  up from what the cap held back, so a run with fewer series than rows still fills the card.
   *  `skip` runs inside both passes: a card that filters afterwards renders short while the pool
   *  still holds rows it could have shown. */
  take(skip?: (row: SummaryRow) => boolean): SummaryRow[] {
    const out: SummaryRow[] = [];
    const seen = new Set<string>();
    for (const row of this.rows) {
      if (out.length === SUM_ROWS)
        break;
      if (seen.has(row.matrix.id) || skip?.(row))
        continue;
      seen.add(row.matrix.id);
      out.push(row);
    }
    for (const row of this.rows) {
      if (out.length === SUM_ROWS)
        break;
      if (!out.includes(row) && !skip?.(row))
        out.push(row);
    }
    return out;
  }
}

/** One measured cell of the row being walked, kept only long enough to enumerate that row's pairs. */
export interface RowCell {
  ci: number;
  value: number;
  molIdx: number;
}

/** A measured swap pooled across series: one R-group exchanged inside one row, everything else
 *  identical. `from`/`to` are the lexically ordered fragment pair, so a swap cannot split into two
 *  mirror pools oriented by an arbitrary column index; the display direction is chosen at render. */
export interface SwapPool {
  from: string;
  to: string;
  n: number;
  sum: number;
  min: number;
  max: number;
  /** Instances where `from → to` improved potency; the reverse count is `n - nUp`. */
  nUp: number;
  roots: Set<string>;
  /** Unordered molIdx pairs already counted — one compound pair recurs in every tier that folded it. */
  seen: Set<string>;
  best: {matrix: SarMatrix, ri: number, ciFrom: number, ciTo: number, delta: number} | null;
  sampled: boolean;
}

/** One matrix where a substituent came first at its position, under the fitted model. */
export interface RGroupWin {
  matrix: SarMatrix;
  ri: number;
  ci: number;
  root: string;
  /** Measured cells behind the column — the tiebreak where no reference is recorded to compare on. */
  n: number;
  /** `dir * (colEffect[c] − colEffect[reference])` — gauge-free, so this one pools. */
  refDelta: number | null;
  fitHolds: boolean;
}

export interface RGroupRow {
  subst: string;
  win: RGroupWin;
  /** Lineages won after the trust gate, and before it. */
  k: number;
  m: number;
  /** Lineages that carried this substituent as a comparable column at all — the denominator a reader
   *  hears in "best in 5 of 6". Wins alone would read as attempts and flatter every row. */
  tried: number;
  magnitude: number | null;
  lo: number;
  hi: number;
  magnitudeSeries: number;
}

/** One matrix's columns keyed by substituent, for the per-core strip. Only the columns the pooled
 *  comparison also counts are recorded, so the strip and the coverage sentence beneath it describe one
 *  set. `refDelta` is null only where the series records no reference substituent, and 0 where the
 *  column IS that reference. */
export type StripCols = Map<string, {ci: number, n: number, refDelta: number | null}>;

/** What the R-group card accumulates over the one fit per matrix. */
export interface RGroupAcc {
  wins: Map<string, RGroupWin[]>;
  losses: Map<string, RGroupWin[]>;
  tried: Map<string, Set<string>>;
}

/** Everything one matrix contributes, extracted from its fit and its cells; the fit arrays are
 *  discarded with it so 345 of them are never live at once. */
export interface SeriesStat {
  matrix: SarMatrix;
  root: string;
  tier: number;
  /** Distinct measured molecules from the walk — never `realCount`, which is captured pre-prune. */
  cpd: number;
  realCells: number;
  totalCells: number;
  /** Cells the decomposition cannot express, so the fill denominator is not `rows × columns`. Zero in
   *  a matrix whose assembly never marked them, which is NOT the same as a complete grid. */
  impossibleCells: number;
  lo: number;
  hi: number;
  best: {ri: number, ci: number, value: number} | null;
  /** The fit's own count-weighted mean. Unlike the centred effects it is one activity column under one
   *  scaling, so it is the one fitted quantity comparable between matrices. */
  typical: number;
  /** Most potent prediction that passes the trust gate — the model's own ceiling for this series. */
  bestVirtual: {ri: number, ci: number, value: number} | null;
  /** Unfilled cells this series carries a prediction for. With the matrix's own R2 it says WHY there
   *  is nothing to aim at. */
  virtualCells: number;
  trusted: number;
  converged: boolean;
  /** Filled only for the matrices a strip can draw; null everywhere else. */
  stripCols: StripCols | null;
  /** Range of the fitted substituent and core effects — centred, same units, same fit, so the two
   *  are comparable within this matrix and the verdict pools as a count. */
  colRange: number | null;
  rowRange: number | null;
  bestRow: {ri: number, effect: number, n: number} | null;
}

export interface StartRow {
  primary: SeriesStat;
  /** Indices into the fixed reason-lane slots, not phrases: the lane's shape is what a reader
   *  compares down the column. */
  reasons: number[];
}

export interface SummaryData {
  compounds: number;
  untested: number;
  measuredCells: number;
  minObserved: number | null;
  maxObserved: number | null;
  trustedCells: number;
  trustedStructures: number;
  /** Predictions that pass the same gate but carry no structure: the core holds an attachment point
   *  none of the picked fragment columns fills, so nothing can be completed over it. */
  trustedNoStructure: number;
  /** Typical leave-one-out prediction error across series — model error, never an assay σ. Null when
   *  no series has a cross-validated fit, so it is never printed as zero. */
  modelError: number | null;
  fitHolds: number;
  unchecked: number;
  lowR2: SarMatrix[];
  lowR2Virtual: number;
  families: number;
  axisRole: string | null;
  coreRole: string | null;
  coresAreSeries: boolean;
  /** One additive fit over every role column at once; null outside fragment-columns mode, where a
   *  substituent label is local to its own series and one pooled offset would average unrelated
   *  quantities. */
  roleFit: RoleFit | null;
  /** Best measured cell per value of each role column, so a leaderboard row lands on a compound.
   *  Keyed by role name and in the order the fit was given the roles, which the fit's own ranking
   *  discards. */
  roleBest: Map<string, Map<string, MatrixCellRef & {value: number}>>;
  tierCounts: {tier: number, n: number}[];
  /** Series whose additive fit stopped short of its tolerance. They are left out of every pooled
   *  R-group comparison, so an empty leaderboard has to be able to name this as the cause. */
  nonConverged: number;
  series: SeriesStat[];
  startHere: StartRow[];
  swaps: SwapPool[];
  /** The same, per component column, grouped straight off the component values rather than off one
   *  matrix's columns — so every component's measured pairs are in hand at once and none of them costs
   *  a rebuild. Empty outside fragment-columns mode, where a substituent label is local to its series. */
  swapsByRole: Map<string, SwapPool[]>;
  /** Whether any swap candidate existed at all — "none clears the gate" is a different sentence from
   *  "no row carries three measured substituents". */
  swapCandidates: number;
  /** The largest single measured move anywhere, so "nothing pooled clears the gate" can still name
   *  what the data does hold. */
  swapBest: SwapPool | null;
  rgroups: RGroupRow[];
  rgroupsThin: RGroupRow[];
  rgroupLosers: RGroupRow[];
  test: TopList;
  /** Virtual analogs ranked on predicted gain in their own series' leave-one-out error, deduplicated
   *  on structure — tiers overlap, so one molecule is proposed by several matrices. */
  analogs: TopList;
  /** Same gates bar the cross-validatable-cell floor: a tiny fit's error is not a scale to divide by,
   *  so these are listed without a rank rather than mixed into one. */
  analogsThin: TopList;
  /** The best of what the trust gate turned down, by raw gain over its own series' best measured
   *  compound. Shown only when the two lists above are empty, so the screen always names a next
   *  structure rather than reporting that nothing qualifies. */
  analogsAny: TopList;

  /** Why a prediction that passes support and fit is still not offered. Never one total: they call for
   *  different things from the reader, and only `withheldFitFails` has a list to open — a fit that was
   *  never checked and a fit that stopped short of its tolerance appear nowhere. */
  withheldBelowError: number;
  withheldThinSupport: number;
  withheldFitFails: number;
  withheldUnchecked: number;
  withheldNotConverged: number;
  /** Predicted structures the dataset already holds — an assay plate, not a synthesis. Structures, so
   *  it stands beside the analog count rather than multiplying it by the tiers that folded the core. */
  alreadyHeld: number;
  /** Structures that cleared every gate and are past the list's row cap, so "N worth making" can say
   *  what it is the top N of. */
  analogOverflow: number;
}

/** One half of a measured pair: the component value it carries, what it measured, and which compound it
 *  is — the three things a pool needs from each side, whether the side came off a matrix cell or off a
 *  row of component columns. */
export interface SwapSide {
  value: string;
  activity: number;
  mol: number;
}

/** The most potent measured cell of one column, so a leaderboard row lands on a compound. */
export function bestMeasuredRow(matrix: SarMatrix, ci: number, dir: number): number {
  let ri = -1;
  for (let r = 0; r < matrix.rows.length; r++) {
    const cell = matrix.cells[r][ci];
    if (cell.kind === 'real' && cell.value !== null &&
      (ri < 0 || dir * cell.value > dir * matrix.cells[ri][ci].value!))
      ri = r;
  }
  return ri;
}

/** One pass over the matrices at one fold tier. */
export class SummaryCollector {
  constructor(private readonly host: SummaryHost, private readonly tierFilter: number | null) {}

  /**
   * One pass over every cell of every matrix, feeding the totals, the two pools and the per-series
   * statistics.
   *
   * The observed cells are collected as triples on the way past and the additive fit is run from
   * those, so the fit does not walk the grid a second time. The buffers are reused across matrices:
   * 345 fresh sets would trade the scan cost for GC cost, and nothing caps the SUM of cells across
   * matrices — only each matrix.
   */
  collect(): SummaryData {
    const t0 = performance.now();
    const host = this.host;
    const dir = host.higherIsBetter ? 1 : -1;
    const log = host.activityIsLog;
    const roots = host.matrixRoots;
    const tiers = host.matrixTiers;

    const compounds = new Set<number>();
    const untested = new Set<number>();
    const trustedStructures = new Set<string>();
    // Assembled from the source table's molecule column while a virtual SMILES comes from RDKit's
    // linker, so the two canonicalisations need not agree: this catches matches, and is not a proof of
    // novelty — which is why a listed analog says "not in this dataset" and never "never made".
    const measuredStructures = new Set<string>();
    const test = new TopList(betterSupported);
    const analogs = new TopList(betterEvidenced, ANALOG_LIST_MAX);
    const analogsThin = new TopList(betterEvidenced, ANALOG_LIST_MAX);
    // Everything the gate turned down, ranked on raw gain. A screen that answers "nothing qualifies"
    // and stops leaves the reader with no next structure at all, when what the gate actually withheld
    // is the confidence, not the ranking.
    const analogsAny = new TopList(betterEvidenced, ANALOG_LIST_MAX);
    const swaps = new Map<string, SwapPool>();
    const acc: RGroupAcc = {wins: new Map(), losses: new Map(), tried: new Map()};
    const series: SeriesStat[] = [];
    const lowR2: SarMatrix[] = [];
    // A strip needs a fixed partner axis and few enough slots to read as a mark; fragmented matrices
    // overlap through tiers and have neither.
    const wantStrip = host.axisRole !== null && host.matrices.length <= STRIP_SLOTS;

    let measuredCells = 0;
    let trustedNoStructure = 0;
    let trustedCells = 0;
    let withheldBelowError = 0;
    let withheldThinSupport = 0;
    let withheldFitFails = 0;
    let withheldUnchecked = 0;
    let withheldNotConverged = 0;
    const alreadyHeld = new Set<string>();
    let lowR2Virtual = 0;
    let unchecked = 0;
    let fitHolds = 0;
    let nonConverged = 0;
    let sampledRows = 0;
    let swapCandidates = 0;
    let minObserved = Infinity;
    let maxObserved = -Infinity;

    const obsRow: number[] = [];
    const obsCol: number[] = [];
    const obsVal: number[] = [];
    const mols = new Set<number>();
    const rowCells: RowCell[] = [];
    // Admitting an analog on its gain needs this matrix's best measured cell, which is final only when
    // its row loop ends — so the candidates wait here rather than being offered inside it. Reused
    // across matrices: one triple of arrays, never 345 of them.
    const riBuf: number[] = [];
    const ciBuf: number[] = [];
    const valBuf: number[] = [];
    /** The same, for the candidates the trust gate turned down — they still need this matrix's best
     *  measured cell to carry a gain. */
    const riAlt: number[] = [];
    const ciAlt: number[] = [];
    const valAlt: number[] = [];

    // `roleColumns` excludes the core column, so the core role is prepended: left out, the report omits
    // the one role the request names first. Empty outside fragment-columns mode, which is what keeps the
    // fit off a decomposition whose labels are local to one series.
    const roleNames = host.axisRole === null ? [] : [host.coreRole!, ...host.roleColumns];
    // Where each role's value lives, decided once: the core on the row, the axis on the column, and
    // every other role folded into the row identity under its own column name.
    const roleKind = roleNames.map((name) => name === host.coreRole ? 0 : name === host.axisRole ? 1 : 2);
    const roleValues: string[][] = roleNames.map(() => []);
    const roleActivity: number[] = [];
    const roleMol: number[] = [];
    const roleBest = roleNames.map(() => new Map<string, MatrixCellRef & {value: number}>());

    for (let mi = 0; mi < host.matrices.length; mi++) {
      // Reading at one tier: skipping the others here is what keeps every ranking below — components,
      // swaps, cores, Start here, Worth making — on the same set of series, since all of them are built
      // from this one walk.
      if (this.tierFilter !== null && tiers[mi] !== this.tierFilter)
        continue;
      const matrix = host.matrices[mi];
      const nRows = matrix.rows.length;
      const nCols = matrix.columns.length;
      const r2 = matrix.confidence?.r2 ?? null;
      const holds = r2 !== null && r2 >= TRUST_R2;

      obsRow.length = 0;
      obsCol.length = 0;
      obsVal.length = 0;
      riBuf.length = 0;
      ciBuf.length = 0;
      valBuf.length = 0;
      riAlt.length = 0;
      ciAlt.length = 0;
      valAlt.length = 0;
      mols.clear();
      let realCells = 0;
      let impossibleCells = 0;
      let trusted = 0;
      let lo = Infinity;
      let hi = -Infinity;
      let best: {ri: number, ci: number, value: number} | null = null;
      let bestVirtual: {ri: number, ci: number, value: number} | null = null;
      let virtualCells = 0;

      for (let ri = 0; ri < nRows; ri++) {
        rowCells.length = 0;
        for (let ci = 0; ci < nCols; ci++) {
          const cell = matrix.cells[ri][ci];
          const value = cell.value;
          if (cell.kind === 'real' && value !== null && cell.molIdx !== null) {
            compounds.add(cell.molIdx);
            mols.add(cell.molIdx);
            measuredCells++;
            realCells++;
            obsRow.push(ri);
            obsCol.push(ci);
            obsVal.push(value);
            rowCells.push({ci, value, molIdx: cell.molIdx});
            if (roleNames.length > 0) {
              for (let r = 0; r < roleNames.length; r++) {
                const level = roleKind[r] === 0 ? matrix.rows[ri].coreSmiles :
                  roleKind[r] === 1 ? matrix.columns[ci].substSmiles :
                    matrix.rows[ri].foldedValues[roleNames[r]] ?? '';
                roleValues[r].push(level);
                const held = roleBest[r].get(level);
                if (held === undefined || dir * value > dir * held.value)
                  roleBest[r].set(level, {matrix, ri, ci, value});
              }
              roleActivity.push(value);
              roleMol.push(cell.molIdx);
            }
            if (value < lo)
              lo = value;
            if (value > hi)
              hi = value;
            if (value < minObserved)
              minObserved = value;
            if (value > maxObserved)
              maxObserved = value;
            if (best === null || dir * value > dir * best.value)
              best = {ri, ci, value};
            if (cell.smiles !== null)
              measuredStructures.add(cell.smiles);
          } else if (cell.kind === 'unmeasured' && cell.molIdx !== null) {
            untested.add(cell.molIdx);
            if (cell.smiles !== null)
              measuredStructures.add(cell.smiles);
            if (value !== null)
              test.offer(String(cell.molIdx), dir * value, matrix, ri, ci);
          } else if (cell.kind === 'virtual' && value !== null) {
            const support = cell.support ?? 0;
            virtualCells++;
            const passes = support >= MIN_SUPPORT && holds;
            if (passes && (bestVirtual === null || dir * value > dir * bestVirtual.value))
              bestVirtual = {ri, ci, value};
            // A structureless prediction is neither actionable nor withheld for want of evidence: the
            // core carries an attachment point no picked column fills, which is a settings fix.
            if (cell.smiles === null) {
              if (passes)
                trustedNoStructure++;
            } else if (passes) {
              trustedCells++;
              trusted++;
              trustedStructures.add(cell.smiles);
              if (measuredStructures.has(cell.smiles))
                alreadyHeld.add(cell.smiles);
              else {
                riBuf.push(ri);
                ciBuf.push(ci);
                valBuf.push(value);
              }
            } else {
              if (!measuredStructures.has(cell.smiles)) {
                riAlt.push(ri);
                ciAlt.push(ci);
                valAlt.push(value);
              }
              if (support < MIN_SUPPORT)
                withheldThinSupport++;
              else if (r2 === null)
                withheldUnchecked++;
              else
                withheldFitFails++;
            }
          } else if (cell.kind === 'impossible')
            impossibleCells++;
        }
        if (rowCells.length >= 2) {
          swapCandidates++;
          if (rowCells.length > SWAP_ROW_CAP)
            sampledRows++;
          this.poolSwaps(swaps, matrix, ri, rowCells, dir, log, roots[mi]);
        }
      }

      const fit = fitAdditiveFromTriples(obsRow, obsCol, obsVal, nRows, nCols);
      if (!fit.converged)
        nonConverged++;
      // Shrinkage toward the zero initialisation differs per design, so a gain read off a fit that
      // stopped short is not comparable with one read off a converged fit.
      const conf = matrix.confidence ?? null;
      if (!fit.converged || best === null || conf === null)
        withheldNotConverged += riBuf.length;
      else {
        // A fit whose error rests on a handful of cross-validatable cells is not a scale to divide by:
        // a plane through four points has a tiny RMSE and would top every ranking.
        const thin = conf.n < BEST_FIT_MIN_N || conf.rmse <= 0;
        for (let k = 0; k < riBuf.length; k++) {
          const gain = dir * (valBuf[k] - best.value);
          if (thin) {
            if (gain > 0) {
              analogsThin.offer(matrix.cells[riBuf[k]][ciBuf[k]].smiles!, gain,
                matrix, riBuf[k], ciBuf[k]);
            }
          } else if (gain < conf.rmse)
            withheldBelowError++;
          else {
            analogs.offer(matrix.cells[riBuf[k]][ciBuf[k]].smiles!, gain / conf.rmse,
              matrix, riBuf[k], ciBuf[k]);
          }
        }
      }
      // Ranked on raw gain, never on gain over error: the fits behind these are the ones that were not
      // trusted, so their error is not a scale to divide by.
      if (best !== null) {
        for (let k = 0; k < riAlt.length; k++) {
          const gain = dir * (valAlt[k] - best.value);
          if (gain > 0) {
            analogsAny.offer(matrix.cells[riAlt[k]][ciAlt[k]].smiles!, gain,
              matrix, riAlt[k], ciAlt[k]);
          }
        }
      }
      const stat: SeriesStat = {
        matrix, root: roots[mi], tier: tiers[mi], cpd: mols.size, realCells, impossibleCells,
        totalCells: nRows * nCols, lo: realCells ? lo : 0, hi: realCells ? hi : 0, best, bestVirtual,
        typical: fit.grandMean, trusted, converged: fit.converged, virtualCells,
        stripCols: this.recordRGroupExtremes(acc, matrix, fit, dir, roots[mi], holds, wantStrip),
        ...this.extractRanges(matrix, fit, dir),
      };
      series.push(stat);

      if (r2 === null)
        unchecked++;
      else if (r2 < TRUST_R2) {
        lowR2.push(matrix);
        lowR2Virtual += matrix.virtualCount;
      } else
        fitHolds++;
    }

    if (nonConverged > 0) {
      _package.logger.warning(`SAR Matrix | the additive fit of ${nonConverged} series did not reach ` +
        'its tolerance; they are left out of the pooled R-group comparison');
    }
    if (sampledRows > 0) {
      _package.logger.warning(`SAR Matrix | ${sampledRows} rows carry more than ${SWAP_ROW_CAP} measured ` +
        'cells; the swap pool keeps their most and least potent halves, so mid-range pairs are dropped');
    }

    // The walk offers a candidate against the structures seen so far, so a molecule first measured in a
    // later matrix is only caught here.
    const dropHeld = (key: string): boolean => {
      if (!measuredStructures.has(key))
        return false;
      alreadyHeld.add(key);
      return true;
    };
    analogs.prune(dropHeld);
    analogsThin.prune(dropHeld);
    analogsAny.prune(dropHeld);

    // Over the matrices this walk read, not over every matrix: at a tier filter none of the others
    // reaches any other number on the tab, and a model error pooled across tiers is the error of no
    // ranking shown — it sets the band every effect is read against and the Start-here threshold.
    const walked = host.matrices.filter((_, mi) =>
      this.tierFilter === null || tiers[mi] === this.tierFilter);
    const confN = walked.map((m) => m.confidence?.n ?? null)
      .filter((n): n is number => n !== null);
    const rmses = walked.map((m) => m.confidence?.rmse ?? null)
      .filter((r): r is number => r !== null);
    const families = new Set(roots).size;
    const roleFit = fitRoleEffects({names: roleNames, values: roleValues, activity: roleActivity,
      molIdx: roleMol, minSupport: MIN_SUPPORT, higherIsBetter: host.higherIsBetter});
    const data: SummaryData = {
      compounds: compounds.size,
      untested: untested.size,
      measuredCells,
      minObserved: measuredCells ? minObserved : null,
      maxObserved: measuredCells ? maxObserved : null,
      trustedCells,
      trustedStructures: trustedStructures.size,
      trustedNoStructure,
      modelError: rmses.length ? median(rmses) : null,
      fitHolds, unchecked, lowR2, lowR2Virtual,
      families,
      axisRole: host.axisRole,
      coreRole: host.coreRole,
      coresAreSeries: host.coresAreSeries,
      roleFit,
      roleBest: new Map(roleNames.map((name, r) => [name, roleBest[r]])),
      tierCounts: this.countTiers(tiers),
      nonConverged,
      series,
      startHere: this.mergeStartHere(series, confN, dir),
      swaps: this.rankSwaps(swaps),
      swapsByRole: this.poolRoleSwaps(roleNames, roleValues, roleActivity, roleMol, dir, log),
      swapCandidates,
      swapBest: this.largestSwap(swaps),
      ...this.rankRGroups(acc),
      test, analogs, analogsThin, analogsAny,
      withheldBelowError, withheldThinSupport, withheldFitFails, withheldUnchecked, withheldNotConverged,
      alreadyHeld: alreadyHeld.size,
      analogOverflow: analogs.overflow + analogsThin.overflow,
    };
    logSarTime('summary', t0);
    return data;
  }

  private countTiers(tiers: number[]): {tier: number, n: number}[] {
    const byTier = new Map<number, number>();
    for (const tier of tiers)
      byTier.set(tier, (byTier.get(tier) ?? 0) + 1);
    return [...byTier.entries()].sort((a, b) => a[0] - b[0]).map(([tier, n]) => ({tier, n}));
  }

  /** The scalars the screen needs out of one fit, so the Float64Arrays can be dropped with it. Both
   *  ranges are of centred effects, so they are comparable to each other inside this matrix — and to
   *  nothing outside it. One support floor on both axes: an effect estimated from two observations is
   *  noisier than one estimated from three, so admitting columns at a lower floor than rows would
   *  widen the column range before any chemistry entered. */
  private extractRanges(matrix: SarMatrix, fit: ReturnType<typeof fitAdditiveFromTriples>, dir: number):
    {colRange: number | null, rowRange: number | null, bestRow: {ri: number, effect: number, n: number} | null} {
    let colLo = Infinity;
    let colHi = -Infinity;
    let colN = 0;
    for (let ci = 0; ci < matrix.columns.length; ci++) {
      if (fit.colN[ci] < MIN_SUPPORT)
        continue;
      colN++;
      colLo = Math.min(colLo, fit.colEffect[ci]);
      colHi = Math.max(colHi, fit.colEffect[ci]);
    }
    let rowLo = Infinity;
    let rowHi = -Infinity;
    let rowN = 0;
    let bestRow: {ri: number, effect: number, n: number} | null = null;
    for (let ri = 0; ri < matrix.rows.length; ri++) {
      if (fit.rowN[ri] < MIN_SUPPORT)
        continue;
      rowN++;
      rowLo = Math.min(rowLo, fit.rowEffect[ri]);
      rowHi = Math.max(rowHi, fit.rowEffect[ri]);
      if (bestRow === null || dir * fit.rowEffect[ri] > dir * bestRow.effect)
        bestRow = {ri, effect: fit.rowEffect[ri], n: fit.rowN[ri]};
    }
    return {
      colRange: colN >= 2 ? colHi - colLo : null,
      rowRange: rowN >= 2 ? rowHi - rowLo : null,
      bestRow,
    };
  }

  /**
   * Measured pairs inside ONE row index, pooled on the unordered fragment pair.
   *
   * The row index, not the core: a row is keyed on `[coreSmiles, ...foldedValues]`, so two rows can
   * share a core and differ at a folded position — pooling per core would break "everything else
   * identical", which is the whole claim of the card.
   */
  private poolSwaps(pools: Map<string, SwapPool>, matrix: SarMatrix, ri: number, cells: RowCell[],
    dir: number, log: boolean, root: string): void {
    const sampled = cells.length > SWAP_ROW_CAP;
    const kept = keepExtremes(cells, (cell) => cell.value);
    for (let a = 0; a < kept.length; a++) {
      for (let b = a + 1; b < kept.length; b++) {
        const added = this.addSwapPair(pools,
          {value: matrix.columns[kept[a].ci].substSmiles, activity: kept[a].value, mol: kept[a].molIdx},
          {value: matrix.columns[kept[b].ci].substSmiles, activity: kept[b].value, mol: kept[b].molIdx},
          root, sampled, dir, log);
        if (added === null)
          continue;
        const lower = added.flip ? kept[b] : kept[a];
        const upper = added.flip ? kept[a] : kept[b];
        if (added.pool.best === null || Math.abs(added.delta) > Math.abs(added.pool.best.delta))
          added.pool.best = {matrix, ri, ciFrom: lower.ci, ciTo: upper.ci, delta: added.delta};
      }
    }
  }

  /**
   * Fold one measured pair into its pool, or decline it. Null when the two sides carry the same value,
   * when a fold cannot be taken, or when this compound pair is already counted.
   *
   * The caller records `best` rather than this: a pair found through a matrix can name the two cells it
   * came from, and a pair found by grouping the component columns cannot.
   */
  private addSwapPair(pools: Map<string, SwapPool>, a: SwapSide, b: SwapSide,
    root: string, sampled: boolean, dir: number, log: boolean):
    {pool: SwapPool, delta: number, flip: boolean} | null {
    const {value: sa, activity: va, mol: molA} = a;
    const {value: sb, activity: vb, mol: molB} = b;
    if (sa === sb)
      return null;
    const flip = sa > sb;
    const lowerVal = flip ? vb : va;
    const upperVal = flip ? va : vb;
    // A fold multiplies, so its statistics are geometric: pooling the log of the ratio makes the mean,
    // min and max the fold statistics a chemist would quote, and lets one accumulator serve both
    // scales. A non-positive raw value has no fold at all.
    if (!log && (lowerVal <= 0 || upperVal <= 0))
      return null;
    const delta = log ? dir * (upperVal - lowerVal) :
      dir * (Math.log10(upperVal) - Math.log10(lowerVal));
    const key = `${flip ? sb : sa}\0${flip ? sa : sb}`;
    let pool = pools.get(key);
    if (pool === undefined) {
      pool = {from: flip ? sb : sa, to: flip ? sa : sb, n: 0, sum: 0, min: Infinity, max: -Infinity,
        nUp: 0, roots: new Set(), seen: new Set(), best: null, sampled: false};
      pools.set(key, pool);
    }
    const pair = `${Math.min(molA, molB)}:${Math.max(molA, molB)}`;
    if (pool.seen.has(pair))
      return null;
    pool.seen.add(pair);
    pool.n++;
    pool.sum += delta;
    pool.min = Math.min(pool.min, delta);
    pool.max = Math.max(pool.max, delta);
    if (delta > 0)
      pool.nUp++;
    pool.roots.add(root);
    pool.sampled = pool.sampled || sampled;
    return {pool, delta, flip};
  }

  /**
   * Measured swaps for every component at once, read off the component values of each measured cell
   * rather than off one matrix's columns.
   *
   * Nothing is rebuilt to reach another component. A swap is two compounds alike in every component but
   * one; where the components are given as columns that is a grouping, not a decomposition — group the
   * measured cells on every component except the one being swapped, and each pair inside a group that
   * differs in it is a matched pair. The matrix columns were only ever one route to the same thing, and
   * being one route is what made the other components cost a rebuild.
   */
  private poolRoleSwaps(roleNames: string[], roleValues: string[][], roleActivity: number[],
    roleMol: number[], dir: number, log: boolean): Map<string, SwapPool[]> {
    const out = new Map<string, SwapPool[]>();
    for (let r = 0; r < roleNames.length; r++) {
      const buckets = new Map<string, number[]>();
      // One compound occupies a cell in every tier that folded it, so the same partner tuple arrives
      // several times; the first occurrence can form every pair the later ones could.
      const placed = new Map<string, Set<number>>();
      for (let k = 0; k < roleActivity.length; k++) {
        let key = '';
        for (let j = 0; j < roleNames.length; j++) {
          if (j !== r)
            key += `${roleValues[j][k]}\u0001`;
        }
        let bucket = buckets.get(key);
        if (bucket === undefined) {
          bucket = [];
          buckets.set(key, bucket);
          placed.set(key, new Set());
        }
        const seen = placed.get(key)!;
        if (seen.has(roleMol[k]))
          continue;
        seen.add(roleMol[k]);
        bucket.push(k);
      }
      const pools = new Map<string, SwapPool>();
      for (const [key, bucket] of buckets) {
        const sampled = bucket.length > SWAP_ROW_CAP;
        const kept = keepExtremes(bucket, (i) => roleActivity[i]);
        for (let a = 0; a < kept.length; a++) {
          for (let b = a + 1; b < kept.length; b++) {
            this.addSwapPair(pools,
              {value: roleValues[r][kept[a]], activity: roleActivity[kept[a]], mol: roleMol[kept[a]]},
              {value: roleValues[r][kept[b]], activity: roleActivity[kept[b]], mol: roleMol[kept[b]]},
              key, sampled, dir, log);
          }
        }
      }
      out.set(roleNames[r], this.rankSwaps(pools));
    }
    return out;
  }

  /** The swap's worth in its better direction: what it bought in EVERY pair we have. Ranking a mean
   *  over many small pools selects the noisiest pool instead. */
  private swapScore(pool: SwapPool): number {
    return Math.max(pool.min, -pool.max);
  }

  /** The single widest measured move anywhere, gate or no gate — so a card with no qualifying pool can
   *  still say what the data does hold rather than only what it does not. */
  private largestSwap(pools: Map<string, SwapPool>): SwapPool | null {
    let best: SwapPool | null = null;
    for (const pool of pools.values()) {
      const reach = Math.max(pool.max, -pool.min);
      if (best === null || reach > Math.max(best.max, -best.min))
        best = pool;
    }
    return best;
  }

  private rankSwaps(pools: Map<string, SwapPool>): SwapPool[] {
    return [...pools.values()]
      .filter((pool) => pool.n >= SWAP_MIN_PAIRS && pool.roots.size >= SWAP_MIN_SERIES)
      .sort((a, b) => this.swapScore(b) - this.swapScore(a) ||
        (a.from < b.from ? -1 : a.from > b.from ? 1 : a.to < b.to ? -1 : 1))
      .slice(0, SUM_ROWS);
  }

  /**
   * The substituents the fitted model puts first and last at this matrix's position, recorded once per
   * lineage, plus every substituent this series tried and — where a strip will draw — its columns.
   *
   * The fit, not the raw column mean: a column mean confounds a substituent with whichever cores it
   * happened to be made on, which is the normal case in a modular library. The magnitude is a
   * within-series difference against that series' own reference substituent, so it survives the fact
   * that two matrices centre their effects over different substituent menus. An ordinal survives it
   * too: adding a constant to every `colEffect` of a matrix does not reorder them, which is what lets
   * the last place pool the same way the first does.
   */
  private recordRGroupExtremes(acc: RGroupAcc, matrix: SarMatrix, fit: AdditiveFit, dir: number,
    root: string, fitHolds: boolean, strip: boolean): StripCols | null {
    const reference = matrix.refValues[matrix.positions[0] ?? ''];
    const found = reference ? matrix.columns.findIndex((c) => c.substSmiles === reference) : -1;
    // A reference measured once has a fitted effect made of one residual, so a difference against it
    // is noise wearing a comparator's name.
    const refCi = found >= 0 && fit.colN[found] >= 2 ? found : -1;
    // A series whose fit stopped short contributes to no pool, so it takes no slot either: an empty
    // strip entry would draw the "never tried" mark over a series that tried the group and was dropped.
    const cols: StripCols | null = strip && fit.converged ? new Map() : null;
    let bestCi = -1;
    let worstCi = -1;
    for (let c = 0; c < matrix.columns.length; c++) {
      const subst = matrix.columns[c].substSmiles;
      // Only a converged fit may be compared with another matrix's, and a column measured once has a
      // fitted effect made of one residual. The strip is filled under the same test, so a square and
      // the coverage sentence under it can never describe different sets of columns.
      if (fit.colN[c] < 2 || !fit.converged)
        continue;
      if (cols !== null) {
        // 0 where the column IS the reference — the honest difference against itself, which is not the
        // "no comparator at all" the null stands for.
        cols.set(subst, {ci: c, n: fit.colN[c],
          refDelta: refCi < 0 ? null : refCi === c ? 0 :
            dir * (fit.colEffect[c] - fit.colEffect[refCi])});
      }
      let attempts = acc.tried.get(subst);
      if (attempts === undefined)
        acc.tried.set(subst, attempts = new Set());
      attempts.add(root);
      if (bestCi < 0 || dir * fit.colEffect[c] > dir * fit.colEffect[bestCi] ||
        (fit.colEffect[c] === fit.colEffect[bestCi] && subst < matrix.columns[bestCi].substSmiles))
        bestCi = c;
      if (worstCi < 0 || dir * fit.colEffect[c] < dir * fit.colEffect[worstCi] ||
        (fit.colEffect[c] === fit.colEffect[worstCi] && subst < matrix.columns[worstCi].substSmiles))
        worstCi = c;
    }
    if (bestCi >= 0)
      this.pushExtreme(acc.wins, matrix, fit, dir, root, fitHolds, bestCi, refCi, 1);
    // One comparable column is first and last at once, which is not a loss.
    if (worstCi >= 0 && worstCi !== bestCi)
      this.pushExtreme(acc.losses, matrix, fit, dir, root, fitHolds, worstCi, refCi, -1);
    return cols;
  }

  /** Record one column as this lineage's extreme, keeping the occurrence with the stronger
   *  within-series margin. `keep` is +1 for the winner pool and −1 for the loser pool, so one dedup
   *  serves both. */
  private pushExtreme(pool: Map<string, RGroupWin[]>, matrix: SarMatrix, fit: AdditiveFit, dir: number,
    root: string, fitHolds: boolean, ci: number, refCi: number, keep: number): void {
    // The most potent measured cell of the column, so the row lands on a compound rather than on a
    // hole the reader has to hunt through.
    const ri = bestMeasuredRow(matrix, ci, dir);
    if (ri < 0)
      return;
    const win: RGroupWin = {
      matrix, ri, ci, root, n: fit.colN[ci],
      refDelta: refCi >= 0 && refCi !== ci ? dir * (fit.colEffect[ci] - fit.colEffect[refCi]) : null,
      fitHolds,
    };
    const subst = matrix.columns[ci].substSmiles;
    const held = pool.get(subst);
    if (held === undefined) {
      pool.set(subst, [win]);
      return;
    }
    // One entry per lineage: a substituent extreme in a matrix and again in the tier that folded it has
    // been confirmed once, over one set of compounds, and must not read as two series.
    const at = held.findIndex((other) => other.root === root);
    if (at < 0)
      held.push(win);
    else if (this.strongerWin(win, held[at], keep))
      held[at] = win;
  }

  /**
   * Orders two occurrences of one substituent, best of `keep`'s pool first.
   *
   * Never on the fitted level itself: `colEffect` is centred over whichever substituents its own
   * matrix happened to make, so two matrices' levels have different zeros and picking a maximum over
   * them is arithmetic on two different gauges. The difference against each matrix's own reference has
   * no gauge in it; where a matrix records no reference, the measured cells behind the column decide.
   */
  private strongerWin(a: RGroupWin, b: RGroupWin, keep: number): boolean {
    if (a.refDelta !== null && b.refDelta !== null && a.refDelta !== b.refDelta)
      return keep * a.refDelta > keep * b.refDelta;
    if ((a.refDelta === null) !== (b.refDelta === null))
      return b.refDelta === null;
    return a.n !== b.n ? a.n > b.n : a.matrix.id < b.matrix.id;
  }

  private rankRGroups(acc: RGroupAcc):
    {rgroups: RGroupRow[], rgroupsThin: RGroupRow[], rgroupLosers: RGroupRow[]} {
    const winners = this.rankPool(acc.wins, acc.tried, 1);
    const losers = this.rankPool(acc.losses, acc.tried, -1);
    return {
      rgroups: winners.filter((row) => row.k >= RGROUP_MIN_SERIES).slice(0, SUM_ROWS + 2),
      rgroupsThin: winners.filter((row) => row.k < RGROUP_MIN_SERIES).slice(0, SUM_ROWS),
      // No thin tier for losers: a loss over two lineages is not a rule, and striking a building block
      // from the next enumeration is acted on without being checked.
      rgroupLosers: losers.filter((row) => row.k >= RGROUP_MIN_SERIES).slice(0, LOSER_ROWS),
    };
  }

  /** One pool ordered best-first, `keep` being +1 for wins and −1 for losses so the strongest of each
   *  pool comes out on top of its own list. */
  private rankPool(pool: Map<string, RGroupWin[]>, tried: Map<string, Set<string>>,
    keep: number): RGroupRow[] {
    const rows: RGroupRow[] = [];
    for (const [subst, entries] of pool) {
      const trusted = entries.filter((win) => win.fitHolds);
      if (trusted.length < SWAP_MIN_SERIES)
        continue;
      const deltas = trusted.map((win) => win.refDelta).filter((d): d is number => d !== null);
      const ranked = [...trusted].sort((a, b) => this.strongerWin(a, b, keep) ? -1 : 1);
      rows.push({
        subst, win: ranked[0], k: trusted.length, m: entries.length,
        tried: tried.get(subst)?.size ?? entries.length,
        magnitude: deltas.length ? median(deltas) : null,
        lo: deltas.length ? Math.min(...deltas) : 0,
        hi: deltas.length ? Math.max(...deltas) : 0,
        magnitudeSeries: deltas.length,
      });
    }
    // Rows with no comparator sort last in either pool: they still rank on how often the group placed.
    const strength = (row: RGroupRow): number =>
      row.magnitude === null ? Number.NEGATIVE_INFINITY : keep * row.magnitude;
    rows.sort((a, b) => b.k - a.k || strength(b) - strength(a) || (a.subst < b.subst ? -1 : 1));
    return rows;
  }

  /**
   * Five independent argmaxes, then merged by lineage root.
   *
   * Independent because a series holding the most compounds, the widest range AND the best compound is
   * the answer, and pushing it out of a taken slot would put a weaker series under a label it does not
   * earn. Merged because those are one place to look, not three.
   */
  private mergeStartHere(series: SeriesStat[], confN: number[], dir: number): StartRow[] {
    if (series.length === 0)
      return [];
    // Floored well above MIN_CV_POINTS: an R² ranked at the cross-validation floor is an extreme-value
    // draw from the noisiest end, on a row captioned "best validated".
    const nFloor = Math.max(BEST_FIT_MIN_N, median(confN));
    const pick = (score: (s: SeriesStat) => number | null): SeriesStat | null => {
      let winner: SeriesStat | null = null;
      let bestScore = 0;
      for (const stat of series) {
        const value = score(stat);
        if (value === null)
          continue;
        if (winner === null || value > bestScore ||
          (value === bestScore && stat.matrix.id < winner.matrix.id)) {
          winner = stat;
          bestScore = value;
        }
      }
      return winner;
    };
    // Index-aligned with REASON_GLYPHS / REASON_WORDS: the lane is a fixed set of slots, so the order
    // here is the order on screen.
    const slots: (SeriesStat | null)[] = [
      pick((s) => s.cpd || null),
      pick((s) => s.realCells ? s.hi - s.lo : null),
      pick((s) => s.best === null ? null : dir * s.best.value),
      pick((s) => s.trusted || null),
      pick((s) => s.matrix.confidence && s.matrix.confidence.n >= nFloor ? s.matrix.confidence.r2 : null),
    ];
    const byRoot = new Map<string, StartRow>();
    for (let slot = 0; slot < slots.length; slot++) {
      const stat = slots[slot];
      if (stat === null)
        continue;
      const held = byRoot.get(stat.root);
      if (held === undefined)
        byRoot.set(stat.root, {primary: stat, reasons: [slot]});
      else
        held.reasons.push(slot);
    }
    return [...byRoot.values()];
  }
}
