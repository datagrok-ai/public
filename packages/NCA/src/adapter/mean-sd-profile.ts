/**
 * Mean ± SD concentration-time profile builder (VR-03 spike).
 *
 * Pure `DataFrame → DataFrame` aggregation: pools an observations frame into a
 * mean ± sample-SD profile per (group × time), the data primitive behind the
 * "mean ± SD by dose group" overlay (PRD VR-03). For SERIAL designs this is the
 * group analogue of the per-subject overlay; for SPARSE/destructive designs the
 * sparse path (`sparse-nca.ts`) already computes the same mean/SD per nominal
 * time — this helper covers the serial case without duplicating that math.
 *
 * Emits one row per (group × time) with TWO uncertainty representations:
 *
 *  - **Arithmetic** (for a LINEAR axis): `Mean`, `SD`, `N`, `Upper` = Mean+SD,
 *    `Lower` = max(Mean−SD, 0). Computed over all finite `c`, with BLQ-as-0.
 *  - **Geometric** (for a SEMI-LOG axis): `geoMean` = exp(mean(ln c)),
 *    `geoUpper` = geoMean·gSD, `geoLower` = geoMean/gSD where gSD = exp(sd(ln c))
 *    is the geometric SD factor, plus `nPos` (count of positive observations).
 *    Computed over `c > 0` ONLY.
 *
 * **Why geometric for log.** Arithmetic mean−SD goes ≤ 0 at high CV / late
 * timepoints, which is undefined on a log axis (the band would clip at the axis
 * floor). The geometric mean ×/÷ gSD is the standard PK measure for a semi-log
 * profile: always positive, symmetric in log space, defined only over measurable
 * (c>0) samples. The chart binds `Mean`/`Lower`/`Upper` on the linear axis and
 * `geoMean`/`geoLower`/`geoUpper` on the semi-log axis.
 *
 * **BLQ convention differs by representation — caller must surface it.** The
 * arithmetic mean includes BLQ as 0 (pulls the mean DOWN); the geometric mean
 * excludes BLQ entirely (log undefined for ≤0). So at a partial-BLQ timepoint the
 * geometric mean is computed from the positives only and is biased HIGH — the
 * magnitude scales with the BLQ fraction and can be substantial (>10%) once fewer
 * than half the group is above LLOQ (routine at late timepoints). The two axes can
 * legitimately disagree on central tendency there. Both `nPos`
 * (positives) and `N` (total) ride the output so the consumer can disclose the
 * exclusion (e.g. tooltip) when `nPos < N` — required for honest uncertainty
 * (CLAUDE.md §4). No LLOQ substitution (LLOQ/2 etc.) is applied here: BLQ policy
 * lives upstream (the observations frame already reflects it), not in this summary.
 *
 * Sample SD (n−1) for both. `N`/`nPos = 1` → SD undefined: the arithmetic band
 * collapses (SD=0) and the geometric BAND is FLOAT_NULL (`geoMean` kept; a
 * collapsed geo band would read as false zero-uncertainty). All-BLQ cells
 * (nPos = 0) → all geometric columns FLOAT_NULL (no positive datum; arithmetic
 * Mean still 0). Both variances use Welford's online algorithm (numerically
 * stable — no catastrophic cancellation).
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

import {f64Column} from './adapter-utils';

export interface MeanSdProfileSpec {
  /** Time column (observed or nominal). */
  readonly timeCol: string;
  /** Concentration column. */
  readonly concCol: string;
  /** Grouping column (dose group / sex / day / compound) — one curve per value. */
  readonly groupCol: string;
}

/** One pooled (group × time) cell. Variance accumulated via Welford's online
 *  algorithm (numerically stable — no `sumSq − n·mean²` catastrophic
 *  cancellation for near-equal concentrations). */
interface Cell {
  readonly group: string;
  readonly time: number;
  mean: number;
  /** Sum of squared deviations from the running mean (Welford M2). */
  m2: number;
  n: number;
  /** Geometric accumulators: Welford over ln(c) for c > 0 only. */
  glnMean: number;
  glnM2: number;
  nPos: number;
}

/**
 * Pool `obs` into a mean ± SD profile per (group × time). Rows with a non-finite
 * time or concentration are skipped (they cannot anchor a composite point).
 */
export function buildMeanSdProfile(obs: DG.DataFrame, spec: MeanSdProfileSpec): DG.DataFrame {
  const timeCol = mustCol(obs, spec.timeCol);
  const concCol = mustCol(obs, spec.concCol);
  const groupCol = mustCol(obs, spec.groupCol);

  const cells = new Map<string, Cell>();
  const order: string[] = [];
  for (let i = 0; i < obs.rowCount; i++) {
    if (concCol.isNone(i) || timeCol.isNone(i) || groupCol.isNone(i)) continue;
    const t = timeCol.get(i) as number;
    const c = concCol.get(i) as number;
    if (!Number.isFinite(t) || !Number.isFinite(c)) continue;
    const g = String(groupCol.get(i));
    const key = `${g}|${t}`;
    let cell = cells.get(key);
    if (cell === undefined) {
      cell = {group: g, time: t, mean: 0, m2: 0, n: 0, glnMean: 0, glnM2: 0, nPos: 0};
      cells.set(key, cell);
      order.push(key);
    }
    // Arithmetic Welford — stable running mean + sum-of-squared-deviations (all
    // finite c, including BLQ-as-0).
    cell.n++;
    const delta = c - cell.mean;
    cell.mean += delta / cell.n;
    cell.m2 += delta * (c - cell.mean);
    // Geometric Welford — over ln(c) for measurable (c > 0) samples only. BLQ/zero
    // contribute nothing to the geometric mean (log undefined for ≤ 0).
    if (c > 0) {
      cell.nPos++;
      const lc = Math.log(c);
      const gd = lc - cell.glnMean;
      cell.glnMean += gd / cell.nPos;
      cell.glnM2 += gd * (lc - cell.glnMean);
    }
  }

  const group: string[] = [];
  const time: number[] = []; const mean: number[] = []; const sd: number[] = [];
  const n: number[] = []; const upper: number[] = []; const lower: number[] = [];
  const geoMean: number[] = []; const geoUpper: number[] = []; const geoLower: number[] = [];
  const nPos: number[] = [];
  for (const key of order) {
    const cell = cells.get(key)!;
    const m = cell.mean;
    // Arithmetic sample (n−1) SD from Welford M2; single-observation cell → 0
    // spread (flagged by N=1). M2 is non-negative by construction.
    const variance = cell.n > 1 ? cell.m2 / (cell.n - 1) : 0;
    const s = Math.sqrt(variance);
    const up = m + s; const lo = Math.max(0, m - s);
    group.push(cell.group); time.push(cell.time); mean.push(m); sd.push(s); n.push(cell.n);
    upper.push(up); lower.push(lo);
    // Geometric: exp(mean ± sd of ln c) over c>0.
    //  - nPos=0 → no positive datum → all FLOAT_NULL (omitted on the log axis).
    //  - nPos=1 → SD undefined → mean kept but the BAND is FLOAT_NULL (no estimable
    //    spread; a collapsed band would read as false zero-uncertainty — WinNonlin/
    //    SAS convention is a null bound at n=1).
    //  - nPos<N (some BLQ) → geoMean is over positives only and biases slightly high
    //    at partial-BLQ timepoints; nPos + N both ride the frame so the consumer can
    //    surface the exclusion (CLAUDE.md §4 honest uncertainty).
    nPos.push(cell.nPos);
    if (cell.nPos === 0) {
      geoMean.push(DGNS.FLOAT_NULL); geoUpper.push(DGNS.FLOAT_NULL); geoLower.push(DGNS.FLOAT_NULL);
    } else if (cell.nPos === 1) {
      geoMean.push(Math.exp(cell.glnMean)); geoUpper.push(DGNS.FLOAT_NULL); geoLower.push(DGNS.FLOAT_NULL);
    } else {
      const glnSd = Math.sqrt(cell.glnM2 / (cell.nPos - 1));
      geoMean.push(Math.exp(cell.glnMean));
      geoUpper.push(Math.exp(cell.glnMean + glnSd));
      geoLower.push(Math.exp(cell.glnMean - glnSd));
    }
  }

  const df = DGNS.DataFrame.fromColumns([
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'Group', group),
    f64Column('Time', time),
    f64Column('Mean', mean),
    f64Column('SD', sd),
    DGNS.Column.fromInt32Array('N', Int32Array.from(n)),
    f64Column('Upper', upper),
    f64Column('Lower', lower),
    f64Column('geoMean', geoMean),
    f64Column('geoUpper', geoUpper),
    f64Column('geoLower', geoLower),
    DGNS.Column.fromInt32Array('nPos', Int32Array.from(nPos)),
  ]);

  // Freeze a NATURAL category order on `Group` so a consuming viewer's colour
  // legend reads dose labels in magnitude order (10, 30, 100, 300), not the
  // string-lexical default (10, 100, 30, 300). This is the platform's own
  // "Sort → Natural" mechanism: `setCategoryOrder` persists the `.category-order`
  // column tag (DG.Tags.CategoryOrder) that viewers honour — not a chart-side hack.
  // Presentation metadata only: no pooled statistic or plotted value changes.
  // Natural (numeric-aware) compare degrades to plain lexical for non-numeric
  // groups (sex F/M, compound IDs, route IV/PO). Clones for the semi-log frame
  // inherit the tag, so both axes are covered.
  const groupOut = df.col('Group')!;
  groupOut.setCategoryOrder(
    groupOut.categories.slice().sort((a, b) => a.localeCompare(b, undefined, {numeric: true})));
  return df;
}

function mustCol(df: DG.DataFrame, name: string): DG.Column {
  const col = df.col(name);
  if (col === null) throw new Error(`buildMeanSdProfile: column '${name}' not found.`);
  return col;
}
