import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {parseSpectronautText} from '../parsers/spectronaut-parser';
import {getGroups} from '../analysis/experiment-setup';
import {runDifferentialExpression} from '../analysis/differential-expression';
import {_package} from '../package-test';

/**
 * HYE ground-truth regression invariants.
 *
 * The HYE three-proteome spike-in standard (files/demo/spectronaut-hye-mix.tsv,
 * HYE mix A vs B, 4 vs 4) is one of the few proteomics inputs where the answer
 * is KNOWN from organism identity, independent of any method:
 *   - Human  = null      → expected log2FC ≈ 0   (should NOT be called significant)
 *   - Yeast  = 2:1  true → |log2FC| ≈ 1.0        (should be recovered)
 *   - E. coli = 3:1 true → |log2FC| ≈ 1.7        (should be recovered)
 *
 * This drives the ACTUAL report path — parse (PG.IBAQ pivot, log2) then the
 * client-side differential-expression (Welch's t-test + BH FDR, the guaranteed
 * fallback that runs without an R environment) — and asserts three invariants
 * against ground truth. It is a calibration guard: a regression in log2 handling,
 * group assignment, the FDR correction, or precursor→protein aggregation moves
 * these numbers off their known-good operating point.
 *
 * Thresholds are deliberately loose floors/ceilings around the measured operating
 * point (recall 0.900, human FPR 0.033, null median log2FC +0.023, at a permissive
 * |log2FC| >= 0.58 cut) so ordinary numerical drift does not flake, while a real
 * breakage trips them. These reproduce the independent report-path HYE validation.
 */

const FIXTURE_PATH = 'demo/spectronaut-hye-mix.tsv';

// Invariant bounds (measured operating point in parentheses; see file header).
const RECALL_MIN = 0.80; // (0.900) true positives — yeast + E. coli — called significant
const HUMAN_FPR_MAX = 0.10; // (0.033) ≈ 2× nominal 5%; human nulls wrongly called significant
const NULL_MEDFC_MAX = 0.10; // (0.023) |median log2FC| of the human null must sit near 0

/** Ground-truth class from the organism string: null (human) vs true-positive
 *  (yeast / E. coli). Anything else (multi-organism groups) is excluded. */
function truthClass(organism: string): 'null' | 'true' | 'other' {
  const o = organism.toLowerCase();
  if (o.includes('sapiens')) return 'null';
  if (o.includes('cerevisiae') || o.includes('coli')) return 'true';
  return 'other';
}

function median(values: number[]): number {
  if (values.length === 0) return NaN;
  const sorted = [...values].sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
}

category('HYE ground truth', () => {
  test('report path recovers known organism truth (recall / null FPR / null bias)', async () => {
    const text = await _package.files.readAsText(FIXTURE_PATH);
    const df = await parseSpectronautText(text);

    // The mix has exactly two conditions, so the parser auto-populates groups.
    const groups = getGroups(df);
    expect(groups !== null, true);

    // Probe at a permissive |log2FC| >= 0.58 cut (p <= 0.05). The tool's DEFAULT
    // cut is 1.0, but yeast's true |log2FC| ~= 1.0 sits exactly on it, so
    // default-threshold recall straddles the boundary and is intrinsically noisy.
    // 0.58 clears both true classes (yeast 1.0, E. coli 1.7) so the recall floor
    // reflects the method's discriminative power, not a threshold knife-edge.
    const FC_CUT = 0.58; const P_CUT = 0.05;
    runDifferentialExpression(df, groups!.group1.columns, groups!.group2.columns,
      groups!.group1.name, groups!.group2.name, FC_CUT, P_CUT);

    const orgCol = df.col('PG.Organisms') ??
      df.columns.toList().find((c) => c.name.toLowerCase().includes('organism')) ?? null;
    expect(orgCol !== null, true);
    const sigCol = df.col('significant')!;
    const fcCol = df.col('log2FC')!;

    let tp = 0; let truePositives = 0; // recovered / total true-difference proteins
    let fp = 0; let nulls = 0; // human nulls wrongly called / total human nulls
    const nullFcs: number[] = [];

    for (let i = 0; i < df.rowCount; i++) {
      const org = orgCol!.get(i);
      if (typeof org !== 'string') continue;
      const cls = truthClass(org);
      const sig = sigCol.get(i) === true;
      const fc = fcCol.get(i);
      if (cls === 'true') {
        truePositives++;
        if (sig) tp++;
      } else if (cls === 'null') {
        nulls++;
        if (sig) fp++;
        if (typeof fc === 'number' && Number.isFinite(fc)) nullFcs.push(fc);
      }
    }

    // Sanity: the three-proteome standard must actually contain all three classes.
    expect(truePositives > 40, true);
    expect(nulls > 20, true);

    const recall = tp / truePositives;
    const humanFpr = fp / nulls;
    const nullMedFc = Math.abs(median(nullFcs));

    // Invariant 1 — recall: the tool recovers the bulk of the true differences.
    expect(recall >= RECALL_MIN, true,
      `recall ${recall.toFixed(3)} below ${RECALL_MIN} (${tp}/${truePositives})`);
    // Invariant 2 — null false-positive rate stays near nominal.
    expect(humanFpr <= HUMAN_FPR_MAX, true,
      `human FPR ${humanFpr.toFixed(3)} above ${HUMAN_FPR_MAX} (${fp}/${nulls})`);
    // Invariant 3 — the null is unbiased: its median fold change sits at ~0.
    expect(nullMedFc <= NULL_MEDFC_MAX, true,
      `|null median log2FC| ${nullMedFc.toFixed(3)} above ${NULL_MEDFC_MAX}`);
  });
});
