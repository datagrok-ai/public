/**
 * Example: sequential threshold test for proportions (Young 1985).
 *
 * Walks from the lowest dose, pooling each non-significant group into a
 * cumulative "control". Stops at the first significant comparison —
 * marking that group as the Effect Level (EL) and prior groups as
 * No-Observed-Effect Levels (NOELs).
 *
 * Optional Šidák alpha-adjustment controls the family-wise error rate
 * across the (k - 1) sequential comparisons.
 *
 * References:
 *   - Young SS (1985). "The Cochran-Armitage test for trends or thresholds
 *     in proportions." Society for Risk Analysis, Tables 5(c) and 6.
 *   - Šidák Z (1967). "Rectangular confidence regions for the means of
 *     multivariate normal distributions." J Am Stat Assoc 62, 626-633.
 *
 * Run: npx tsx src/stats/examples/threshold-test.ts
 */

import {thresholdTest} from '..';
import {fmt, printCase, printKv, printProblem, printSection} from './_helpers';

/* ================================================================== */
/*  Dataset                                                           */
/* ================================================================== */

// Young (1985) Table 5(c): liver carcinoma incidence in male B6C3F1 mice
// (NCI study on DDT/DDE/TDE), first 4 dose groups.
// Expected per Young (1985):
//   group 0 vs 1: Z ≈ -1.22, NS
//   groups 0,1 vs 2: Z ≈ -0.79, NS
//   groups 0..2 vs 3: Z ≈ 2.99, significant → EL = group 3
const DDT_COUNTS = [4, 1, 1, 7];
const DDT_TOTALS = [56, 49, 48, 41];

/* ================================================================== */
/*  Demo                                                              */
/* ================================================================== */

printSection('SEQUENTIAL THRESHOLD TEST (YOUNG 1985)');

printProblem([
  'Young (1985) Table 5(c): NCI mice liver carcinoma by DDT dose.',
  '4 dose groups: incidences 4/56, 1/49, 1/48, 7/41.',
  'goal: identify the first dose group whose incidence differs',
  '      significantly from the pooled lower-dose "control".',
  '',
]);

/* ----- Without alpha adjustment ----- */

printCase('Per-comparison α = 0.05 (no Šidák adjustment)');
const stepsRaw = thresholdTest(DDT_COUNTS, DDT_TOTALS, {
  alpha: 0.05,
  adjustAlpha: false,
});
printSteps(stepsRaw);

/* ----- With Šidák adjustment ----- */

printCase('FWER α = 0.05, Šidák-adjusted');
const stepsAdj = thresholdTest(DDT_COUNTS, DDT_TOTALS, {
  alpha: 0.05,
  adjustAlpha: true,
});
printSteps(stepsAdj);

if (stepsAdj.length > 0) {
  printKv('per-comparison α_adj', fmt(stepsAdj[0].alphaAdj, 5));
  const last = stepsAdj[stepsAdj.length - 1];
  printKv('Effect Level group', String(last.effectGroup ?? 'none'));
  printKv('NOEL groups', JSON.stringify(last.noelGroups ?? []));
}

function printSteps(steps: ReturnType<typeof thresholdTest>): void {
  console.log('  step  control     treated     z       p       sig');
  console.log('  ----  ----------  ----------  ------  ------  ---');
  for (let i = 0; i < steps.length; i++) {
    const s = steps[i];
    const idx = String(i + 1).padEnd(4);
    const ctrl = `${s.controlCount}/${s.controlTotal}`.padEnd(10);
    const treat = `${s.treatedCount}/${s.treatedTotal}`.padEnd(10);
    const z = fmt(s.z).padEnd(6);
    const p = fmt(s.p, 4).padEnd(6);
    const sig = s.significant ? 'YES' : ' no';
    console.log(`  ${idx}  ${ctrl}  ${treat}  ${z}  ${p}  ${sig}`);
  }
}
