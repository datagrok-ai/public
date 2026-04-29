/**
 * Example: Williams' step-down test for dose-response trend, with PAVA
 * isotonic regression as the inner ingredient.
 *
 * Step 1 — PAVA (Pool-Adjacent-Violators) gives the maximum-likelihood
 *          estimate of the dose-group means under the monotone-trend
 *          constraint.
 *
 * Step 2 — Williams' procedure walks from the highest dose down,
 *          rejecting H₀ at each step until a non-significant result is
 *          met. The lowest dose still rejected is the Minimum Effective
 *          Dose (MED).
 *
 * References:
 *   - Williams DA (1971). "A test for differences between treatment
 *     means when several dose levels are compared with a zero-dose
 *     control." Biometrics 27(1), 103-117.
 *   - Williams DA (1972). "The comparison of several dose levels with
 *     a zero-dose control." Biometrics 28(2), 519-531.
 *   - Barlow RE et al. (1972). Statistical Inference Under Order
 *     Restrictions. Wiley (PAVA reference).
 *
 * Run: npx tsx src/stats/examples/williams-step-down.ts
 */

import {pavaIncreasing, williamsTest} from '..';
import {fmt, printCase, printKv, printProblem, printSection} from './_helpers';

/* ================================================================== */
/*  PAVA demonstration                                                */
/* ================================================================== */

printSection('PAVA — POOL-ADJACENT-VIOLATORS');

// Williams (1971) §2 p.106 numerical example: 7 group means with one
// non-monotone dip at index 1. PAVA pools the violating block.
const PAVA_VALUES = [10.4, 9.9, 10.0, 10.6, 11.4, 11.9, 11.7];
const PAVA_WEIGHTS = [1, 1, 1, 1, 1, 1, 1];

printProblem([
  'Williams (1971) §2 p.106 example: 7 group means.',
  '  raw:     [10.4, 9.9, 10.0, 10.6, 11.4, 11.9, 11.7]',
  '  expected (after PAVA, increasing):',
  '           [10.1, 10.1, 10.1, 10.6, 11.4, 11.8, 11.8]',
  '',
]);

printCase('PAVA increasing');
const pava = pavaIncreasing(PAVA_VALUES, PAVA_WEIGHTS);
console.log(`  result: [${[...pava].map((v) => v.toFixed(2)).join(', ')}]`);

/* ================================================================== */
/*  Williams' step-down test                                          */
/* ================================================================== */

printSection('WILLIAMS\' STEP-DOWN TEST');

// Williams (1971) §5 p.117 worked example: 6 dose groups + control,
// r = 8 per group, pooled within-group variance s² = 1.16.
// Expected MED = group 4 (groups 6, 5, 4 reject; group 3 does not).
const MEANS = [10.4, 9.9, 10.0, 10.6, 11.4, 11.9, 11.7];
const SDS = Array(7).fill(Math.sqrt(1.16));
const NS = [8, 8, 8, 8, 8, 8, 8];
const LABELS = ['ctrl', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6'];

printProblem([
  'Williams (1971) §5 p.117: 6 dose groups, r = 8 each, s² = 1.16.',
  '  index 0 = control, indices 1..6 = increasing doses.',
  '  expected: MED = d4 (doses 6, 5, 4 reject; dose 3 does not).',
  '',
]);

printCase('Williams\' step-down (α = 0.05, increasing)');
const w = williamsTest(MEANS, SDS, NS, LABELS, {direction: 'increase', alpha: 0.05});

printKv('direction', w.direction);
printKv('pooled variance', fmt(w.pooledVariance));
printKv('pooled df', String(w.pooledDf));
printKv('MED', w.minimumEffectiveDose ?? 'none');
printKv('MED index', w.minimumEffectiveIndex !== null ? String(w.minimumEffectiveIndex) : 'none');

console.log('\n  step-down trace (high dose to low):');
console.log('    label  i   mean_PAVA   t̄       cv       sig');
console.log('    -----  -   ---------   ------  ------   ---');
for (const r of w.stepDownResults) {
  const lbl = r.doseLabel.padEnd(5);
  const idx = String(r.doseIndex).padEnd(2);
  const cm = fmt(r.constrainedMean).padEnd(10);
  const t = fmt(r.testStatistic).padEnd(7);
  const cv = fmt(r.criticalValue).padEnd(7);
  const sig = r.significant ? 'YES' : ' no';
  console.log(`    ${lbl}  ${idx}  ${cm}  ${t} ${cv}  ${sig}`);
}

/* ================================================================== */
/*  Auto-direction detection                                          */
/* ================================================================== */

printSection('AUTO-DIRECTION');

printProblem([
  'A decreasing dose-response: highest dose < control.',
  '  means: [50.0, 45.0, 40.0]; ns: [10, 10, 10]; sds: [2, 2, 2]',
  '',
]);

printCase('direction: \'auto\'');
const auto = williamsTest([50.0, 45.0, 40.0], [2, 2, 2], [10, 10, 10],
  ['ctrl', 'd1', 'd2'], {direction: 'auto', alpha: 0.05});
printKv('detected direction', auto.direction);
printKv('MED', auto.minimumEffectiveDose ?? 'none');
