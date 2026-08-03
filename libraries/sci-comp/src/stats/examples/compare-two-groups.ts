/**
 * Example: comparing two groups with Welch's t-test, Mann-Whitney U, and
 * Hedges' g (effect size).
 *
 * The three methods together give a complete picture: significance under
 * Gaussian assumptions (Welch), distribution-free significance (MWU), and
 * the *magnitude* of the difference on a scale-free scale (Hedges' g).
 *
 * References:
 *   - Welch BL (1947). "The generalization of Student's problem when several
 *     different population variances are involved." Biometrika 34, 28-35.
 *   - Mann HB, Whitney DR (1947). "On a test of whether one of two random
 *     variables is stochastically larger than the other." Ann Math Stat
 *     18, 50-60.
 *   - Hedges LV (1981). "Distribution theory for Glass's estimator of effect
 *     size and related estimators." J Educ Stat 6, 107-128.
 *   - Body-temperature dataset: McDonald JH (2014). Handbook of Biological
 *     Statistics, 3rd ed., Sparky House Publishing, pp. 128-129.
 *
 * Run: npx tsx src/stats/examples/compare-two-groups.ts
 */

import {hedgesG, mannWhitneyU, welchTTest} from '..';
import {fmt, printCase, printKv, printProblem, printSection} from './_helpers';

/* ================================================================== */
/*  Datasets                                                          */
/* ================================================================== */

// McDonald (2014) pp. 128-129: human body temperatures (°F) recorded
// at 2pm and 5pm. Different individuals, different sample sizes.
const TEMP_2PM = [69, 70, 66, 63, 68, 70, 69, 67, 62, 63,
  76, 59, 62, 62, 75, 62, 72, 63];
const TEMP_5PM = [68, 62, 67, 68, 69, 67, 61, 59, 62, 61,
  69, 66, 62, 62, 61, 70];

// A clearly separated synthetic example to show what a "significant"
// result looks like.
const CONTROL = [10.2, 11.1, 9.8, 10.5, 10.0];
const TREATED = [15.3, 14.8, 16.1, 15.0, 15.5];

/* ================================================================== */
/*  Demo                                                              */
/* ================================================================== */

printSection('TWO-GROUP COMPARISON');

printProblem([
  'Body temperature 2pm vs 5pm (McDonald 2014, pp. 128-129)',
  'n1 = 18, n2 = 16. R t.test reports t = 1.3109, p = 0.1995 (two-sided).',
  '',
]);

printCase('Welch t-test');
const tw1 = welchTTest(TEMP_2PM, TEMP_5PM);
printKv('t-statistic', fmt(tw1.statistic));
printKv('p-value', fmt(tw1.pValue));

printCase('Mann-Whitney U');
const mw1 = mannWhitneyU(TEMP_2PM, TEMP_5PM, 'two-sided');
printKv('U-statistic', fmt(mw1.statistic, 1));
printKv('p-value', fmt(mw1.pValue));

printCase('Hedges\' g (effect size)');
const g1 = hedgesG(TEMP_2PM, TEMP_5PM);
printKv('g', fmt(g1));
printKv('interpretation', interpret(g1));

printSection('CLEARLY SEPARATED GROUPS');

printProblem([
  'Synthetic control vs treated. Means differ by ~5 SD.',
  '',
]);

printCase('Welch t-test');
const tw2 = welchTTest(CONTROL, TREATED);
printKv('t-statistic', fmt(tw2.statistic));
printKv('p-value', fmt(tw2.pValue, 6));

printCase('Mann-Whitney U');
const mw2 = mannWhitneyU(CONTROL, TREATED, 'two-sided');
printKv('U-statistic', fmt(mw2.statistic, 1));
printKv('p-value', fmt(mw2.pValue, 6));

printCase('Hedges\' g (effect size)');
const g2 = hedgesG(CONTROL, TREATED);
printKv('g', fmt(g2));
printKv('interpretation', interpret(g2));

/**
 * Cohen's rough cutoffs (1988): |g| < 0.2 trivial, 0.2-0.5 small,
 * 0.5-0.8 medium, > 0.8 large.
 */
function interpret(g: number | null): string {
  if (g === null) return 'undefined';
  const a = Math.abs(g);
  if (a < 0.2) return `|g|=${a.toFixed(2)} (trivial)`;
  if (a < 0.5) return `|g|=${a.toFixed(2)} (small)`;
  if (a < 0.8) return `|g|=${a.toFixed(2)} (medium)`;
  return `|g|=${a.toFixed(2)} (large)`;
}
