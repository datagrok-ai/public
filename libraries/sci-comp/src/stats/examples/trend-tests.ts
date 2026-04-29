/**
 * Example: testing for an ordered (dose-response) trend across groups.
 *
 *   - `jonckheere` is a non-parametric trend test for *continuous* data
 *     across ordered groups. It is built from pairwise Mann-Whitney U
 *     counts and standardised using the moments under H₀.
 *   - `cochranArmitage` is the analogous test for *binary* (incidence)
 *     data: events / total per dose group with chosen scores.
 *
 * Both reduce to a single Z-statistic and two-sided p-value.
 *
 * References:
 *   - Jonckheere AR (1954). "A distribution-free k-sample test against
 *     ordered alternatives." Biometrika 41, 133-145.
 *   - Terpstra TJ (1952). "The asymptotic normality and consistency of
 *     Kendall's test against trend, when ties are present in one ranking."
 *     Indagationes Mathematicae 14, 327-333.
 *   - Cochran WG (1954). "Some methods for strengthening the common
 *     chi-squared tests." Biometrics 10, 417-451.
 *   - Armitage P (1955). "Tests for linear trends in proportions and
 *     frequencies." Biometrics 11, 375-386.
 *   - Tang ML et al. (2006). "Exact Cochran-Armitage trend tests."
 *     J Stat Comput Simul 76(10), 847-859, Table 1 (rat neurotoxicity).
 *
 * Run: npx tsx src/stats/examples/trend-tests.ts
 */

import {cochranArmitage, jonckheere} from '..';
import {fmt, printCase, printKv, printProblem, printSection} from './_helpers';

/* ================================================================== */
/*  Continuous trend — Jonckheere-Terpstra                            */
/* ================================================================== */

printSection('JONCKHEERE-TERPSTRA (CONTINUOUS DATA)');

// Three dose groups with monotonically increasing values: a clear trend.
const G1 = [4.2, 5.1, 5.5, 5.8, 6.0];
const G2 = [5.6, 6.0, 6.3, 6.7, 7.0];
const G3 = [6.4, 7.1, 7.5, 7.9, 8.2];

printProblem([
  'Three dose groups, n = 5 each, with monotonically increasing values.',
  'goal: detect ordered alternative G1 ≤ G2 ≤ G3 (with ≥ 1 strict).',
  '',
]);

printCase('Three ordered groups');
const jt1 = jonckheere([G1, G2, G3]);
printKv('Z-statistic', fmt(jt1.statistic));
printKv('p-value (two-sided)', fmt(jt1.pValue, 6));

printProblem([
  '',
  'Same three groups in REVERSE order — Z should flip sign, p unchanged.',
  '',
]);

printCase('Reversed group order');
const jt2 = jonckheere([G3, G2, G1]);
printKv('Z-statistic', fmt(jt2.statistic));
printKv('p-value (two-sided)', fmt(jt2.pValue, 6));

/* ================================================================== */
/*  Incidence trend — Cochran-Armitage                                */
/* ================================================================== */

printSection('COCHRAN-ARMITAGE (INCIDENCE / BINARY DATA)');

// Tang et al. (2006) Table 1: 90-day neurotoxicity study, stained-face
// incidence by cyclohexane concentration. 4 dose levels, n=12 each.
const TANG_COUNTS = [0, 0, 1, 3];
const TANG_TOTALS = [12, 12, 12, 12];

printProblem([
  'Tang (2006) Table 1: rat neurotoxicity, 4 dose levels, n = 12 each.',
  'incidence: 0/12, 0/12, 1/12, 3/12',
  '',
]);

printCase('Default (scores = 0..k-1, two-sided, binomial variance)');
const ca1 = cochranArmitage(TANG_COUNTS, TANG_TOTALS);
printKv('Z-statistic', fmt(ca1.zStatistic));
printKv('chi² (= Z²)', fmt(ca1.chi2Statistic));
printKv('p-value', fmt(ca1.pValue, 6));
printKv('p̄ (pooled)', fmt(ca1.pBar));

printCase('One-sided (alternative = "increasing")');
const ca2 = cochranArmitage(TANG_COUNTS, TANG_TOTALS, {alternative: 'increasing'});
printKv('p-value', fmt(ca2.pValue, 6));

printCase('Hypergeometric variance (R prop.trend.test convention)');
const ca3 = cochranArmitage(TANG_COUNTS, TANG_TOTALS, {variance: 'hypergeometric'});
printKv('Z-statistic', fmt(ca3.zStatistic));
printKv('chi²', fmt(ca3.chi2Statistic));

printCase('With Buonaccorsi-style modified statistic');
const ca4 = cochranArmitage(TANG_COUNTS, TANG_TOTALS, {modified: true});
printKv('Z (standard)', fmt(ca4.zStatistic));
printKv('Z (modified)', fmt(ca4.zModified ?? null));
printKv('p (modified)', fmt(ca4.pValueModified ?? null, 6));
