/**
 * Example: Boschloo's unconditional exact test for a 2×2 contingency table.
 *
 * Boschloo uses Fisher's one-sided p-value as the test statistic and takes
 * the supremum of the rejection probability over the nuisance parameter
 * π = P(success). Uniformly more powerful than Fisher's conditional exact
 * test on the same data, at the cost of higher computational expense.
 *
 * For preclinical incidence comparisons (one fixed margin: group sizes),
 * Boschloo is the recommended primary test; Fisher is the alternative.
 *
 * References:
 *   - Boschloo R. D. (1970). "Raised conditional level of significance for
 *     the 2×2-table when testing the equality of two probabilities."
 *     Statistica Neerlandica 24(1), 1–9.
 *   - SciPy documentation for `scipy.stats.boschloo_exact`.
 *   - Saari L. M. et al. (2004). "Employee attitudes and job satisfaction."
 *     Human Resource Management 43(4), 395–407 (used in the SciPy docs).
 *
 * Run: npx tsx src/stats/examples/boschloo-exact.ts
 */

import {boschlooExact, incidenceExactBoth, fisherExact2x2} from '..';
import {fmt, printCase, printKv, printProblem, printSection} from './_helpers';

/* ================================================================== */
/*  Demo                                                              */
/* ================================================================== */

printSection('BOSCHLOO 2×2 EXACT TEST');

/* ----- SciPy doc example (Saari 2004) ----- */

printProblem([
  'SciPy doc example (Saari 2004): job-satisfaction survey.',
  'Table: [[74, 31], [43, 32]]',
  'Expected (alternative = greater): statistic ≈ 0.04831, p ≈ 0.03556.',
  '',
]);

printCase('Boschloo, alternative = greater');
const saari = boschlooExact([[74, 31], [43, 32]], {alternative: 'greater'});
printKv('statistic', fmt(saari.statistic, 6));
printKv('p-value', fmt(saari.pValue, 6));

/* ----- Preclinical 2/3 vs 0/3 ----- */

printProblem([
  '',
  'Preclinical 2/3 vs 0/3 incidence (group 1: 2 successes in 3 animals,',
  'group 2: 0 successes in 3). Cols-as-groups layout: col 0 = (2 success,',
  '1 failure), col 1 = (0 success, 3 failures).',
  '  Fisher (conditional)  p = 0.40',
  '  Boschloo (unconditional) p ≈ 0.22 — power gain is real, but',
  '  2/3 vs 0/3 is still not significant at α = 0.05.',
  'Table: [[2, 0], [1, 3]]',
  '',
]);

printCase('Fisher exact (conditional)');
const fisherSmall = fisherExact2x2([[2, 0], [1, 3]]);
printKv('odds ratio', fmt(fisherSmall.oddsRatio));
printKv('p (two-sided)', fmt(fisherSmall.pValue, 6));

printCase('Boschloo exact (unconditional)');
const boschlooSmall = boschlooExact([[2, 0], [1, 3]]);
printKv('statistic', fmt(boschlooSmall.statistic, 6));
printKv('p (two-sided)', fmt(boschlooSmall.pValue, 6));

/* ----- incidenceExactBoth: primary + alternative together ----- */

printProblem([
  '',
  '`incidenceExactBoth` returns Boschloo (primary) and Fisher (alternative)',
  'p-values together with the sample odds ratio. SEND policy: a degenerate',
  'Boschloo p (zero column → NaN) is reported as 1.0.',
  'Table: [[10, 3], [5, 12]]',
  '',
]);

printCase('Boschloo + Fisher');
const both = incidenceExactBoth([[10, 3], [5, 12]]);
printKv('odds ratio', fmt(both.oddsRatio, 4));
printKv('p (Boschloo)', fmt(both.pValue, 6));
printKv('p (Fisher)', fmt(both.pValueFisher, 6));

/* ----- Edge cases ----- */

printSection('EDGE CASES');

printProblem([
  'Zero column total → both p-values are NaN by SciPy convention.',
  '`incidenceExactBoth` reports the degenerate Boschloo p as 1.0.',
  'Table: [[0, 5], [0, 5]]',
  '',
]);

printCase('Degenerate column (raw boschlooExact)');
const deg = boschlooExact([[0, 5], [0, 5]]);
printKv('statistic', String(deg.statistic));
printKv('p-value', String(deg.pValue));

printCase('Degenerate column (incidenceExactBoth)');
const degBoth = incidenceExactBoth([[0, 5], [0, 5]]);
printKv('odds ratio', String(degBoth.oddsRatio));
printKv('p (Boschloo)', fmt(degBoth.pValue));
printKv('p (Fisher)', String(degBoth.pValueFisher));
