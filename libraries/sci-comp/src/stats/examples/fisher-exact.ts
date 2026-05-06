/**
 * Example: Fisher's exact test for a 2×2 contingency table.
 *
 * Conditions on row and column totals; the null distribution of the
 * top-left cell is hypergeometric. Returns the two-sided "minlike"
 * p-value (matches scipy default), both one-sided p-values, and the
 * sample odds ratio.
 *
 * References:
 *   - Fisher RA (1935). The Design of Experiments. Oliver and Boyd,
 *     Chapter 2 ("Lady Tasting Tea").
 *   - Agresti A (2002). Categorical Data Analysis, 2nd ed. Wiley, §3.5.
 *   - Saari LM et al. (2004). "Employee attitudes and job satisfaction."
 *     Human Resource Management 43(4), 395-407 (used in the SciPy docs).
 *
 * Run: npx tsx src/stats/examples/fisher-exact.ts
 */

import {fisherExact2x2} from '..';
import {fmt, printCase, printKv, printProblem, printSection} from './_helpers';

/* ================================================================== */
/*  Demo                                                              */
/* ================================================================== */

printSection('FISHER 2×2 EXACT TEST');

/* ----- Lady Tasting Tea ----- */

printProblem([
  'Fisher (1935) Chapter 2: the lady-tasting-tea experiment.',
  '8 cups, 4 milk-first, 4 tea-first; the lady identifies 3 of 4 correctly.',
  'Table layout (2×2):           [[3, 1],',
  '                                [1, 3]]',
  'Expected: two-sided p = 34/70 = 17/35 ≈ 0.4857 (exact rational),',
  '          OR = 9.0.',
  '',
]);

printCase('Fisher exact');
const tea = fisherExact2x2([[3, 1], [1, 3]]);
printKv('odds ratio', fmt(tea.oddsRatio, 6));
printKv('p (two-sided)', fmt(tea.pValue, 6));
printKv('p (greater)', fmt(tea.pGreater, 6));
printKv('p (less)', fmt(tea.pLess, 6));

/* ----- Saari (2004) employee survey ----- */

printProblem([
  '',
  'Saari et al. (2004): job-satisfaction survey, professors vs scientists.',
  'Table: [[74, 31], [43, 32]]',
  '',
]);

printCase('Fisher exact');
const saari = fisherExact2x2([[74, 31], [43, 32]]);
printKv('odds ratio', fmt(saari.oddsRatio));
printKv('p (two-sided)', fmt(saari.pValue, 6));
printKv('p (greater)', fmt(saari.pGreater, 6));

/* ----- Edge cases ----- */

printSection('EDGE CASES');

printProblem([
  'Perfect separation [[5, 0], [0, 5]]: OR = ∞, p (two-sided) = 1/126.',
  '',
]);

printCase('Perfect separation');
const perf = fisherExact2x2([[5, 0], [0, 5]]);
printKv('odds ratio', String(perf.oddsRatio));
printKv('p (two-sided)', fmt(perf.pValue, 8));

printProblem([
  '',
  'Null table [[5, 5], [5, 5]]: OR = 1, p (two-sided) ≈ 1.0.',
  '',
]);

printCase('Null table');
const nul = fisherExact2x2([[5, 5], [5, 5]]);
printKv('odds ratio', fmt(nul.oddsRatio));
printKv('p (two-sided)', fmt(nul.pValue));
