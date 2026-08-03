/**
 * Example: rank-based correlation between two variables and a thin
 * dose-response wrapper (`severityTrend`) that guards against constant
 * severity (an undefined-correlation degeneracy that arises in toxicology
 * when a dose group has no observed effect).
 *
 * References:
 *   - Spearman C (1904). "The proof and measurement of association
 *     between two things." Am J Psychol 15, 72-101.
 *   - IQ vs TV-hours dataset: Wikipedia, "Spearman's rank correlation
 *     coefficient" worked example.
 *   - Clinical-psychology rankings: Armitage P, Berry G (1994). Statistical
 *     Methods in Medical Research, 3rd ed., Blackwell Science, p. 466.
 *
 * Run: npx tsx src/stats/examples/correlation.ts
 */

import {severityTrend, spearman} from '..';
import {fmt, printCase, printKv, printProblem, printSection} from './_helpers';

/* ================================================================== */
/*  Datasets                                                          */
/* ================================================================== */

// Wikipedia "Spearman's rank correlation coefficient" example: IQ vs
// TV viewing hours per week. Expected ρ = -29/165 ≈ -0.1758.
const WIKI_IQ = [106, 86, 100, 101, 99, 103, 97, 113, 112, 110];
const WIKI_TV = [7, 0, 27, 50, 28, 29, 20, 12, 6, 17];

// Armitage & Berry (1994) p. 466: ten students ranked by interest in
// a clinical-psychology career and by aptitude. Σd² = 52, ρ = 113/165 ≈ 0.6848.
const AB94_CAREER = [4, 10, 3, 1, 9, 2, 6, 7, 8, 5];
const AB94_PSYCH = [5, 8, 6, 2, 10, 3, 9, 4, 7, 1];

// Toxicology-like dose-response: organ severity grade increases with dose
// but plateaus near the maximum.
const TOX_DOSES = [0, 50, 100, 200, 500];
const TOX_SEVERITIES = [0.0, 0.2, 0.8, 1.5, 1.6];

/* ================================================================== */
/*  Demo                                                              */
/* ================================================================== */

printSection('SPEARMAN RANK CORRELATION');

printProblem([
  'Wikipedia IQ vs TV-hours: 10 subjects, two ranking variables.',
  'Expected: ρ = -29/165 ≈ -0.1758, non-significant.',
  '',
]);

printCase('Spearman ρ');
const wiki = spearman(WIKI_IQ, WIKI_TV);
printKv('rho', fmt(wiki.rho));
printKv('p-value', fmt(wiki.pValue));

printProblem([
  '',
  'Armitage & Berry (1994) p.466: clinical-psychology rankings.',
  'Expected: ρ = 113/165 ≈ 0.6848, p ≈ 0.029.',
  '',
]);

printCase('Spearman ρ');
const ab94 = spearman(AB94_CAREER, AB94_PSYCH);
printKv('rho', fmt(ab94.rho));
printKv('p-value', fmt(ab94.pValue));

printSection('SEVERITY TREND (DOSE-RESPONSE)');

printProblem([
  'Toxicology dose-response: average severity by dose level.',
  'doses:      [0, 50, 100, 200, 500]',
  'severities: [0.0, 0.2, 0.8, 1.5, 1.6]',
  'Expected: strong positive trend.',
  '',
]);

printCase('severityTrend');
const tox = severityTrend(TOX_DOSES, TOX_SEVERITIES);
printKv('rho', fmt(tox.rho));
printKv('p-value', fmt(tox.pValue));

printProblem([
  '',
  'Edge case: severity is constant across all doses (no observed effect).',
  'Expected: rho = null (correlation undefined).',
  '',
]);

printCase('severityTrend with constant severity');
const flat = severityTrend([0, 10, 50, 100], [1.0, 1.0, 1.0, 1.0]);
printKv('rho', fmt(flat.rho));
printKv('p-value', fmt(flat.pValue));
