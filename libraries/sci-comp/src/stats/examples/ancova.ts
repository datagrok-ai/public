/**
 * Example: ANCOVA (analysis of covariance) for adjusting group means by
 * a continuous covariate.
 *
 * Fits the model `y ~ C(group) + covariate` via OLS and reports:
 *   - LS-adjusted means (group means at the overall covariate mean)
 *   - pairwise contrasts vs control (with t-statistic and p-value)
 *   - slope estimate and slope-homogeneity F-test
 *   - effect decomposition (total / direct / indirect with Hedges' g)
 *
 * The slope-homogeneity test is a *prerequisite check*: if slopes differ
 * across groups, the basic ANCOVA model is misspecified.
 *
 * References:
 *   - Montgomery DC (2012). Design and Analysis of Experiments, 8th ed.
 *     Wiley, Table 15.10 — fiber strength by machine, with diameter as
 *     a continuous covariate.
 *   - Lazic SE et al. (2020). "Improving basic and translational science
 *     by accounting for litter-to-litter variation and batch effects."
 *     Sci Rep 10, 6625 (organ-free body-weight covariate).
 *
 * Run: npx tsx src/stats/examples/ancova.ts
 */

import {runAncova} from '..';
import {fmt, printCase, printKv, printProblem, printSection} from './_helpers';

/* ================================================================== */
/*  Dataset — Montgomery Table 15.10                                  */
/* ================================================================== */

// 3 machines × 5 observations: y = fiber strength (lbs),
// x = diameter (covariate). Reproduces SAS PROC GLM output to within
// 10⁻⁵.
const Y = [36, 41, 39, 42, 49, 40, 48, 39, 45, 44, 35, 37, 42, 34, 32];
const X = [20, 25, 24, 25, 32, 22, 28, 22, 30, 28, 21, 23, 26, 21, 15];
const G = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3];

/* ================================================================== */
/*  Demo                                                              */
/* ================================================================== */

printSection('ANCOVA — MONTGOMERY TABLE 15.10');

printProblem([
  'Montgomery DC (2012) Table 15.10: fiber strength by machine, diameter',
  'as covariate. 3 machines × 5 observations each.',
  'Expected per SAS PROC GLM:',
  '  R² = 0.9192, MSE = 2.5442',
  '  LS means: machine 1 = 40.382, machine 2 = 41.419, machine 3 = 38.798',
  '  slope t = 8.365 (highly significant)',
  '  slope homogeneity: F = 0.49, p = 0.6293 (slopes are homogeneous)',
  '',
]);

printCase('runAncova (control = machine 1)');
const r = runAncova(Y, X, G, {controlGroup: 1, alpha: 0.05});
if (r === null) throw new Error('ANCOVA returned null');

printKv('R²', fmt(r.modelRSquared));
printKv('MSE', fmt(r.mse));
printKv('covariate mean', fmt(r.covariateMean));

console.log('\n  Adjusted (LS) means:');
console.log('    group  n   raw mean   adj mean   se');
console.log('    -----  --  --------   --------   ------');
for (const am of r.adjustedMeans) {
  console.log(
    `    ${String(am.group).padEnd(5)}  ${String(am.n).padEnd(2)}  ` +
    `${fmt(am.rawMean).padEnd(9)}  ${fmt(am.adjustedMean).padEnd(9)}  ` +
    `${fmt(am.se)}`,
  );
}

console.log('\n  Pairwise vs control (machine 1):');
console.log('    group  diff       se       t        p          sig');
console.log('    -----  --------   ------   ------   --------   ---');
for (const pw of r.pairwise) {
  const sig = pw.significant ? 'YES' : ' no';
  console.log(
    `    ${String(pw.group).padEnd(5)}  ${fmt(pw.difference).padEnd(9)}  ` +
    `${fmt(pw.se).padEnd(7)}  ${fmt(pw.tStatistic).padEnd(7)}  ` +
    `${fmt(pw.pValue, 6).padEnd(9)}  ${sig}`,
  );
}

printCase('Slope (regression coefficient)');
printKv('estimate', fmt(r.slope.estimate, 6));
printKv('se', fmt(r.slope.se, 6));
printKv('t', fmt(r.slope.tStatistic));
printKv('p', fmt(r.slope.pValue, 6));

printCase('Slope homogeneity F-test (precondition check)');
printKv('F', fmt(r.slopeHomogeneity.fStatistic));
printKv('p', fmt(r.slopeHomogeneity.pValue));
printKv('homogeneous', String(r.slopeHomogeneity.homogeneous));

printCase('Effect decomposition (total = direct + indirect)');
console.log('    group  total      direct     indirect   prop_dir   |g|');
console.log('    -----  --------   --------   --------   --------   ------');
for (const ed of r.effectDecomposition) {
  console.log(
    `    ${String(ed.group).padEnd(5)}  ${fmt(ed.totalEffect).padEnd(9)}  ` +
    `${fmt(ed.directEffect).padEnd(9)}  ${fmt(ed.indirectEffect).padEnd(9)}  ` +
    `${fmt(ed.proportionDirect).padEnd(9)}  ${fmt(ed.directG)}`,
  );
}

/* ================================================================== */
/*  Lazic-style organ-free covariate                                  */
/* ================================================================== */

printSection('LAZIC (2020) ORGAN-FREE COVARIATE');

printProblem([
  'Lazic et al. (2020) recommend using BW − organ_weight as the covariate',
  'when the response itself is part of the body. The flag rewires the',
  'covariate inside the model and reports its mean.',
  '',
]);

printCase('runAncova with useOrganFreeBw = true');
const rLazic = runAncova(Y, X, G, {controlGroup: 1, useOrganFreeBw: true});
if (rLazic !== null) {
  printKv('covariate mean (BW − Y)', fmt(rLazic.covariateMean));
  printKv('useOrganFreeBw', String(rLazic.useOrganFreeBw));
  printKv('R² (organ-free)', fmt(rLazic.modelRSquared));
  printKv('R² (BW only)', fmt(r.modelRSquared));
}
