/**
 * Example: comparing each treated group against a control while
 * controlling the family-wise error rate (FWER).
 *
 * Two strategies on the same data:
 *   1. Welch's t-test for each treated vs control + Bonferroni correction.
 *   2. Dunnett's many-to-one test, which uses the joint distribution of
 *      the contrasts and is uniformly more powerful than Bonferroni-Welch
 *      when the within-group variance is homogeneous.
 *
 * The two methods give similar significance decisions on this dataset,
 * but Dunnett's adjusted p-values are tighter — that is the point of the
 * comparison.
 *
 * References:
 *   - Dunnett CW (1955). "A multiple comparison procedure for comparing
 *     several treatments with a control." J Am Stat Assoc 50, 1096-1121.
 *   - Bonferroni CE (1936). "Teoria statistica delle classi e calcolo
 *     delle probabilità." Pubbl R Ist Super Sci Econ Comm Firenze 8, 3-62.
 *   - Dataset: Dunnett (1955), Section 4 — blood-count measurements
 *     (millions of cells/mm³) from 6 controls and two drugs.
 *
 * Run: npx tsx src/stats/examples/pairwise-vs-control.ts
 */

import {bonferroniCorrect, dunnettPairwise, welchPairwise} from '..';
import {fmt, printCase, printKv, printProblem, printSection} from './_helpers';

/* ================================================================== */
/*  Dataset                                                           */
/* ================================================================== */

// Dunnett (1955) Section 4: blood counts (millions cells/mm³).
// 6 controls, 4 subjects on drug A, 5 on drug B.
const CONTROL = [7.40, 8.50, 7.20, 8.24, 9.84, 8.32];
const DRUG_A = [9.76, 8.80, 7.68, 9.36];
const DRUG_B = [12.80, 9.68, 12.16, 9.20, 10.55];

const TREATED = [
  {doseLevel: 1, values: DRUG_A},
  {doseLevel: 2, values: DRUG_B},
];

/* ================================================================== */
/*  Demo                                                              */
/* ================================================================== */

printSection('PAIRWISE COMPARISONS VS CONTROL (FWER-CONTROLLED)');

printProblem([
  'Dunnett 1955 dataset: blood counts (millions cells/mm³).',
  'Control n = 6, Drug A n = 4, Drug B n = 5.',
  'Goal: identify which drug differs significantly from control while',
  'keeping the family-wise error rate at α = 0.05.',
  '',
]);

/* ----- Strategy 1: Welch + Bonferroni ----- */

printCase('Strategy 1 — Welch t-test + Bonferroni correction');
const welch = welchPairwise(CONTROL, TREATED);
const rawP = welch.map((w) => w.pValueWelch);
const adjP = bonferroniCorrect(rawP);

console.log('  dose      raw p (Welch)     adj p (Bonferroni)');
console.log('  -----     -------------     ------------------');
for (let i = 0; i < welch.length; i++) {
  const dose = String(welch[i].doseLevel).padEnd(8);
  const raw = fmt(rawP[i], 4).padEnd(17);
  const adj = fmt(adjP[i], 4);
  console.log(`  ${dose}  ${raw} ${adj}`);
}

/* ----- Strategy 2: Dunnett ----- */

printCase('Strategy 2 — Dunnett\'s many-to-one test');
const dun = dunnettPairwise(CONTROL, TREATED);
console.log('  dose      t-statistic       adj p (Dunnett)        Hedges g');
console.log('  -----     -----------       ---------------        --------');
for (const r of dun) {
  const dose = String(r.doseLevel).padEnd(8);
  const t = fmt(r.statistic, 4).padEnd(17);
  const p = fmt(r.pValueAdj, 4).padEnd(22);
  const g = fmt(r.effectSize, 4);
  console.log(`  ${dose}  ${t} ${p} ${g}`);
}

/* ----- Summary ----- */

printCase('Decision at α = 0.05');
for (let i = 0; i < dun.length; i++) {
  const dose = dun[i].doseLevel;
  const sigBon = adjP[i] !== null && adjP[i]! < 0.05;
  const sigDun = dun[i].pValueAdj !== null && dun[i].pValueAdj! < 0.05;
  printKv(
    `dose ${dose}`,
    `Bonferroni-Welch ${sigBon ? 'sig' : 'NS'}, ` +
    `Dunnett ${sigDun ? 'sig' : 'NS'}`,
    16,
  );
}
