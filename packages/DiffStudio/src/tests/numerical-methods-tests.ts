// Tests of numerical methods
import {_package} from '../package-test';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {mrt, ros3prw, ros34prw, CorrProblem, ODEs, corrProbs, perfProbs,
  rk3, rk4, rkdp, ab4, ab5, lsoda, cvode} from 'diff-grok';

/** Return numerical solution error: maximum absolute deviation between approximate & exact solutions */
function getError(method: (odes: ODEs) => Float64Array[], corProb: CorrProblem): number {
  let error = 0;

  // Get numerical solution
  const approxSolution = method(corProb.odes);

  const exact = corProb.exact;

  const arg = approxSolution[0];

  const pointsCount = arg.length;
  const funcsCount = approxSolution.length - 1;

  // Compute error
  for (let i = 0; i < pointsCount; ++i) {
    const exactSolution = exact(arg[i]);

    for (let j = 0; j < funcsCount; ++j)
      error = Math.max(error, Math.abs(exactSolution[j] - approxSolution[j + 1][i]));
  }

  return error;
}

const TIMEOUT = 4000;
const TINY = 0.1;
const MIN_ROWS = 1000;

const implicitMethods = new Map([
  ['MRT', mrt],
  ['ROS3PRw', ros3prw],
  ['ROS34PRw', ros34prw],
  ['LSODA', lsoda],
  ['CVODE', cvode],
]);

/** Performance problems to skip per method (known limitations) */
const perfExclusions = new Map<string, Set<string>>([
  ['CVODE', new Set(['E5'])],
]);

const explicitMethods = new Map([
  ['RK3', rk3],
  ['RK4', rk4],
  ['RKDP', rkdp],
  ['AB4', ab4],
  ['AB5', ab5],
]);

const allMethods = new Map([...implicitMethods, ...explicitMethods]);

// Correctness tests
allMethods.forEach((method, name) => {
  category(`Correctness: ${name}`, () => {
    corrProbs.forEach((problem) => test(problem.odes.name, async () => {
      const error = getError(method, problem);
      expect(
        error < TINY,
        true,
        `The ${name} method failed to solve "${problem.odes.name}", too big error: ${error}; expected: < ${TINY}`,
      );
    }, {timeout: TIMEOUT}));
  });
}); // Correctness

// Performance tests
implicitMethods.forEach((method, methodName) => {
  category(`Performance: ${methodName}`, () => {
    const excluded = perfExclusions.get(methodName);
    perfProbs.forEach((odes) => {
      if (excluded?.has(odes.name))
        return;

      test(odes.name, async () => {
        const rows = method(odes)[0].length;
        expect(
          rows > MIN_ROWS,
          true,
          `The ${methodName} method failed, solution DF rows: ${rows}`,
        );
      }, {benchmark: true});
    });
  });
}); // Performance
