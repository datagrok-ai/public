// Tests of numerical methods

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/utils/src/test';

import {mrt} from '../solver-tools/mrt-method';
import {ros3prw} from '../solver-tools/ros3prw-method';
import {ros34prw} from '../solver-tools/ros34prw-method';

import {evaluateMethod} from './testing-solvers-utils';
import {problems} from './performance-problems';

const TIMEOUT = 4000;
const TINY = 0.1;
const MIN_ROWS = 1000;

const methods = new Map([
  ['MRT', mrt],
  ['ROSP3PRw', ros3prw],
  ['ROSP34PRw', ros34prw],
]);

// Correctness tests
category('Correctness', () => {
  methods.forEach((method, name) => test(`The ${name} method`, async () => {
    const error = evaluateMethod(method);
    expect(
      error < TINY,
      true,
      `The ${name} method failed, too big error: ${error}; expected: < ${TINY}`,
    );
  }, {timeout: TIMEOUT}));
}); // Correctness

// Performance tests
category('Performance', () => {
  methods.forEach((method, methodName) => {
    problems.forEach((odes, odesName) => {
      test(`Method: ${methodName}, problem: ${odesName}`, async () => {
        const rows = method(odes).rowCount;
        expect(
          rows > MIN_ROWS,
          true,
          `The ${name} method failed, solution DF rows: ${rows}`,
        );
      }, {benchmark: true});
    });
  });
}); // Performance
