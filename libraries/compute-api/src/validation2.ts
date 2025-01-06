/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import type {
  makeValidationResult as makeValidationResultType,
  makeAdvice as makeAdviceType,
  mergeValidationResults as mergeValidationResultsType,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver';

export function makeValidationResult2(
  ...args: Parameters<typeof makeValidationResultType>
): ReturnType<typeof makeValidationResultType> {
  return window.compute.makeValidationResult2(...args);
}

export function makeAdvice2(
  ...args: Parameters<typeof makeAdviceType>
): ReturnType<typeof makeAdviceType> {
  return window.compute.makeAdvice2(...args);
}

export function mergeValidationResults2(
  ...args: Parameters<typeof mergeValidationResultsType>
): ReturnType<typeof mergeValidationResultsType> {
  return window.compute.mergeValidationResults2(...args);
}
