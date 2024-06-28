/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {
  makeValidationResult as makeValidationResultType,
  makeAdvice as makeAdviceType,
  makeRevalidation as makeRevalidationType,
  mergeValidationResults as mergeValidationResultsType,
  ValidationInfo as ValidationInfoType,
} from '@datagrok-libraries/compute-utils';

export function makeValidationResult(
  ...args: Parameters<typeof makeValidationResultType>
): ReturnType<typeof makeValidationResultType> {
  return window.compute.makeValidationResult(...args);
}

export function makeAdvice(
  ...args: Parameters<typeof makeAdviceType>
): ReturnType<typeof makeAdviceType> {
  return window.compute.makeAdvice(...args);
}

export function makeRevalidation(
  ...args: Parameters<typeof makeRevalidationType>
): ReturnType<typeof makeRevalidationType> {
  return window.compute.makeRevalidation(...args);
}

export function mergeValidationResults(
  ...args: Parameters<typeof mergeValidationResultsType>
): ReturnType<typeof mergeValidationResultsType> {
  return window.compute.mergeValidationResults(...args);
}

export type ValidationInfo = ValidationInfoType;

