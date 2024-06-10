/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {
  makeValidationResult as makeValidationResultType,
  makeAdvice as makeAdviceType,
  makeRevalidation as makeRevalidationType,
  ValidationInfo as ValidationInfoType,
} from '@datagrok-libraries/compute-utils';

//@ts-ignore
export function makeValidationResult(
  ...args: Parameters<typeof makeValidationResultType>
): ReturnType<typeof makeValidationResultType> {
  //@ts-ignore
  return window.compute.makeValidationResult(...args);
}

//@ts-ignore
export function makeAdvice(
  ...args: Parameters<typeof makeAdviceType>
): ReturnType<typeof makeAdviceType> {
  //@ts-ignore
  return window.compute.makeAdvice(...args);
}

//@ts-ignore
export function makeRevalidation(
  ...args: Parameters<typeof makeRevalidationType>
): ReturnType<typeof makeRevalidationType> {
  //@ts-ignore
  return window.compute.makeRevalidation(...args);
}

export type ValidationInfo = ValidationInfoType;

