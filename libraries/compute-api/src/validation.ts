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
export const makeValidationResult = (window.compute.makeValidationResult) as
  ((...args: Parameters<typeof makeValidationResultType>) => ReturnType<typeof makeValidationResultType>);

//@ts-ignore
export const makeAdvice = (window.compute.makeAdvice) as
  ((...args: Parameters<typeof makeAdviceType>) => ReturnType<typeof makeAdviceType>);

//@ts-ignore
export const makeRevalidation = (window.compute.makeRevalidation) as
  ((...args: Parameters<typeof makeRevalidationType>) => ReturnType<typeof makeRevalidationType>);

export type ValidationInfo = ValidationInfoType;

