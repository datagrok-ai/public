import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export interface SubscriptionLike {
  unsubscribe(): void;
}

export interface FuncCallInput<T = any> {
  root: HTMLElement;
  value: T | null | undefined;
  notify: boolean;
  enabled: boolean;
  onInput: (cb: Function) => SubscriptionLike;
}

export interface FuncCallInputValidated<T = any> extends FuncCallInput<T> {
  setValidation: (messages?: ValidationResult) => void;
}

export type InputFactory = (params: any) => { input: FuncCallInput | FuncCallInputValidated };

// runtime checkers
export function isFuncCallInput<T = any>(arg: any): arg is FuncCallInput<T> {
  return arg && arg.root && arg.onInput;
}

export function isFuncCallInputValidated<T = any>(arg: any): arg is FuncCallInputValidated<T> {
  return arg?.setValidation && isFuncCallInput(arg);
}

// validation
export interface ValidationResult {
  warnings?: string[];
  errors?: string[];
}

export function isValidationPassed(result?: ValidationResult) {
  return !result || !result.errors?.length;
}

export function getErrorMessage(result?: ValidationResult) {
  if (result?.errors)
    return result.errors.join('; ');
}

export function getWarningMessage(result?: ValidationResult) {
  if (result?.warnings)
    return result.warnings.join('; ');
}

export function makeValidationResult(errors?: string[], warnings?: string[]) {
  return {errors, warnings};
}

export function mergeValidationResults(results: ValidationResult[] = []) {
  const errors = results.flatMap((res) => res.errors);
  const warnings = results.flatMap((res) => res.warnings);
  return {errors, warnings};
}

export type Validator = (val: any, context?: any) => ValidationResult | void;

export type ValidatorFactory = (params: any) => { validator: Validator };

export const nonNullValidator: Validator = (value: any) => {
  if (value == null)
    return makeValidationResult(['Missing value']);
};
