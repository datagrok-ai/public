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
  setValidation: (messages?: ValidationResultBase) => void;
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
export interface ValidationResultBase {
  warnings?: string[];
  errors?: string[];
  // TODO: some implementation
  actions?: Record<string, Function>;
}

export interface ValidationResult extends ValidationResultBase {
  // awaiting for validation results
  pending?: boolean;
  // revalidation request
  revalidate?: string[];
  // revalidations context
  context?: any;
}

export function isValidationPassed(result?: ValidationResult) {
  return !result?.errors?.length && !result?.pending;
}

export function getErrorMessage(result?: ValidationResult) {
  if (result?.errors)
    return result.errors.join('; ');
}

export function getWarningMessage(result?: ValidationResult) {
  if (result?.warnings)
    return result.warnings.join('; ');
}

export function makeValidationResult(errors?: string[], warnings?: string[]): ValidationResult {
  return {errors, warnings};
}

export function makePendingValidationResult(): ValidationResult {
  return {pending: true};
}

export function makeRevalidation(revalidate: string[], context?: any): ValidationResult {
  return {revalidate, context};
}

export function mergeValidationResults(results: ValidationResult[] = []) {
  const errors = results.flatMap((res) => res.errors);
  const warnings = results.flatMap((res) => res.warnings);
  return {errors, warnings};
}

export type Validator =
  (val: any, param: string, fc: DG.FuncCall, isRevalidation: boolean, context?: any)
  => Promise<ValidationResult | undefined>;

export type ValidatorFactory = (params: any) => { validator: Validator, isAdvisory: boolean };

export const nonNullValidator: Validator = async (value: any) => {
  if (value == null)
    return makeValidationResult(['Missing value']);
};
