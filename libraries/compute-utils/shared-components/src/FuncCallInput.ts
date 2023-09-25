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

// validation/advisory system
export interface ActionItems {
  actionName: string;
  action: Function;
}

export interface Advice {
  description: string;
  actions?: ActionItems[];
}

export interface ValidationResultBase {
  // awaiting for validation results
  pending?: boolean;
  errors?: Advice[];
  warnings?: Advice[];
  notifications?: Advice[];
}

export interface ValidationResult extends ValidationResultBase {
  // revalidation request
  revalidate?: string[];
  // revalidations context
  context?: any;
}

export function isValidationPassed(result?: ValidationResult) {
  return !result?.errors?.length && !result?.pending;
}

export function makeAdvice(description: string, actions?: ActionItems[]) {
  return {description, actions};
}

export function getErrorMessage(result?: ValidationResult) {
  if (result?.errors)
    return result.errors.map((err) => err.description).join('; ');
}

export interface ValidationPayload {
  errors?: (string | Advice)[],
  warnings?: (string | Advice)[],
  notifications?: (string | Advice)[],
}

export function makeValidationResult(payload?: ValidationPayload): ValidationResultBase {
  const wrapper = (item: string | Advice) => typeof item === 'string' ? makeAdvice(item) : item;
  return {
    errors: payload?.errors?.map((err) => wrapper(err)),
    warnings: payload?.warnings?.map((warn) => wrapper(warn)),
    notifications: payload?.notifications?.map((note) => wrapper(note)),
  };
}

export function makePendingValidationResult(): ValidationResult {
  return {pending: true};
}

export function makeRevalidation(revalidate: string[], context?: any, result?: ValidationResultBase): ValidationResult {
  return {revalidate, context, ...result};
}

export interface ValidationInfo {
  param: string,
  funcCall: DG.FuncCall,
  lastCall?:DG.FuncCall,
  isRevalidation: boolean,
  isNewOutput: boolean,
  context?: any
}

export type Validator = (val: any, info: ValidationInfo)
  => Promise<ValidationResult | undefined>;

export type ValidatorFactory = (params: any) => { validator: Validator };

export const nonNullValidator: Validator = async (value: any) => {
  if (value == null)
    return makeValidationResult({errors: ['Missing value']});
};
