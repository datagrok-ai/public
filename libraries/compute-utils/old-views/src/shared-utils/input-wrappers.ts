import {Observable} from 'rxjs';
import {ValidationResultBase} from './validation';

export interface SubscriptionLike {
  unsubscribe(): void;
}

export interface FuncCallInput<T = any> {
  root: HTMLElement;
  value: T | null | undefined;
  notify: boolean;
  enabled: boolean;
  onInput: ((cb: Function) => SubscriptionLike) | Observable<T>;
}

export interface FuncCallInputValidated<T = any> extends FuncCallInput<T> {
  setValidation: (messages?: ValidationResultBase) => void;
}

export interface FuncCallInputLockable<T = any> extends FuncCallInput<T> {
  setDisabled: () => void;
  setRestricted: () => void;
  setRestrictedUnlocked: () => void;
  setInconsistent: () => void;
  setUserInput: () => void;
}

export type InputFactory = (params: any) => { input: FuncCallInput | FuncCallInputValidated };

// runtime checkers
export function isFuncCallInput<T = any>(arg: any): arg is FuncCallInput<T> {
  return arg && arg.root && arg.onInput;
}

export function isFuncCallInputValidated<T = any>(arg: any): arg is FuncCallInputValidated<T> {
  return arg?.setValidation && isFuncCallInput(arg);
}

export function isInputLockable(arg: any): arg is FuncCallInputLockable {
  return arg?.setDisabled &&
    arg?.setRestricted &&
    arg?.setRestrictedUnlocked &&
    arg?.setInconsistent &&
    arg?.setUserInput &&
    isFuncCallInput(arg);
}
