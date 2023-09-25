import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export interface SubscriptionLike {
  unsubscribe(): void;
}

// an arbitrary input single input bindable form
export interface FuncCallInput<T = any> {
  root: HTMLElement;
  value: T | null | undefined;
  notify: boolean;
  enabled: boolean;
  onInput: (cb: Function) => SubscriptionLike;
}

// runtime checker
export function isFuncCallInput<T = any>(arg: any): arg is FuncCallInput<T> {
  return arg && arg.root && arg.onInput && arg.hasOwnProperty('value');
}

// for an input wrapper with a single standard input
export interface InputWrapper<T = any> {
  primaryInput: DG.InputBase<T>;
}
