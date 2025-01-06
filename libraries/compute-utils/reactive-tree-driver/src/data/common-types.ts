
export type ItemId = string;
export type LinkSpecString = string | string[];
export type ItemPathArray = string[];
export type NqName = string;
export type TypeKey = string;
export type RestrictionType = 'disabled' | 'restricted' | 'info' | 'none';
export type StateType = 'input' | 'output' | 'state';
export type Constructor = new (...args: any[]) => {};
export type GConstructor<T = {}> = new (...args: any[]) => T;
export type TraverseHandler<R, I, A> = (acc: R, item: I, pathAddress: A, stop: () => void) => R;

export interface ActionItem {
  actionName: string;
  action: string;
  additionalParams?: Record<string, any>;
}

export interface Advice {
  description: string;
  actions?: ActionItem[];
}

export interface ValidationResult {
  errors?: Advice[];
  warnings?: Advice[];
  notifications?: Advice[];
}

export interface ValidationPayload {
  errors?: (string | Advice)[],
  warnings?: (string | Advice)[],
  notifications?: (string | Advice)[],
}
