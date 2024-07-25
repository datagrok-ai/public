
export type ItemId = string;
export type ItemPath = string;
export type ItemPathArray = string[];
export type ItemAddress = string[];
export type NqName = string;
export type TypeKey = string;
export type InputState = 'disabled' | 'restricted' | 'user input';
export type StateType = 'input' | 'output' | 'state';
export type Constructor = new (...args: any[]) => {};
export type GConstructor<T = {}> = new (...args: any[]) => T;
export type TraverseHandler<R, I, A> = (acc: R, item: I, pathAddress: A, stop: () => void) => R;
