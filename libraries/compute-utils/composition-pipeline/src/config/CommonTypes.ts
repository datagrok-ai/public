
export type ItemName = string;
export type ItemPath = ItemName[];
export type NqName = string;
export type TypeKey = string;
export type InputState = 'disabled' | 'restricted' | 'user input';
export type GroupState<T = any> = Record<string, T>;
export type StateType = 'input' | 'output' | 'state';
