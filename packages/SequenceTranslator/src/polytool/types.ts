import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export enum PolyToolEnumeratorTypes {
  Single = 'single',
  Matrix = 'matrix',
}

export type PolyToolEnumeratorType = typeof PolyToolEnumeratorTypes[keyof typeof PolyToolEnumeratorTypes];

export type PolyToolPlaceholders = { position: number, monomers: string[] } [];

export type PolyToolPlaceholdersBreadth = { start: number, end: number, monomers: string[] }[];

export type PolyToolEnumeratorParams = {
  type: PolyToolEnumeratorType;
  /** position key is zero-based */
  placeholders?: PolyToolPlaceholders;
  placeholdersBreadth?: PolyToolPlaceholdersBreadth;
  keepOriginal?: boolean;
  trivialName?: boolean;
}
