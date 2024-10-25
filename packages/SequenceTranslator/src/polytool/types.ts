import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {PolymerType} from '@datagrok-libraries/bio/src/helm/types';

export enum PolyToolEnumeratorTypes {
  Single = 'single',
  Matrix = 'matrix',
}

export type PolyToolEnumeratorType = typeof PolyToolEnumeratorTypes[keyof typeof PolyToolEnumeratorTypes];

export type PolyToolPlaceholder = { position: number, monomers: string[] };

export type PolyToolBreadthPlaceholder = { start: number, end: number, monomers: string[] };

export type PolyToolEnumeratorParams = {
  type: PolyToolEnumeratorType;
  /** position key is zero-based */
  placeholders?: PolyToolPlaceholder[];
  breadthPlaceholders?: PolyToolBreadthPlaceholder[];
  keepOriginal?: boolean;
  trivialName?: boolean;
}

export class MonomerNotFoundError extends Error {
  public type = 'MonomerNotFoundError';

  constructor(polymerType: PolymerType, symbol: string, options?: ErrorOptions) {
    super(`Monomer '${symbol}' of polymer type '${polymerType}' not found`, options);
  }
}

export class InvalidReactionError extends Error {
  public type = 'InvalidReactionError';

  constructor(reaction: string, options?: ErrorOptions) {
    super(`Invalid reaction '${reaction}'.`);
  }
}
