/** Typed socket for FuncFlow. Wraps a Datagrok type string and routes
 * compatibility checks through `areTypesCompatible` from `type-map`. */
import {ClassicPreset} from 'rete';
import {areTypesCompatible} from '../types/type-map';

export class TypedSocket extends ClassicPreset.Socket {
  constructor(public dgType: string) {
    super(dgType);
  }

  /** Whether this socket (output side) can connect to `other` (input side). */
  isCompatibleWith(other: TypedSocket): boolean {
    return areTypesCompatible(this.dgType, other.dgType);
  }
}

/** Cache of TypedSocket instances by dgType — sockets are pure value-objects
 * so we share one per type. Keeps reference equality where possible. */
const _socketCache = new Map<string, TypedSocket>();

export function getSocket(dgType: string): TypedSocket {
  let s = _socketCache.get(dgType);
  if (!s) {
    s = new TypedSocket(dgType);
    _socketCache.set(dgType, s);
  }
  return s;
}
