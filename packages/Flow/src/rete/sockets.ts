/** Typed socket wrapping a Datagrok type string; compatibility via `areTypesCompatible`. */
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

/** One shared TypedSocket per dgType so reference equality holds. */
const _socketCache = new Map<string, TypedSocket>();

export function getSocket(dgType: string): TypedSocket {
  let s = _socketCache.get(dgType);
  if (!s) {
    s = new TypedSocket(dgType);
    _socketCache.set(dgType, s);
  }
  return s;
}
