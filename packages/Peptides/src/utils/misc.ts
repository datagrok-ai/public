import * as DG from 'datagrok-api/dg';

import * as C from './constants';

export function stringToBool(str: string): boolean {
  return str === 'true' ? true : false;
}

export function getTypedArrayConstructor(
  maxNum: number): Uint8ArrayConstructor | Uint16ArrayConstructor | Uint32ArrayConstructor {
  return maxNum < 256 ? Uint8Array :
    maxNum < 65536 ? Uint16Array :
      Uint32Array;
}
