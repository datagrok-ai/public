export function stringToBool(str: string): boolean {
  switch (str) {
    case 'true':
      return true;
    case 'false':
      return false;
    default:
      throw new Error(`StrToBoolError: cannot convert string '${str}' to boolean`);
  }
}

export function getTypedArrayConstructor(
  maxNum: number): Uint8ArrayConstructor | Uint16ArrayConstructor | Uint32ArrayConstructor {
  return maxNum < 256 ? Uint8Array :
    maxNum < 65536 ? Uint16Array :
      Uint32Array;
}
