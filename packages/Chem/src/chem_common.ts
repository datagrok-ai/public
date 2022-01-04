import BitArray from '@datagrok-libraries/utils/src/bit-array';

export const defaultMorganFpRadius = 2;
export const defaultMorganFpLength = 2048;

const lockPromiseForKey: any = {};
const unlockFunctionForKey: any = {};

// By https://github.com/mistval/locko

export async function criticalSectionBegin(key: string) {
  if (!lockPromiseForKey[key]) {
    lockPromiseForKey[key] = Promise.resolve();
  }
  const takeLockPromise = lockPromiseForKey[key];
  lockPromiseForKey[key] = takeLockPromise.then(() => new Promise((fulfill) => {
    unlockFunctionForKey[key] = fulfill;
  }));
  return takeLockPromise;
}

export function criticalSectionEnd(key: string) {
  if (unlockFunctionForKey[key]) {
    unlockFunctionForKey[key]();
    delete unlockFunctionForKey[key];
  }
}

let _chemLocked = false;

export async function chemLock(token: string | null = null) {
  if (_chemLocked) {
    throw (`RdKit Service usage locked:\n${(new Error()).stack}`);
  }
  _chemLocked = true;
}

export async function chemUnlock(token: string | null = null) {
  _chemLocked = false;
}

export function tanimoto(x: BitArray, y: BitArray) {
  const total = x.trueCount() + y.trueCount();
  if (total == 0) {
    return 1.0;
  }
  const common = x.andWithCountBits(y, true);
  return common / (total - common);
}

export function rdKitFingerprintToBitArray(fp: string, fpLength: number) {
  const arr = new BitArray(fpLength);
  for (let j = 0; j < fpLength; ++j) {
    if (fp[j] === '1') {
      arr.setTrue(j);
    }
  }
  return arr;
}
