import BitArray from '@datagrok-libraries/utils/src/bit-array';

export const defaultMorganFpRadius = 2;
export const defaultMorganFpLength = 2048;

const lockPromiseForKey: any = {};
const unlockFunctionForKey: any = {};

// By https://github.com/mistval/locko

export async function criticalSectionBegin(key: string) {
  if (!lockPromiseForKey[key])
    lockPromiseForKey[key] = Promise.resolve();
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

const CHEM_TOKEN = 'CHEM_TOKEN';

export async function chemBeginCriticalSection(token = CHEM_TOKEN) {
  let warned = false;
  if (unlockFunctionForKey[token]) {
    console.warn('Chem | Is already in a critical section, waiting...');
    warned = true;
  }
  await criticalSectionBegin(token);
  if (warned) {
    console.warn('Chem | Left the critical section');
  }
}

export function chemEndCriticalSection(token = CHEM_TOKEN) {
  criticalSectionEnd(token);
}

export function tanimoto(x: BitArray, y: BitArray) {
  const total = x.trueCount() + y.trueCount();
  if (total == 0)
    return 1.0;
  const common = x.andWithCountBits(y, true);
  return common / (total - common);
}

export function rdKitFingerprintToBitArray(fp: string) {
  const length = fp.length;
  const arr = new BitArray(length);
  for (let j = 0; j < length; ++j) {
    if (fp[j] === '1')
      arr.setTrue(j);
  }
  return arr;
}
