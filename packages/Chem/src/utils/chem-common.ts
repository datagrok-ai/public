import BitArray from '@datagrok-libraries/utils/src/bit-array';

export const defaultMorganFpRadius = 2;
export const defaultMorganFpLength = 2048;

const lockPromiseForKey: any = {};
const unlockFunctionForKey: any = {};

export enum Fingerprint {
  Morgan = 'Morgan',
  RDKit = 'RDKit',
  Pattern = 'Pattern'
}

// By https://github.com/mistval/locko

export async function criticalSectionBegin(key: string): Promise<any> {
  if (!lockPromiseForKey[key])
    lockPromiseForKey[key] = Promise.resolve();
  const takeLockPromise = lockPromiseForKey[key];
  lockPromiseForKey[key] = takeLockPromise.then(() => new Promise((fulfill) => {
    unlockFunctionForKey[key] = fulfill;
  }));
  return takeLockPromise;
}

export function criticalSectionEnd(key: string): void {
  if (unlockFunctionForKey[key]) {
    unlockFunctionForKey[key]();
    delete unlockFunctionForKey[key];
  }
}

const CHEM_TOKEN = 'CHEM_TOKEN';

export async function chemBeginCriticalSection(token = CHEM_TOKEN): Promise<void> {
  let warned = false;
  if (unlockFunctionForKey[token]) {
    console.warn('Chem | Is already in a critical section, waiting...');
    warned = true;
  }
  await criticalSectionBegin(token);
  if (warned)
    console.warn('Chem | Left the critical section');
}

export function chemEndCriticalSection(token = CHEM_TOKEN): void {
  criticalSectionEnd(token);
}

export function rdKitFingerprintToBitArray(fp: Uint8Array): BitArray | null {
  return fp ? BitArray.fromBytes(fp) : fp;
}

export function isMolBlock(s: string) {
  return s.includes('M  END');
}
