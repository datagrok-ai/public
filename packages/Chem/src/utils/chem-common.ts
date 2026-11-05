/** Chem's common helpers, engine-neutral: those that take a module (`RDModule`, RDKit.js's shape) work on either
 * engine's, as their callers pass it. */
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {getMolSafe} from './mol-creation_rdkit';
import {MAX_SMILES_LENGTH} from './chem-constants';
import {MolList, RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';

export const defaultMorganFpRadius = 2;
export const defaultMorganFpLength = 2048;

/**
 * Per-key promise-chain mutex (adapted from https://github.com/mistval/locko).
 * `begin(key)` returns a promise that resolves when the previous holder of `key`
 * calls `end(key)`; `end(key)` hands the lock to the next waiter. State is
 * confined to the class so callers can't reach into it from anywhere in the file.
 */
class LockManager {
  private lockPromiseForKey = new Map<string, Promise<void>>();
  private unlockFunctionForKey = new Map<string, () => void>();

  async begin(key: string): Promise<void> {
    const takeLockPromise = this.lockPromiseForKey.get(key) ?? Promise.resolve();
    this.lockPromiseForKey.set(key, takeLockPromise.then(() => new Promise<void>((fulfill) => {
      this.unlockFunctionForKey.set(key, fulfill);
    })));
    return takeLockPromise;
  }

  end(key: string): void {
    const unlock = this.unlockFunctionForKey.get(key);
    if (unlock) {
      unlock();
      this.unlockFunctionForKey.delete(key);
    }
  }

  has(key: string): boolean {
    return this.unlockFunctionForKey.has(key);
  }
}

const chemLockManager = new LockManager();

export enum Fingerprint {
  Morgan = 'Morgan',
  RDKit = 'RDKit',
  Pattern = 'Pattern',
  //Avalon = 'Avalon',
  MACCS = 'MACCS',
  AtomPair = 'AtomPair',
  TopologicalTorsion = 'TopologicalTorsion'
}

export async function criticalSectionBegin(key: string): Promise<void> {
  return chemLockManager.begin(key);
}

export function criticalSectionEnd(key: string): void {
  chemLockManager.end(key);
}

const CHEM_TOKEN = 'CHEM_TOKEN';

// opId of the operation currently holding the section, or null.
let runningOpId: string | null = null;

// Monotonic counter → unique operation ids.
let chemOpCounter = 0;

export function newChemOpId(prefix = 'op'): string {
  return `${prefix}-${++chemOpCounter}`;
}

// True for the op holding the section, not ones merely queued.
export function isChemOpRunning(opId: string): boolean {
  return runningOpId === opId;
}

export async function chemBeginCriticalSection(opId?: string, token = CHEM_TOKEN): Promise<void> {
  let warned = false;
  if (chemLockManager.has(token)) {
    console.warn('Chem | Is already in a critical section, waiting...');
    warned = true;
  }
  await criticalSectionBegin(token);
  // Mark this operation as the section holder so a cancel can tell it apart from queued operations.
  if (opId)
    runningOpId = opId;
  if (warned)
    console.warn('Chem | Left the critical section');
}

export function chemEndCriticalSection(opId?: string, token = CHEM_TOKEN): void {
  // Clear only if this op is the registered holder, so a non-cancellable / different op can't wipe a running op's id.
  if (opId && runningOpId === opId)
    runningOpId = null;
  criticalSectionEnd(token);
}

// Runs `work` inside the Chem critical section, always releasing it — so a concurrent cancel (which restarts
// a worker) can't kill the operation mid-flight. Pass `opId` for a cancellable operation (see runCancellableChemOp).
export async function withChemCriticalSection<T>(work: () => Promise<T>, opId?: string): Promise<T> {
  await chemBeginCriticalSection(opId);
  try {
    return await work();
  } finally {
    chemEndCriticalSection(opId);
  }
}

// Cancellable wrapper: undefined if cancelled, else the result.
export async function runCancellableChemOp<T>(
  opId: string, isCanceled: () => boolean, work: () => Promise<T>): Promise<T | undefined> {
  return withChemCriticalSection(async (): Promise<T | undefined> => {
    if (isCanceled())
      return undefined;
    try {
      return await work();
    } catch (e) {
      if (isCanceled())
        return undefined;
      throw e;
    }
  }, opId);
}

/** A fingerprint's bytes (either engine's) as a BitArray. */
export function rdKitFingerprintToBitArray(fp: Uint8Array): BitArray | null {
  return fp ? BitArray.fromBytes(fp) : null;
}

export function isMolBlock(s: string): boolean {
  return s.includes('M  END');
}

/** `m` as Chem hands it to Crux: a MOL block or a SMILES within MAX_SMILES_LENGTH (a longer one can run the wasm out of
 * memory); null for anything else and an empty cell. */
export function withinSmilesLimit(m: unknown): string | null {
  return typeof m === 'string' && m !== '' && (isMolBlock(m) || m.length <= MAX_SMILES_LENGTH) ? m : null;
}

export function hasNewLines(s: string): boolean {
  return s.includes('\n') || s.includes('\r');
}

export function hexToPercentRgb(hex: string): number[] | null {
  const result = hex.length === 7 ? /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex) :
    /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  return result ? [
    parseInt(result[1], 16) / 256,
    parseInt(result[2], 16) / 256,
    parseInt(result[3], 16) / 256,
    result.length > 4 ? parseInt(result[4], 16) / 256 : 0.3,
  ] : null;
}

export function stringArrayToMolList(molecules: string[], module: RDModule): MolList {
  const mols = new module.MolList();
  for (let i = 0; i < molecules.length; i++)
    appendMolList(molecules[i], mols, module);
  return mols;
}

export function appendMolList(molString: string | null, mols: MolList, module: RDModule) {
  let molSafe;
  try {
    molSafe = getMolSafe(molString ?? '', {}, module);
    if (molSafe.mol === null || molSafe.isQMol || !molSafe.kekulize) {
      molSafe.mol?.delete();
      molSafe = getMolSafe('', {}, module);
    }
    mols.append(molSafe.mol!);
  } finally {
    molSafe?.mol?.delete();
  }
}

export function getSigFigs(n: number, sig: number) {
  const mult = Math.pow(10, sig - Math.floor(Math.log(n) / Math.LN10) - 1);
  return Math.round(n * mult) / mult;
}

export function getUncommonAtomsAndBonds(molecule: string, mcsMol: RDMol | null, module: RDModule,
  col?: string): ISubstruct | null {
  const mol = getMolSafe(molecule, {}, module).mol;
  let uncommonAtomsBest: number[] = [];
  let uncommonBondsBest: number[] = [];
  const highlightAtomColors: { [key: string]: number[] } = {};
  const highlightBondColors: { [key: string]: number[] } = {};
  try {
    if (mol) {
      let rgb: number[] | null = null;
      if (col)
        rgb = hexToPercentRgb(col);
      const uncommonAtoms = [...Array(mol.get_num_atoms()).keys()];
      const uncommonBonds = [...Array(mol.get_num_bonds()).keys()];
      if (mcsMol) {
        const matchedAtomsAndBondsList: ISubstruct[] = JSON.parse(mol!.get_substruct_matches(mcsMol!));
        let errorRateTotal = Number.MAX_SAFE_INTEGER;
        for (let i = 0; i < matchedAtomsAndBondsList.length; i++) {
          let uncommonAtomsSingle = [...uncommonAtoms];
          let uncommonBondsSingle = [...uncommonBonds];
          if (matchedAtomsAndBondsList[i].atoms)
            uncommonAtomsSingle = getArraysDifference(uncommonAtomsSingle, matchedAtomsAndBondsList[i].atoms!.sort((a, b) => {return a - b;}));
          if (matchedAtomsAndBondsList[i].bonds)
            uncommonBondsSingle = getArraysDifference(uncommonBondsSingle, matchedAtomsAndBondsList[i].bonds!.sort((a, b) => {return a - b;}));

          let errorRate = 0;
          for (let j = 0; j < uncommonBondsSingle.length - 1; j++) {
            if (Math.abs(uncommonBondsSingle[j] - uncommonBondsSingle[j + 1]) > 1)
              errorRate++;
          }
          if (errorRate < errorRateTotal) {
            errorRateTotal = errorRate;
            uncommonAtomsBest = uncommonAtomsSingle;
            uncommonBondsBest = uncommonBondsSingle;
          }
        }
      }
      if (rgb) {
        for (let i = 0; i < uncommonAtomsBest.length; i++)
          highlightAtomColors[uncommonAtomsBest[i].toString()] = rgb;
      }
      if (rgb) {
        for (let i = 0; i < uncommonBondsBest.length; i++)
          highlightBondColors[uncommonBondsBest[i].toString()] = rgb;
      }
      return {atoms: uncommonAtomsBest, bonds: uncommonBondsBest, highlightAtomColors, highlightBondColors};
    }
  } catch {
    //do nothing, if molecule is invalid - return empty substruct obj
  } finally {
    mol?.delete();
  }
  return null;
}

/**
 * {@link getUncommonAtomsAndBonds} from a pairwise MCES's own matches (one molecule's side of RDKit Rascal's
 * `atomMatches` and `bondMatches`): the atoms and bonds of `mol` outside them.
 */
export function getUncommonAtomsAndBondsOfMatches(mol: RDMol, atoms: number[], bonds: number[]): ISubstruct {
  const [matchedAtoms, matchedBonds] = [new Set(atoms), new Set(bonds)];
  return {atoms: [...Array(mol.get_num_atoms()).keys()].filter((atom) => !matchedAtoms.has(atom)),
    bonds: [...Array(mol.get_num_bonds()).keys()].filter((bond) => !matchedBonds.has(bond))};
}

function getArraysDifference(largerArr: number[], smallerArr: number[]) {
  const diffArr: number[] = [];
  let counter = 0;
  for (let i = 0; i < largerArr.length; i++)
    largerArr[i] === smallerArr[counter] ? counter++ : diffArr.push(largerArr[i]);
  return diffArr;
}
