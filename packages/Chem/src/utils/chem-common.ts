import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {getMolSafe} from './mol-creation_rdkit';
import {MolList, RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import { ISubstruct } from '@datagrok-libraries/chem-meta/src/types';

export const defaultMorganFpRadius = 2;
export const defaultMorganFpLength = 2048;

const lockPromiseForKey: any = {};
const unlockFunctionForKey: any = {};

export enum Fingerprint {
  Morgan = 'Morgan',
  RDKit = 'RDKit',
  Pattern = 'Pattern',
  //Avalon = 'Avalon',
  MACCS = 'MACCS',
  AtomPair = 'AtomPair',
  TopologicalTorsion = 'TopologicalTorsion'
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
  return fp ? BitArray.fromBytes(fp) : null;
}

export function isMolBlock(s: string): boolean {
  return s.includes('M  END');
}

/** @deprecated Use DG.Color.hexToPercentRgb */
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

export function stringArrayToMolList(molecules: string[], rdkit: RDModule): MolList {
  const mols = new rdkit.MolList();
  for (let i = 0; i < molecules.length; i++)
    appendMolList(molecules[i], mols, rdkit);
  return mols;
}

export function appendMolList(molString: string | null, mols: MolList, rdkit: RDModule) {
  let molSafe;
  try {
    molSafe = getMolSafe(molString ?? '', {}, rdkit);
    if (molSafe.mol === null || molSafe.isQMol || !molSafe.kekulize) {
      molSafe.mol?.delete();
      molSafe = getMolSafe('', {}, rdkit);
    }
    mols.append(molSafe.mol!);
  } finally {
    molSafe?.mol?.delete();
  }
}

export function getSigFigs(n: number, sig: number) {
  var mult = Math.pow(10, sig - Math.floor(Math.log(n) / Math.LN10) - 1);
  return Math.round(n * mult) / mult;
}


export function getFirstNSymbols(number: number, digits: number): string {
  const str = number.toFixed(digits).toString().slice(0, digits);
  return str[str.length - 1] === '.' ? number.toFixed(1).toString() : str;
}

export function getUncommonAtomsAndBonds(molecule: string, mcsMol: RDMol | null, rdkit: RDModule,
  col?: string, addHs?: boolean): ISubstruct | null {
  let mol = getMolSafe(molecule, {}, rdkit).mol;
  const substruct: ISubstruct | null = null;
  let uncommonAtomsBest: number[] = [];
  let uncommonBondsBest: number[] = [];
  const highlightAtomColors: { [key: string]: number[] } = {};
  const highlightBondColors: { [key: string]: number[] } = {};
  try {
    if (mol) {
      if (addHs) {
        const molblockWithHs = mol.add_hs();
        mol = getMolSafe(molblockWithHs, {}, rdkit).mol;
      }
      if (mol) {
        let rdKol: number[] | null = null;
        if (col)
          rdKol = hexToPercentRgb(col);
        const uncommonAtoms = [...Array(mol.get_num_atoms()).keys()];
        const uncommonBonds = [...Array(mol.get_num_bonds()).keys()];
        if (mcsMol) {
          const matchedAtomsAndBondsList: ISubstruct[] = JSON.parse(mol!.get_substruct_matches(mcsMol!));
          let errorRateTotal = Number.MAX_SAFE_INTEGER;
          for (let i = 0; i < matchedAtomsAndBondsList.length; i++) {
            let uncommonAtomsSingle = [...uncommonAtoms];
            let uncommonBondsSingle = [...uncommonBonds];
            if (matchedAtomsAndBondsList[i].atoms)
              uncommonAtomsSingle = getArraysDifference(uncommonAtomsSingle, matchedAtomsAndBondsList[i].atoms!.sort((a, b) => { return a - b; }));
            if (matchedAtomsAndBondsList[i].bonds)
              uncommonBondsSingle = getArraysDifference(uncommonBondsSingle, matchedAtomsAndBondsList[i].bonds!.sort((a, b) => { return a - b; }));

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
        if (rdKol) {
          for (let i = 0; i < uncommonAtomsBest.length; i++)
            highlightAtomColors[uncommonAtomsBest[i].toString()] = rdKol;
        }
        if (rdKol) {
          for (let i = 0; i < uncommonBondsBest.length; i++)
            highlightBondColors[uncommonBondsBest[i].toString()] = rdKol;
        }
        return {atoms: uncommonAtomsBest, bonds: uncommonBondsBest, highlightAtomColors, highlightBondColors};
      }
    }
  } catch {
    //do nothing, if molecule is invalid - return empty substruct obj
  } finally {
    mol?.delete();
  }
  return substruct;
}

function getArraysDifference(largerArr: number[], smallerArr: number[]) {
  const diffArr: number[] = [];
  let counter = 0;
  for (let i = 0; i < largerArr.length; i++)
    largerArr[i] === smallerArr[counter] ? counter++ : diffArr.push(largerArr[i]);
  return diffArr;
}
