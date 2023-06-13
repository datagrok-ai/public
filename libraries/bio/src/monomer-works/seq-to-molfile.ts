import {HELM_POLYMER_TYPE} from '../utils/const';
import {ALPHABET} from '../utils/macromolecule';
import {MolGraph} from './types';

export type SeqToMolfileResult = {
    molfileWarningList: string[]
    molfileList: string[]
}

export async function seqToMolFileWorker(monomerSequencesArray: string[][], monomersDict: Map<string, MolGraph>,
  alphabet: ALPHABET, polymerType: HELM_POLYMER_TYPE, srcColLength: number
): Promise<SeqToMolfileResult> {
  const worker = new Worker(new URL('./seq-to-molfile-worker', import.meta.url));
  worker.postMessage({monomerSequencesArray, monomersDict, alphabet, polymerType, srcColLength});
  return new Promise((resolve) => {
    worker.onmessage = ({data: {molfileList, molfileWarningList}}) => {
      resolve({molfileList, molfileWarningList});
      setTimeout(() => worker.terminate(), 0);
    };
  });
}
