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
  const threadCount = Math.max(navigator.hardwareConcurrency - 2, 1);
  const workers = new Array(threadCount).fill(null)
    .map(() => new Worker(new URL('./seq-to-molfile-worker', import.meta.url)));
  const chunkSize = srcColLength / threadCount;
  let res: string[] = [];
  let warnings: string[] = [];
  const promises = new Array<Promise<SeqToMolfileResult>>(threadCount);
  for (let i = 0; i < threadCount; i++) {
    const start = Math.floor(i * chunkSize);
    const end = (i === threadCount - 1) ? srcColLength : Math.floor((i + 1) * chunkSize);
    workers[i].postMessage({monomerSequencesArray, monomersDict, alphabet, polymerType, start, end});
    promises[i] = new Promise<SeqToMolfileResult>((resolveWorker) => {
      workers[i].onmessage = ({data: {molfileList, molfileWarningList}}): void => {
        resolveWorker({molfileList, molfileWarningList});
      };
    });
  }
  const resArray = await Promise.all(promises);

  resArray.forEach((resItem) => {
    res = res.concat(...resItem.molfileList);
    warnings = warnings.concat(...resItem.molfileWarningList);
  });

  setTimeout(() => {
    workers.forEach((worker) => {
      worker.terminate();
    });
  }, 0);
  return {molfileList: res, molfileWarningList: warnings};
}
