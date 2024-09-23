import {monomerSeqToMolfile} from './to-atomic-level-utils';
import {MolfileWithMap, SeqToMolfileWorkerData, SeqToMolfileWorkerRes} from './types';

onmessage = (event) => {
  const {seqList, monomersDict, alphabet, polymerType, start, end}: SeqToMolfileWorkerData = event.data;
  const resMolList: MolfileWithMap[] = new Array<MolfileWithMap>(end - start);
  const molfileWarningList = new Array<string>(0);
  for (let rowI = start; rowI < end; ++rowI) {
    try {
      const seq = seqList[rowI];
      resMolList[rowI - start] = monomerSeqToMolfile(seq, monomersDict, alphabet, polymerType);
    } catch (err: any) {
      const errMsg: string = err instanceof Error ? err.message : err.toString();
      const msg: string = `Cannot get molfile of row #${rowI}: ${errMsg}.`;
      molfileWarningList.push(msg);
    }
  }
  postMessage({molfiles: resMolList, warnings: molfileWarningList} as SeqToMolfileWorkerRes);
};
