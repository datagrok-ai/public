import {monomerSeqToMolfile} from './to-atomic-level-utils';
onmessage = (event) => {
  const {monomerSequencesArray, monomersDict, alphabet, polymerType, start, end} = event.data;
  const molfileList: string[] = new Array<string>(end - start);
  const molfileWarningList = new Array<string>(0);
  for (let rowI = start; rowI < end; ++rowI) {
    try {
      const monomerSeq = monomerSequencesArray[rowI];
      molfileList[rowI - start] = monomerSeqToMolfile(monomerSeq, monomersDict, alphabet, polymerType);
    } catch (err: any) {
      const errMsg: string = err instanceof Error ? err.message : err.toString();
      const msg: string = `Cannot get molfile of row #${rowI}: ${errMsg}.`;
      molfileWarningList.push(msg);
    }
  }
  postMessage({molfileList, molfileWarningList});
};
