import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {HELM_POLYMER_TYPE} from '../utils/const';
import {ALPHABET} from '../utils/macromolecule';
import {
  MolfileWithMap, MolGraph, MonomerMap, SeqToMolfileWorkerData, SeqToMolfileWorkerRes
} from './types';
import {ToAtomicLevelRes} from '../utils/seq-helper';
import {getMolColName, getMolHighlightColName} from './utils';
import {SeqHandler} from '../utils/seq-handler';
import {IMonomerLib} from '../types';
import {PolymerType} from '../helm/types';
import {buildMonomerHoverLink} from './monomer-hover';


export type SeqToMolfileResult = {
  mols: {
    molfile: string,
    monomers: MonomerMap
  }[],
  warnings: string[]
}

export async function seqToMolFileWorker(seqCol: DG.Column<string>, monomersDict: Map<string, MolGraph>,
  alphabet: ALPHABET, polymerType: PolymerType, monomerLib: IMonomerLib, rdKitModule: RDModule
): Promise<ToAtomicLevelRes> {
  const srcColLength = seqCol.length;
  const df: DG.DataFrame | undefined = seqCol.dataFrame;
  const threadCount = Math.max(navigator.hardwareConcurrency - 2, 1);
  const workers = new Array(threadCount).fill(null)
    .map(() => new Worker(new URL('./seq-to-molfile-worker', import.meta.url)));
  const chunkSize = srcColLength / threadCount;
  const promises = new Array<Promise<SeqToMolfileWorkerRes>>(threadCount);
  const seqSH = SeqHandler.forColumn(seqCol);
  const canonicalSeqList: string[][] = wu.count(0).take(seqCol.length).map((rowIdx) => {
    const seqSS = seqSH.getSplitted(rowIdx);
    return wu.count(0).take(seqSS.length).map((posIdx) => seqSS.getCanonical(posIdx)).toArray();
  }).toArray();
  for (let i = 0; i < threadCount; i++) {
    const worker = workers[i];
    const start = Math.floor(i * chunkSize);
    const end = (i === threadCount - 1) ? srcColLength : Math.floor((i + 1) * chunkSize);
    promises[i] = new Promise<SeqToMolfileWorkerRes>((resolve) => {
      worker.onmessage = (res: { data: SeqToMolfileWorkerRes }): void => {
        resolve(res.data);
      };
    });
    worker.postMessage({canonicalSeqList, monomersDict, alphabet, polymerType, start, end} as SeqToMolfileWorkerData);
  }

  let molList: MolfileWithMap[] = [];
  let warnings: string[] = [];
  await Promise.all(promises).then((resArray: SeqToMolfileWorkerRes[]) => {
    resArray.forEach((resItem) => {
      molList = molList.concat(...resItem.molfiles);
      warnings = warnings.concat(...resItem.warnings);
    });
  });

  setTimeout(() => {
    workers.forEach((worker) => {
      worker.terminate();
    });
  }, 0);
  const molColName = getMolColName(df, seqCol.name);
  const molHlColName = getMolHighlightColName(df, molColName);
  const molCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, molColName, molList.length)
    .init((rowIdx) => molList[rowIdx].molfile);
  const molHlCol = DG.Column.fromType(DG.COLUMN_TYPE.OBJECT, molHlColName)
    .init((rowIdx) => null); // Highlight nothing
  molCol.semType = DG.SEMTYPE.MOLECULE;
  molCol.meta.units = DG.UNITS.Molecule.MOLBLOCK;

  const monomerHoverHandler = buildMonomerHoverLink(
    seqCol, molCol, monomerLib, rdKitModule);

  return {
    mol: {
      col: molCol, highlightCol: molHlCol,
      monomerHoverLink: monomerHoverHandler
    }, warnings: warnings
  };
}
