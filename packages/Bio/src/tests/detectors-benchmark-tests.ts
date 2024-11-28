/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, category, test, expect} from '@datagrok-libraries/utils/src/test';
import {ALPHABET, getAlphabet, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqHelper, getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {_package} from '../package-test';


category('detectorsBenchmark', () => {
  let seqHelper: ISeqHelper;
  let detectFunc: DG.Func;

  before(async () => {
    seqHelper = await getSeqHelper();
    const funcList: DG.Func[] = DG.Func.find({package: 'Bio', name: 'detectMacromolecule'});
    detectFunc = funcList[0];

    // warm up the detector function
    const col: DG.Column = DG.Column.fromStrings('seq', ['ACGT', 'ACGT', 'ACGT']);
    await detectFunc.prepare({col: col}).call();
  });

  // -- fasta --

  test('fastaDnaShorts50Few50', async () => {
    await detectMacromoleculeBenchmark(20, NOTATION.FASTA, ALPHABET.DNA, 50, 50);
  });

  test('fastaDnaShorts50Many1E6', async () => {
    await detectMacromoleculeBenchmark(20, NOTATION.FASTA, ALPHABET.DNA, 50, 1E6);
  });

  test('fastaDnaLong1e6Few50', async () => {
    await detectMacromoleculeBenchmark(20, NOTATION.FASTA, ALPHABET.DNA, 1E6, 50);
  });

  // -- separator --

  test('separatorDnaShorts50Few50', async () => {
    await detectMacromoleculeBenchmark(20, NOTATION.SEPARATOR, ALPHABET.DNA, 50, 50, '/');
  });

  test('separatorDnaShorts50Many1E6', async () => {
    await detectMacromoleculeBenchmark(20, NOTATION.SEPARATOR, ALPHABET.DNA, 50, 1E6, '/');
  });

  test('separatorDnaLong1e6Few50', async () => {
    await detectMacromoleculeBenchmark(20, NOTATION.SEPARATOR, ALPHABET.DNA, 1E6, 50, '/');
  });

  async function detectMacromoleculeBenchmark(
    maxET: number, notation: NOTATION, alphabet: ALPHABET, length: number, count: number, separator?: string,
  ): Promise<number> {
    return await benchmark<DG.FuncCall, DG.Column>(maxET,
      /* prepare */ async (): Promise<DG.FuncCall> => {
        const col: DG.Column = generate(notation, [...getAlphabet(alphabet)], length, count, separator);
        const funcCall: DG.FuncCall = detectFunc.prepare({col: col});
        // warm-up Bio
        testDetector(funcCall);
        return funcCall;
      },
      /* test */ (funcCall: DG.FuncCall): DG.Column => { // sync call for stability
        return testDetector(funcCall);
      },
      /* check */ (col: DG.Column) => {
        checkDetectorRes(col, {
          semType: DG.SEMTYPE.MACROMOLECULE,
          notation: notation,
          alphabet: alphabet,
          separator: separator,
        });
      });
  }

  function generate(
    notation: NOTATION, alphabet: string[], length: number, count: number, separator?: string,
  ): DG.Column {
    let seqMerger: (seqMList: string[], separator?: string) => string;

    switch (notation) {
    case NOTATION.FASTA:
      seqMerger = (seqMList: string[]): string => {
        let res: string = '';
        for (let j = 0; j < seqMList.length; j++) {
          const m = seqMList[j];
          res += m.length == 1 ? m : `[${m}]`;
        }
        return res;
      };
      break;
    case NOTATION.SEPARATOR:
      seqMerger = (seqMList: string[], separator?: string): string => {
        return seqMList.join(separator);
      };
      break;
    default:
      throw new Error(`Not supported notation '${notation}'.`);
    }

    const buildSeq = (alphabet: string[], length: number): string => {
      const seqMList = new Array<string>(length);
      for (let j = 0; j < length; j++)
        seqMList[j] = alphabet[Math.floor(Math.random() * alphabet.length)];

      return seqMerger(seqMList, separator);
    };

    const seqList: string[] = Array(count);
    for (let i = 0; i < count; i++)
      seqList[i] = buildSeq(alphabet, length);


    return DG.Column.fromStrings('seq', seqList);
  }

  type TgtType = { semType: string, notation: NOTATION, alphabet: ALPHABET, separator?: string };

  function testDetector(funcCall: DG.FuncCall): DG.Column {
    funcCall.callSync();
    const semType = funcCall.getOutputParamValue() as string;

    const col: DG.Column = funcCall.inputs.col as unknown as DG.Column;
    if (semType) col.semType = semType;
    return col;
  }

  function checkDetectorRes(col: DG.Column, tgt: TgtType): void {
    const sh = seqHelper.getSeqHandler(col);
    expect(col.semType === tgt.semType, true);
    expect(sh.notation === tgt.notation, true);
    expect(sh.alphabet === tgt.alphabet, true);
    expect(sh.separator === tgt.separator, true);
  }
});


//Returns ET [ms] of test()
async function benchmark<TData, TRes>(
  maxET: number, prepare: () => Promise<TData>, test: (data: TData) => TRes, check: (res: TRes) => void,
): Promise<number> {
  const data: TData = await prepare();

  const t1: number = Date.now();
  // console.profile();
  const res: TRes = test(data); // sync call for stability
  //console.profileEnd();
  const t2: number = Date.now();

  check(res);

  const resET = t2 - t1;
  if (resET > maxET) {
    const errMsg = `ET ${resET} ms is more than max allowed ${maxET} ms.`;
    console.error(errMsg);
    throw new Error(errMsg);
  } else
    console.log(`ET ${resET} ms is OK.`);

  return resET;
}
