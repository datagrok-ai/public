import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';
import {ALPHABET, getAlphabet, NOTATION, UnitsHandler} from '@datagrok-libraries/bio';
import {Column} from 'datagrok-api/dg';

category('detectorsBenchmark', () => {

  let detectFunc: DG.Func;

  before(async () => {
    const funcList: DG.Func[] = DG.Func.find({package: 'Bio', name: 'detectMacromolecule'});
    detectFunc = funcList[0];

    // warm up the detector function
    const col: DG.Column = DG.Column.fromStrings('seq', ['ACGT', 'ACGT', 'ACGT']);
    await detectFunc.prepare({col: col}).call();
  });

  // -- fasta --

  test('fastaDnaShorts50Few50', async () => {
      const et: number = await detectMacromoleculeBenchmark(10, NOTATION.FASTA, ALPHABET.DNA, 50, 50);
    },
    {skipReason: '#1192'});

  test('fastaDnaShorts50Many1E6', async () => {
      const et: number = await detectMacromoleculeBenchmark(10, NOTATION.FASTA, ALPHABET.DNA, 50, 1E6);
    },
    {skipReason: '#1192'});

  test('fastaDnaLong1e6Few50', async () => {
      const et: number = await detectMacromoleculeBenchmark(10, NOTATION.FASTA, ALPHABET.DNA, 1E6, 50);
    },
    {skipReason: '#1192'});

  // -- separator --

  test('separatorDnaShorts50Few50', async () => {
    const et: number = await detectMacromoleculeBenchmark(10, NOTATION.SEPARATOR, ALPHABET.DNA, 50, 50, '/');
  });

  test('separatorDnaShorts50Many1E6', async () => {
      const et: number = await detectMacromoleculeBenchmark(10, NOTATION.SEPARATOR, ALPHABET.DNA, 50, 1E6, '/');
    },
    { /* skipReason: 'slow transmit large dataset to detector' */});

  test('separatorDnaLong1e6Few50', async () => {
      const et: number = await detectMacromoleculeBenchmark(10, NOTATION.SEPARATOR, ALPHABET.DNA, 1E6, 50, '/');
    },
    {skipReason: '#1192'});

  async function detectMacromoleculeBenchmark(
    maxET: number, notation: NOTATION, alphabet: ALPHABET, length: number, count: number, separator?: string
  ): Promise<number> {
    return await benchmark<DG.FuncCall, DG.Column>(10,
      (): DG.FuncCall => {
        const col: DG.Column = generate(notation, [...getAlphabet(alphabet)], length, count, separator);
        const funcCall: DG.FuncCall = detectFunc.prepare({col: col});
        return funcCall;
      },
      async (funcCall: DG.FuncCall): Promise<DG.Column> => {
        return testDetector(funcCall);
      },
      (col: DG.Column) => {
        checkDetectorRes(col, {
          semType: DG.SEMTYPE.MACROMOLECULE,
          notation: notation,
          alphabet: alphabet,
          separator: separator
        });
      });
  }

  function generate(
    notation: NOTATION, alphabet: string[], length: number, count: number, separator?: string
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
      for (let j = 0; j < length; j++) {
        seqMList[j] = alphabet[Math.floor(Math.random() * alphabet.length)];
      }
      return seqMerger(seqMList, separator);
    };

    const seqList: string[] = Array(count);
    for (let i = 0; i < count; i++) {
      seqList[i] = buildSeq(alphabet, length);
    }

    return DG.Column.fromStrings('seq', seqList);
  }

  type TgtType = { semType: string, notation: NOTATION, alphabet: ALPHABET, separator?: string };

  function testDetector(funcCall: DG.FuncCall): DG.Column {
    //const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: col});
    funcCall.callSync();
    const semType = funcCall.getOutputParamValue();

    const col: DG.Column = funcCall.inputs.col;
    if (semType) col.semType = semType;
    return col;
  }

  function checkDetectorRes(col: DG.Column, tgt: TgtType): void {
    const uh = new UnitsHandler(col);
    expect(col.semType, tgt.semType);
    expect(uh.notation, tgt.notation);
    expect(uh.alphabet, tgt.alphabet);
    expect(uh.separator, tgt.separator);
  }
});


/** Returns ET [ms] of test() */
async function benchmark<TData, TRes>(
  maxET: number, prepare: () => TData, test: (data: TData) => Promise<TRes>, check: (res: TRes) => void
): Promise<number> {
  const data: TData = prepare();

  const t1: number = Date.now();
  // console.profile();
  const res: TRes = await test(data);
  //console.profileEnd();
  const t2: number = Date.now();

  check(res);

  const resET: number = t2 - t1;
  if (resET > maxET) {
    const errMsg = `ET ${resET} ms is more than max allowed ${maxET} ms.`;
    console.error(errMsg);
    throw new Error(errMsg);
  } else {
    console.log(`ET ${resET} ms is OK.`);
  }

  return resET;
}
