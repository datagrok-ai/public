import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {ALIGNMENT, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {expect} from '@datagrok-libraries/utils/src/test';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';

export type DetectorTestData = { [testName: string]: { csv: string, neg?: string[], pos?: { [colName: string]: PosCol } } };

export type DfReaderFunc = () => Promise<DG.DataFrame>;

export class PosCol {
  constructor(
    public readonly units: NOTATION,
    public readonly aligned: ALIGNMENT | null,
    public readonly alphabet: string | null,
    public readonly alphabetSize: number,
    public readonly alphabetIsMultichar?: boolean,
    public readonly separator?: string,
  ) { };
}

export async function _testNeg(readDf: DfReaderFunc, colName: string) {
  const df: DG.DataFrame = await readDf();
  const col: DG.Column = df.getCol(colName)!;
  const semType: string = await grok.functions
    .call('Bio:detectMacromolecule', {col: col}) as unknown as string;
  if (semType)
    col.semType = semType;

  if (col.semType === DG.SEMTYPE.MACROMOLECULE) {
    const msg = `Negative test detected semType='${col.semType}', units='${col.meta.units}'.`;
    throw new Error(msg);
  }
}

export async function _testPos(
  readDf: DfReaderFunc, colName: string, units: string, aligned: string | null,
  alphabet: string | null, alphabetSize: number, alphabetIsMultichar?: boolean,
  separator: string | null = null,
) {
  const df: DG.DataFrame = await readDf();
  const col: DG.Column = df.col(colName)!;
  const semType: string = await grok.functions
    .call('Bio:detectMacromolecule', {col: col}) as unknown as string;
  if (semType)
    col.semType = semType;

  expect(col.semType, DG.SEMTYPE.MACROMOLECULE);
  expect(col.meta.units, units);
  expect(col.getTag(bioTAGS.aligned), aligned);
  expect(col.getTag(bioTAGS.alphabet), alphabet);
  if (separator)
    expect(col.getTag(bioTAGS.separator), separator);

  const sh = SeqHandler.forColumn(col);
  expect(sh.getAlphabetSize(), alphabetSize);
  expect(sh.getAlphabetIsMultichar(), alphabetIsMultichar);
  if (!sh.isHelm()) {
    expect(sh.aligned, aligned);
    expect(sh.alphabet, alphabet);
  }
}
