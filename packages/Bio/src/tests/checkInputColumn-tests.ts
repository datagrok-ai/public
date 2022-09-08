import {after, before, category, test, expect, expectArray} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {checkInputColumn, multipleSequenceAlignmentAny} from '../package';
import {UNITS} from 'datagrok-api/dg';
import {ALPHABET, UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';


category('checkInputColumn', () => {

  const csv = `seq
seq1,
seq2,
seq3,
seq4`;

  test('testMsaPos', async () => {
    const func: DG.Func = DG.Func.find({package: 'Bio', name: 'multipleSequenceAlignmentAny'})[0];
    const funcInputColumnProperty: DG.Property = func.inputs.find((i) => i.name == 'sequence')!;

    let k = 11;

    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const col: DG.Column = df.getCol('seq');
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.setTag(DG.TAGS.UNITS, 'fasta');
    col.setTag(UnitsHandler.TAGS.alphabet, ALPHABET.DNA);

    const [res, msg]: [boolean, string] = checkInputColumn(
      col, 'Test', ['fasta',], ['DNA', 'RNA', 'PT']);

    expect(res, true);
  });

  test('testMsaNegHelm', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const col: DG.Column = df.getCol('seq');
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.setTag(DG.TAGS.UNITS, 'helm');
    col.setTag(UnitsHandler.TAGS.alphabetSize, '11');
    col.setTag(UnitsHandler.TAGS.alphabetIsMultichar, 'true');

    const [res, msg]: [boolean, string] = checkInputColumn(
      col, 'Test', ['fasta',], ['DNA', 'RNA', 'PT']);

    expect(res, false);
  });

  test('testMsaNegUN', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const col: DG.Column = df.getCol('seq');
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.setTag(DG.TAGS.UNITS, 'fasta');
    col.setTag(UnitsHandler.TAGS.alphabet, 'UN');
    col.setTag(UnitsHandler.TAGS.alphabetSize, '11');
    col.setTag(UnitsHandler.TAGS.alphabetIsMultichar, 'true');

    const [res, msg]: [boolean, string] = checkInputColumn(
      col, 'Test', ['fasta',], ['DNA', 'RNA', 'PT']);

    expect(res, false);
  });

  test('testGetActionFunctionMeta', async () => {
    const func: DG.Func = DG.Func.find({package: 'Bio', name: 'multipleSequenceAlignmentAny'})[0];
    const sequenceInput: DG.Property = func.inputs.find((i) => i.name == 'sequence')!;
    let k = 11;
  });
});