import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';


import {after, before, category, test, expect, expectArray} from '@datagrok-libraries/utils/src/test';

import {checkInputColumn, multipleSequenceAlignmentAny} from '../package';
import {UNITS} from 'datagrok-api/dg';


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
    col.setTag(DG.TAGS.UNITS, bio.NOTATION.FASTA);
    col.setTag(bio.TAGS.alphabet, bio.ALPHABET.DNA);
    col.setTag(bio.TAGS.aligned, 'SEQ');

    const [res, msg]: [boolean, string] = checkInputColumn(
      col, 'Test', [bio.NOTATION.FASTA,],
      [bio.ALPHABET.DNA, bio.ALPHABET.RNA, bio.ALPHABET.PT]);

    expect(res, true);
  });

  test('testMsaNegHelm', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const col: DG.Column = df.getCol('seq');
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.setTag(DG.TAGS.UNITS, 'helm');
    col.setTag(bio.TAGS.alphabetSize, '11');
    col.setTag(bio.TAGS.alphabetIsMultichar, 'true');

    const [res, msg]: [boolean, string] = checkInputColumn(
      col, 'Test', [bio.NOTATION.FASTA,],
      [bio.ALPHABET.DNA, bio.ALPHABET.RNA, bio.ALPHABET.PT]);

    expect(res, false);
  });

  test('testMsaNegUN', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const col: DG.Column = df.getCol('seq');
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.setTag(DG.TAGS.UNITS, 'fasta');
    col.setTag(bio.TAGS.alphabet, 'UN');
    col.setTag(bio.TAGS.alphabetSize, '11');
    col.setTag(bio.TAGS.alphabetIsMultichar, 'true');
    col.setTag(bio.TAGS.aligned, 'SEQ');

    const [res, msg]: [boolean, string] = checkInputColumn(
      col, 'Test', [bio.NOTATION.FASTA,],
      [bio.ALPHABET.DNA, bio.ALPHABET.RNA, bio.ALPHABET.PT]);

    expect(res, false);
  });

  test('testGetActionFunctionMeta', async () => {
    const func: DG.Func = DG.Func.find({package: 'Bio', name: 'multipleSequenceAlignmentAny'})[0];
    const sequenceInput: DG.Property = func.inputs.find((i) => i.name == 'sequence')!;
    let k = 11;
  });
});