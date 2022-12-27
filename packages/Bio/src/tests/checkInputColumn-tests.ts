import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import {after, before, category, test, expect, expectArray} from '@datagrok-libraries/utils/src/test';

import {checkInputColumn, multipleSequenceAlignmentAny} from '../package';
import {ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio';

category('checkInputColumn', () => {
  const csv = `seq
seq1,
seq2,
seq3,
seq4`;

  test('testMsaPos', async () => {
    const func: DG.Func = DG.Func.find({package: 'Bio', name: 'multipleSequenceAlignmentAny'})[0];
    const funcInputColumnProperty: DG.Property = func.inputs.find((i) => i.name == 'sequence')!;

    const k = 11;

    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const col: DG.Column = df.getCol('seq');
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.setTag(DG.TAGS.UNITS, NOTATION.FASTA);
    col.setTag(bioTAGS.alphabet, ALPHABET.DNA);
    col.setTag(bioTAGS.aligned, 'SEQ');

    const [res, msg]: [boolean, string] = checkInputColumn(
      col, 'Test', [NOTATION.FASTA],
      [ALPHABET.DNA, ALPHABET.RNA, ALPHABET.PT]);

    expect(res, true);
  });

  test('testMsaNegHelm', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const col: DG.Column = df.getCol('seq');
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    // col.setTag(bio.TAGS.alphabetSize, '11');
    col.setTag(bioTAGS.alphabetIsMultichar, 'true');

    const [res, msg]: [boolean, string] = checkInputColumn(
      col, 'Test', [NOTATION.FASTA],
      [ALPHABET.DNA, ALPHABET.RNA, ALPHABET.PT]);

    expect(res, false);
  });

  test('testMsaNegUN', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const col: DG.Column = df.getCol('seq');
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.setTag(DG.TAGS.UNITS, NOTATION.FASTA);
    col.setTag(bioTAGS.alphabet, 'UN');
    col.setTag(bioTAGS.alphabetSize, '11');
    col.setTag(bioTAGS.alphabetIsMultichar, 'true');
    col.setTag(bioTAGS.aligned, 'SEQ');

    const [res, msg]: [boolean, string] = checkInputColumn(
      col, 'Test', [NOTATION.FASTA],
      [ALPHABET.DNA, ALPHABET.RNA, ALPHABET.PT]);

    expect(res, false);
  });

  test('testGetActionFunctionMeta', async () => {
    const func: DG.Func = DG.Func.find({package: 'Bio', name: 'multipleSequenceAlignmentAny'})[0];
    const sequenceInput: DG.Property = func.inputs.find((i) => i.name == 'sequence')!;
    const k = 11;
  });
});
