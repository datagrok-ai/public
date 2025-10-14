import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {checkInputColumn} from '../utils/check-input-column';

category('checkInputColumn', () => {
  let seqHelper: ISeqHelper;

  before(async () => {
    seqHelper = await getSeqHelper();
  });

  const csv = `seq
seq1,
seq2,
seq3,
seq4`;

  test('testMsaPos', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const col: DG.Column = df.getCol('seq');
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = NOTATION.FASTA;
    col.setTag(bioTAGS.alphabet, ALPHABET.DNA);
    col.setTag(bioTAGS.aligned, 'SEQ');

    const [res, _msg]: [boolean, string] = checkInputColumn(col, 'Test', seqHelper,
      [NOTATION.FASTA], [ALPHABET.DNA, ALPHABET.RNA, ALPHABET.PT]);

    expect(res, true);
  });

  test('testMsaNegHelm', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const col: DG.Column = df.getCol('seq');
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = NOTATION.HELM;
    // col.setTag(bio.TAGS.alphabetSize, '11');
    col.setTag(bioTAGS.alphabetIsMultichar, 'true');

    const [res, _msg]: [boolean, string] = checkInputColumn(col, 'Test', seqHelper,
      [NOTATION.FASTA], [ALPHABET.DNA, ALPHABET.RNA, ALPHABET.PT]);

    expect(res, false);
  });

  test('testMsaNegUN', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const col: DG.Column = df.getCol('seq');
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = NOTATION.FASTA;
    col.setTag(bioTAGS.alphabet, 'UN');
    col.setTag(bioTAGS.alphabetSize, '11');
    col.setTag(bioTAGS.alphabetIsMultichar, 'true');
    col.setTag(bioTAGS.aligned, 'SEQ');

    const [res, _msg]: [boolean, string] = checkInputColumn(col, 'Test', seqHelper,
      [NOTATION.FASTA], [ALPHABET.DNA, ALPHABET.RNA, ALPHABET.PT]);

    expect(res, false);
  });

  test('testGetActionFunctionMeta', async () => {
    const func: DG.Func = DG.Func.find({package: 'Bio', name: 'multipleSequenceAlignmentDialog'})[0];
    const _sequenceInput: DG.Property = func.inputs.find((i) => i.name == 'sequence')!;
  });
});
