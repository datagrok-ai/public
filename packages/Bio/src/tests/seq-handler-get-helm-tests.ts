import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {ISeqHelper, getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

category('SeqHandler: getHelm', () => {
  let seqHelper: ISeqHelper;
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings; // backup

  before(async () => {
    seqHelper = await getSeqHelper();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    // Test 'helm' requires default monomer library loaded
    await monomerLibHelper.loadMonomerLibForTests(); // load default libraries
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true); // load user settings libraries
  });

  const tests: {
    [testName: string]: {
      src: {
        seq: string, notation: NOTATION, separator?: string
      }, tgt: { helm: string }
    }
  } = {
    'fasta': {
      src: {seq: 'MDYKETMDYKET', notation: NOTATION.FASTA},
      tgt: {helm: 'PEPTIDE1{M.D.Y.K.E.T.M.D.Y.K.E.T}$$$$'},
    },
    'separator': {
      src: {seq: 'M-D-Y-K-E-T-M-D-Y-K-E-T', notation: NOTATION.SEPARATOR, separator: '-'},
      tgt: {helm: 'PEPTIDE1{M.D.Y.K.E.T.M.D.Y.K.E.T}$$$$'},
    },
    'helm': {
      src: {seq: 'PEPTIDE1{M.D.Y.K.E.T}$$$$', notation: NOTATION.HELM},
      tgt: {helm: 'PEPTIDE1{M.D.Y.K.E.T}$$$$'},
    },
    'helm-cyclic': {
      src: {seq: 'PEPTIDE1{M.D.Y.K.E.T}$PEPTIDE1,PEPTIDE1,6:R2-1:R1$$$V2.0', notation: NOTATION.HELM},
      tgt: {helm: 'PEPTIDE1{M.D.Y.K.E.T}$PEPTIDE1,PEPTIDE1,6:R2-1:R1$$$V2.0'}
    },

    // TODO: Add tests for cyclized and dimerized
    // 'separator-cyclized': {
    //   src: {seq: 'R-F-C(1)-T-G-H-F-Y-P-C(1)-meI', notation: NOTATION.SEPARATOR, separator: '-',},
    //   tgt: {helm: 'PEPTIDE1{[R].[F].[C].[T].[G].[H].[F].[Y].[P].[C].[meI]}$PEPTIDE1,PEPTIDE1,3:R3-10:R3$$$'}
    // },
    // 'separator-dimerized': {
    //   src: {seq: '(#2)C-{R-F-C(2)-T-G-H-F-Y-P-C(2)-Mei}', notation: NOTATION.SEPARATOR, separator: '-',},
    //   tgt: {helm: ''}
    // },
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(testName, async () => {
      await _testSeqHandlerGetHelm(testData.src.seq, testData.src.notation, testData.src.separator, testData.tgt.helm);
    });
  }

  async function _testSeqHandlerGetHelm(
    srcSeq: string, srcNotation: NOTATION, srcSeparator: string | undefined, tgtHelm: string
  ): Promise<void> {
    const seqCol = DG.Column.fromStrings('seq', [srcSeq]);
    const df = DG.DataFrame.fromColumns([seqCol]);
    // seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    // seqCol.setTag(DG.TAGS.UNITS, srcNotation);
    // seqCol.setTag(TAGS.alphabet, 'PT');
    // if (srcSeparator) seqCol.setTag(TAGS.separator, srcSeparator);
    await grok.data.detectSemanticTypes(df);

    const sh = seqHelper.getSeqHandler(seqCol);
    const resSeqValue = await sh.getValue(0);
    const resHelm = resSeqValue.helm;
    expect(resHelm, tgtHelm);
  }
});
