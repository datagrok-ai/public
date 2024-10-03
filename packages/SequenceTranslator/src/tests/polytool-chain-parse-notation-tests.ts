import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, after, category, expect, test, expectArray, testEvent, delay} from '@datagrok-libraries/utils/src/test';
import {Chain} from '../polytool/pt-conversion';
import {getRules} from '../polytool/pt-rules';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

category('PolyTool: Chain: parseNotation', () => {
  let helmHelper: IHelmHelper;

  before(async () => {
    helmHelper = await getHelmHelper(); // initialize JSDraw2 and org
  });

  const tests = {
    'cyclized': {
      src: {seq: 'R-F-C(1)-T-G-H-F-Y-P-C(1)-meI'},
      tgt: {
        //monomerCount: [11], linkageCount: 1,
        helm: 'PEPTIDE1{R.F.[C(1)].T.G.H.F.Y.P.[C(1)].[meI]}$$$$',
      },
    },
    'dimerized1': {
      src: {seq: '(#2)Succ-{A(CHOL)-F-C(2)-T-G-H-Y-P-C(2)-NH2}'},
      tgt: {
        // TODO: Target test data requires clarification
        //monomerCount: [2, 8], linkageCount: 0,
        helm: 'PEPTIDE1{[(#2)Succ]}|PEPTIDE2{A(CHOL).F.[C(1)].T.G.H.Y.P.[C(1)].[NH2]}$PEPTIDE1,PEPTIDE2,1:R1-1:R1$$$',
      }
    },
    'dimerized2': {
      src: {seq: '($2)Succ-{R-F-C(1)-T-G-H-F-P-C(1)-NH2}($2){A(CHOL)-F-C(1)-T-G-H-F-P-C(1)-NH2}'},
      tgt: {
        // TODO: Target test data requires clarification
        //monomerCount: [2, 8], linkageCount: 0,
        helm: 'PEPTIDE1{[($2)Succ]}|PEPTIDE2{R.F.[C(1)].T.G.H.F.P.[C(1)].NH2}|PEPTIDE3{[($2)A(CHOL)].F.[C(1)].T.G.H.F.P.[C(1)].[NH2]}$PEPTIDE1,PEPTIDE2,1:R1-1:R1$$$',
      }
    },
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}`, async () => {
      const rules = await getRules(['rules_example.json']);
      const resChain = await Chain.parseNotation(testData.src.seq);
      //expectArray(resChain.monomers.map((mL) => mL.length), testData.tgt.monomerCount);
      //expect(resChain.linkages.length, testData.tgt.linkageCount);
      // expect(resChain.getNotationHelm(), testData.tgt.helm);
      // expect(resChain.getNotation(), testData.src.seq);

      const hwe = helmHelper.createHelmWebEditor();
      hwe.editor.setMol(resChain.mol!);
      const resHelm = hwe.editor.getHelm();
      expect(resHelm, testData.tgt.helm);
    }, testName == 'reaction2' ? {skipReason: 'reverse reaction'} : undefined);
  }
});
