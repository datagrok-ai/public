import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, after, category, expect, test, expectArray, testEvent, delay} from '@datagrok-libraries/utils/src/test';
import {Chain} from '../polytool/pt-conversion';
import {getRules} from '../polytool/pt-rules';

category('PolyTool: Chain: fromNotation', () => {
  const tests = {
    'cyclized': {
      src: {seq: 'R-F-C(1)-T-G-H-F-Y-P-C(1)-meI'},
      tgt: {
        monomerCount: [11], linkageCount: 1,
        helm: 'PEPTIDE1{R.F.C.T.G.H.F.Y.P.C.[meI]}$PEPTIDE1,PEPTIDE1,3:R3-10:R3$$$',
      },
    },
    'reaction1': {
      src: {seq: 'R-F-azG(3)-T-G-H-F-Y-P-aG(3)-meI'},
      tgt: {
        monomerCount: [9, 1], linkageCount: 2,
        helm: 'PEPTIDE1{R.F.[GGaz].T.G.H.F.Y.P}|PEPTIDE2{[meI]}$PEPTIDE1,PEPTIDE1,3:R3-9:R2|PEPTIDE1,PEPTIDE2,3:R4-1:R1$$$',
      }
    },
    'reaction2': {
      src: {seq: 'R-F-aG(3)-T-G-H-F-Y-P-azG(3)-meI'},
      tgt: {
        // TODO: Target test data requires clarification
        monomerCount: [2, 8], linkageCount: 0,
        helm: 'PEPTIDE1{R.F}|PEPTIDE2{T.G.H.F.Y.P.[GGaz].[meI]}$PEPTIDE1,PEPTIDE2,2:R2-7:R3|PEPTIDE2,PEPTIDE2,1:R1-7:R4,$$$',
      }
    }

  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}`, async () => {
      const rules = await getRules(['rules_example.json']);
      const resChain = Chain.fromNotation(testData.src.seq, rules);
      expectArray(resChain.monomers.map((mL) => mL.length), testData.tgt.monomerCount);
      expect(resChain.linkages.length, testData.tgt.linkageCount);
      expect(resChain.getHelm(), testData.tgt.helm);
    }, testName == 'reaction2' ? {skipReason: 'reverse reaction'} : undefined);
  }
});
