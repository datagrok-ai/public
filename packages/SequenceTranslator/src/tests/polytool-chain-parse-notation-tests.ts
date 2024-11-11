import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, after, category, expect, test, expectArray, testEvent, delay} from '@datagrok-libraries/utils/src/test';
import {Chain} from '../polytool/conversion/pt-chain';
import {getInnerIdx, getOuterIdx} from '../polytool/conversion/pt-misc';
import {getRules} from '../polytool/conversion/pt-rules';
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
        helm: 'PEPTIDE1{R.F.[C(1)].T.G.H.F.Y.P.[C(1)].[meI]}$$$$V2.0',
      }
    },
    'dimerized1': {
      src: {seq: '(#3)Succ-{A(CHOL)-F-C(1)-T-G-H-Y-P-C(1)-NH2}'},
      tgt: {
        // TODO: Target test data requires clarification
        //monomerCount: [2, 8], linkageCount: 0,
        helm: 'PEPTIDE1{[(#3)Succ]}' + '|' +
          'PEPTIDE2{[A(CHOL)].F.[C(1)].T.G.H.Y.P.[C(1)].[NH2]}' + '$' +
          'PEPTIDE1,PEPTIDE2,1:R1-1:R1' + '$$$' + 'V2.0',
      }
    },
    'dimerized2': {
      src: {seq: '($3)Succ-{R-F-C(1)-T-G-H-F-P-C(1)-NH2}($2){A(CHOL)-F-C(1)-T-G-H-F-P-C(1)-NH2}'},
      tgt: {
        // TODO: Target test data requires clarification
        //monomerCount: [2, 8], linkageCount: 0,
        helm: 'PEPTIDE1{[($3)Succ]}' + '|' +
          'PEPTIDE2{R.F.[C(1)].T.G.H.F.P.[C(1)].[NH2]}' + '|' +
          'PEPTIDE3{[($2)A(CHOL)].F.[C(1)].T.G.H.F.P.[C(1)].[NH2]}' + '$' +
          'PEPTIDE1,PEPTIDE2,1:R1-1:R1' + '$$$' + 'V2.0',
      }
    }
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}`, async () => {
      const rules = await getRules(['rules_example.json']);
      const resChain = await Chain.parseNotation(testData.src.seq, helmHelper);
      //expectArray(resChain.monomers.map((mL) => mL.length), testData.tgt.monomerCount);
      //expect(resChain.linkages.length, testData.tgt.linkageCount);
      // expect(resChain.getNotationHelm(), testData.tgt.helm);
      // expect(resChain.getNotation(), testData.src.seq);

      const resMol = resChain.mol!;
      const hwe = helmHelper.createHelmWebEditor();
      hwe.editor.setMol(resMol!);
      const resMolHelm = hwe.editor.getHelm();

      const resHelm = resChain.getNotationHelm();

      expect(resMolHelm, testData.tgt.helm);
      expect(resHelm, testData.tgt.helm);
    }, testName == 'reaction2' ? {skipReason: 'reverse reaction'} : undefined);
  }

  const innerIdxTests = {
    '0-in-4-4-3': {inIdx: 0, spIdx: 0, outIdx: 0, monomers: [['0', '1', '2', '3'], ['4', '5', '6', '7'], ['8', '9', '10']]},
    '3-in-4-4-3': {inIdx: 3, spIdx: 0, outIdx: 3, monomers: [['0', '1', '2', '3'], ['4', '5', '6', '7'], ['8', '9', '10']]},
    '4-in-4-4-3': {inIdx: 0, spIdx: 1, outIdx: 4, monomers: [['0', '1', '2', '3'], ['4', '5', '6', '7'], ['8', '9', '10']]},
    '7-in-4-4-3': {inIdx: 3, spIdx: 1, outIdx: 7, monomers: [['0', '1', '2', '3'], ['4', '5', '6', '7'], ['8', '9', '10']]},
    '8-in-4-4-3': {inIdx: 0, spIdx: 2, outIdx: 8, monomers: [['0', '1', '2', '3'], ['4', '5', '6', '7'], ['8', '9', '10']]},
    // '11-in-4-4-3': {inIdx: 0, spIdx: 3, outIdx: 11, monomers: [['0', '1', '2', '3'], ['4', '5', '6', '7'], ['8', '9', '10']]},
    // '12-in-4-4-3': {inIdx: 1, spIdx: 3, outIdx: 12, monomers: [['0', '1', '2', '3'], ['4', '5', '6', '7'], ['8', '9', '10']]},
    '0-in-1-1-6-3': {inIdx: 0, spIdx: 0, outIdx: 0, monomers: [['0'], ['1'], ['2', '3', '4', '5', '6', '7'], ['8', '9', '10']]},
    '1-in-1-1-6-3': {inIdx: 0, spIdx: 1, outIdx: 1, monomers: [['0'], ['1'], ['2', '3', '4', '5', '6', '7'], ['8', '9', '10']]},
    '2-in-1-1-6-3': {inIdx: 0, spIdx: 2, outIdx: 2, monomers: [['0'], ['1'], ['2', '3', '4', '5', '6', '7'], ['8', '9', '10']]},
    '3-in-1-1-6-3': {inIdx: 1, spIdx: 2, outIdx: 3, monomers: [['0'], ['1'], ['2', '3', '4', '5', '6', '7'], ['8', '9', '10']]},
    '7-in-1-1-6-3': {inIdx: 5, spIdx: 2, outIdx: 7, monomers: [['0'], ['1'], ['2', '3', '4', '5', '6', '7'], ['8', '9', '10']]},
    '8-in-1-1-6-3': {inIdx: 0, spIdx: 3, outIdx: 8, monomers: [['0'], ['1'], ['2', '3', '4', '5', '6', '7'], ['8', '9', '10']]},
//    '11-in-1-1-6-3': {inIdx: 0, spIdx: 4}, src: {outIdx: 11, monomers: [['0'], ['1'], ['2', '3', '4', '5', '6', '7'], ['8', '9', '10']]}},
  };

  for (const [testName, {inIdx, spIdx, outIdx, monomers}] of Object.entries(innerIdxTests)) {
    test(`innerIdx-${testName}`, async () => {
      const [resInIdx, resSpIdx] = getInnerIdx(outIdx, monomers);
      expect(resInIdx, inIdx);
      expect(resSpIdx, spIdx);
    });
  }

  for (const [testName, {inIdx, spIdx, outIdx, monomers}] of Object.entries(innerIdxTests)) {
    test(`outerIdx-${testName}`, async () => {
      const resOutIdx = getOuterIdx(inIdx, spIdx, monomers);
      expect(resOutIdx, outIdx);
    });
  }
});
