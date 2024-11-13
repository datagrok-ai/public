import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, after, category, expect, test, expectArray, testEvent, delay} from '@datagrok-libraries/utils/src/test';
import {Chain} from '../polytool/conversion/pt-chain';
import {getRules} from '../polytool/conversion/pt-rules';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

category('PolyTool: Chain', () => {
  let helmHelper: IHelmHelper;

  before(async () => {
    helmHelper = await getHelmHelper();
  });

  const tests = {
    'cyclized': {
      data: {
        templateSeq: 'R-F-C(1)-T-G-H-F-Y-P-C(1)-meI',
        templateHelm: 'PEPTIDE1{R.F.[C(1)].T.G.H.F.Y.P.[C(1)].[meI]}$$$$V2.0',
        mmHelm: 'PEPTIDE1{R.F.C.T.G.H.F.Y.P.C.[meI]}$PEPTIDE1,PEPTIDE1,3:R3-10:R3$$$V2.0',
      },
      tgt: {
        templateChain: {monomerCount: [11], linkageCount: 0},
        mmChain: {monomerCount: [11], linkageCount: 1}
      },
    },
    'reaction1': {
      data: {
        templateSeq: 'R-F-azG(4)-T-G-H-F-Y-P-aG(4)-meI',
        templateHelm: 'PEPTIDE1{R.F.[azG(4)].T.G.H.F.Y.P.[aG(4)].[meI]}$$$$V2.0',
        mmHelm: 'PEPTIDE1{R.F.[GGaz].T.G.H.F.Y.P}|PEPTIDE2{[meI]}$PEPTIDE1,PEPTIDE1,3:R3-9:R2|PEPTIDE1,PEPTIDE2,3:R4-1:R1$$$V2.0',
      },
      tgt: {
        templateChain: {monomerCount: [11], linkageCount: 0},
        mmChain: {monomerCount: [9, 1], linkageCount: 2}
      }
    },
    'reaction2': {
      data: {
        templateSeq: 'R-F-aG(4)-T-G-H-F-Y-P-azG(4)-meI',
        templateHelm: 'PEPTIDE1{R.F.[aG(4)].T.G.H.F.Y.P.[azG(4)].[meI]}$$$$V2.0',
        mmHelm: 'PEPTIDE1{R.F}|PEPTIDE2{T.G.H.F.Y.P.[GGaz].[meI]}$PEPTIDE1,PEPTIDE2,2:R2-7:R3|PEPTIDE2,PEPTIDE2,1:R1-7:R4,$$$V2.0',
      },
      tgt: {
        templateChain: {monomerCount: [11], linkageCount: 0},
        mmChain: {monomerCount: [2, 8], linkageCount: 2}
      }
    },
    'dimerized1': {
      data: {
        templateSeq: '(#3)Succ-{A(CHOL)-F-C(1)-T-G-H-Y-P-C(1)-NH2}',
        templateHelm: 'PEPTIDE1{[(#3)Succ]}' + '|' +
          'PEPTIDE2{[A(CHOL)].F.[C(1)].T.G.H.Y.P.[C(1)].[NH2]}' + '$' +
          'PEPTIDE1,PEPTIDE2,1:R1-1:R1' + '$$$' + 'V2.0',
        mmHelm: 'PEPTIDE1{[(#3)Succ].[{A(CHOL)].F.C.T.G.H.Y.P.C.[NH2}]}$PEPTIDE1,PEPTIDE1,4:R3-10:R3$$$V2.0',
      },
      tgt: {
        templateChain: {monomerCount: [1, 10], linkageCount: 1},
        mmChain: {monomerCount: [11], linkageCount: 1}
      }
    },
    'dimerized2': {
      data: {
        templateSeq: '($3)Succ-{R-F-C(1)-T-G-H-F-P-C(1)-NH2}($3){A(CHOL)-F-C(1)-Y-H-G-D-N-C(1)-meI}',
        templateHelm: 'PEPTIDE1{[($3)Succ]}' + '|' +
          'PEPTIDE2{R.F.[C(1)].T.G.H.F.P.[C(1)].[NH2]}' + '|' +
          'PEPTIDE3{[($3)A(CHOL)].F.[C(1)].Y.H.G.D.N.[C(1)].[meI]}' + '$' +
          'PEPTIDE1,PEPTIDE2,1:R1-1:R1' + '$$$' + 'V2.0',
        mmHelm: 'PEPTIDE1{[($3)Succ].[{R].F.C.T.G.H.F.P.C.[NH2}($3){A(CHOL)].F.[C(1)].Y.H.G.D.N.[C(1)].[meI}]}$PEPTIDE1,PEPTIDE1,4:R3-10:R3$$$V2.0',
      },
      tgt: {
        templateChain: {monomerCount: [1, 10, 10], linkageCount: 1},
        mmChain: {monomerCount: [20], linkageCount: 1}
      }
    }
  };

  for (const [testName, {data, tgt}] of Object.entries(tests)) {
    test(`fromNotation-${testName}`, async () => {
      const rules = await getRules(['rules_example.json']);
      const resMmChain = Chain.fromSeparator(data.templateSeq, helmHelper);
      resMmChain.applyRules(rules);
      resMmChain.check(true);
      expectArray(resMmChain.monomers.map((mL) => mL.length), tgt.mmChain.monomerCount);
      expect(resMmChain.linkages.length, tgt.mmChain.linkageCount);
      expect(resMmChain.getHelm(), data.mmHelm);
    }, testName == 'reaction2' ? {skipReason: 'reverse reaction'} : undefined);
  }

  for (const [testName, {data, tgt}] of Object.entries(tests)) {
    test(`parseNotation-${testName}`, async () => {
      const rules = await getRules(['rules_example.json']);
      const resTemplateChain = Chain.fromSeparator(data.templateSeq, helmHelper);
      resTemplateChain.check(true);
      expectArray(resTemplateChain.monomers.map((mL) => mL.length), tgt.templateChain.monomerCount);
      expect(resTemplateChain.linkages.length, tgt.templateChain.linkageCount);
      expect(resTemplateChain.getHelm(), data.templateHelm);
      expect(resTemplateChain.getNotation(), data.templateSeq);
    });
  }

  for (const [testName, {data, tgt}] of Object.entries(tests)) {
    test(`parseHelm-${testName}`, async () => {
      const rules = await getRules(['rules_example.json']);
      const resTemplateChain = Chain.fromHelm(data.templateHelm, helmHelper);
      resTemplateChain.check(true);
      expectArray(resTemplateChain.monomers.map((mL) => mL.length), tgt.templateChain.monomerCount);
      expect(resTemplateChain.linkages.length, tgt.templateChain.linkageCount);
      expect(resTemplateChain.getHelm(), data.templateHelm);
      expect(resTemplateChain.getNotation(), data.templateSeq);
    });
  }

  for (const [testName, {data, tgt}] of Object.entries(tests)) {
    test(`applyRules-${testName}`, async () => {
      const rules = await getRules(['rules_example.json']);
      const resTemplateChain = Chain.fromSeparator(data.templateSeq, helmHelper);
      resTemplateChain.applyRules(rules);
      resTemplateChain.check(true);
      expectArray(resTemplateChain.monomers.map((mL) => mL.length), tgt.mmChain.monomerCount);
      expect(resTemplateChain.linkages.length, tgt.mmChain.linkageCount);
      expect(resTemplateChain.getHelm(), data.mmHelm);
    }, {skipReason: 'applyRules is not implemented'});
  }
});
