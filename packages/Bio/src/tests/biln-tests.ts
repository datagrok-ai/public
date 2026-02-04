/* eslint-disable max-len */
/* eslint-disable max-lines */
/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, before} from '@datagrok-libraries/test/src/test';

import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqHelper, getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {_testNeg, _testPos, DetectorTestData, DfReaderFunc, PosCol} from './utils/detectors-utils';
import {SeqTemps} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';

export async function prepareBiln(list: string[], seqHelper: ISeqHelper): Promise<DG.Column> {
  const col: DG.Column = DG.Column.fromList(DG.TYPE.STRING, 'seq', list);
  const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: col});
  if (semType)
    col.semType = semType;
  expect(col.semType, DG.SEMTYPE.MACROMOLECULE);
  const sh = seqHelper.getSeqHandler(col);
  await sh.refinerPromise; // wait for refiner to finish
  const _newSh = seqHelper.getSeqHandler(col);
  return col;
}


export async function _testBilnDetection(list: string[], seqHelper: ISeqHelper, negativeTest = false): Promise<void> {
  const col = await prepareBiln(list, seqHelper);
  expect(col.meta.units === NOTATION.BILN, !negativeTest, `Incorrectly detected as ${col.meta.units}`);
  expect(col.temp[SeqTemps.notationProvider] != null, !negativeTest, `No notation provider for BILN`);
}

export async function _testBilnToHelm(list: string[], expectedHelm: string[], seqHelper: ISeqHelper): Promise<void> {
  const col = await prepareBiln(list, seqHelper);
  const sh = seqHelper.getSeqHandler(col);
  for (let i = 0; i < list.length; i++) {
    const helm = sh.getHelm(i);
    expect(helm === expectedHelm[i], true, `Incorrect HELM conversion for ${list[i]}: Expected ${expectedHelm[i]} \n Got ${helm}`);
  }
  // also test through converter
  const converter = sh.getConverter(NOTATION.HELM);
  for (let i = 0; i < list.length; i++) {
    const helm = converter(list[i]);
    expect(helm === expectedHelm[i], true, `Incorrect HELM conversion for ${list[i]}: Expected ${expectedHelm[i]} \n Got ${helm}`);
  }
}

export async function _testHelmToBiln(helmList: string[], expectedBiln: string[], seqHelper: ISeqHelper): Promise<void> {
  const col: DG.Column = DG.Column.fromList(DG.TYPE.STRING, 'helm', helmList);
  const df = DG.DataFrame.fromColumns([col]);
  await df.meta.detectSemanticTypes();
  await grok.data.detectSemanticTypes(df);
  expect(col.semType === DG.SEMTYPE.MACROMOLECULE, true, `Incorrectly detected as ${col.semType}`);
  expect(col.meta.units === NOTATION.HELM, true, `Incorrectly detected as ${col.meta.units}`);
  const sh = seqHelper.getSeqHandler(col);
  const converter = sh.getConverter(NOTATION.BILN);
  for (let i = 0; i < helmList.length; i++) {
    const biln = converter(helmList[i]);
    expect(biln === expectedBiln[i], true, `Incorrect BILN conversion for ${helmList[i]}: Expected ${expectedBiln[i]} \n Got ${biln}`);
  }
}

export const detectorTestsDataForBiln: {name: string, seqs: string[], negative: boolean}[] = [
  {
    name: 'Valid Biln',
    seqs: [
      'A-C(1,3)-G-G-H-A-V-E-A-K-Y-L-V-C(3,3)-S.G-I-V-E-A-C(2,3)-C(1,3)-T-S-I-C(2,3)-S-L-Y-Q-L-E-N-Y-C(3,3)-Y',
      'A-C(1,3)-G-G-H-A-V-E-A-K-Y-L-V-C(3,3)-S.G-I-V-E-A-C(2,3)-C(1,3)-T-S-I-C(2,3)-S-L-Y-Q-L-E-N-Y-C(3,3)-Y',
      'C-C(1,3)-S-W-P-A-R-C(2,3)-L-H-Q-D-L-C(3,3)-NH2.C(1,1)(2,2)(3,3)',
      'C-C(1,3)-S-W-P-A-R-C(2,3)-L-H-Q-D-L-C(3,3)-NH2.[C](1,1)(2,2)(3,3)',
      'D-T-H-F-P-I-C(1,3)-I-F-C(2,3)-C(3,3)-G-C(2,3)-C(4,3)-H-R-S-K-C(3,3)-G-M-C(4,3)-C(1,3)-K-T'
    ], negative: false
  }, {
    name: 'Simple Separator Neg',
    seqs: [
      'meI/hHis/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr//Phe_4Me',
      'meI/hHis/Aca/Cys_SEt/T/dK/Thr_PO3H2/Aca/Tyr_PO3H2/D-Chg/dV/Phe_ab-dehydro/N/D-Orn/D-aThr//Phe_4Me',
      'Lys_Boc/hHis/Aca/Cys_SEt/T/dK/Thr_PO3H2/Aca/Tyr_PO3H2/D-Chg/dV/Thr_PO3H2/N/D-Orn/D-aThr//Phe_4Me',
      'meI/hHis/Aca/Cys_SEt/T/dK/Thr_PO3H2/Aca/Tyr_PO3H2/D-Chg/dV/Thr_PO3H2/N/D-Orn/D-aThr//Phe_4Me',
      'meI/hHis/Aca/N/T/dK/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/Chg/N/D-Orn/D-aThr//Phe_4Me'
    ], negative: true
  }, {
    name: 'Valid Biln without cyclization',
    seqs: [
      'meI-hHis-Aca-N-T-dE-Thr_PO3H2-Aca-[D-Tyr_Et]-[Tyr_ab-dehydroMe]-dV-E-N-[D-Orn]-[D-aThr]--Phe_4Me',
      'meI-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2-[D-Chg]-dV-[Phe_ab-dehydro]-N-[D-Orn]-[D-aThr]--Phe_4Me',
      'Lys_Boc-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2-[D-Chg]-dV-Thr_PO3H2-N-[D-Orn]-[D-aThr]--Phe_4Me',
      'meI-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2-[D-Chg]-dV-Thr_PO3H2-N-[D-Orn]-[D-aThr]--Phe_4Me',
      'meI-hHis-Aca-N-T-dK-Thr_PO3H2-Aca-[D-Tyr_Et]-[Tyr_ab-dehydroMe]-dV-Chg-N-[D-Orn]-[D-aThr]--Phe_4Me',
      'meI-hHis-Aca-N-T-dK-Thr_PO3H2-Aca-[D-Tyr_Et]-Tyr_Bn-dV-E-N-dV---Phe_4Me',
      'meI-hHis-Aca-N-T-dK-Thr_PO3H2-Aca-[D-Tyr_Et]-Aze-dV-E-N-dV---Phe_4Me',
      'meI-hHis-Aca-N-T-dK-Thr_PO3H2-Aca-[D-Tyr_Et]-meQ-dV-E-N-dV---Phe_4Me'
    ], negative: false
  }
];

export const bilnToHelmTestsData: {name: string, biln: string[], helm: string[]}[] = [
  {
    name: 'Linear',
    biln: ['meI-hHis-Aca-N-T-dE-Thr_PO3H2-Aca-[D-Tyr_Et]-[Tyr_ab-dehydroMe]-dV-E-N-[D-Orn]-[D-aThr]--Phe_4Me',
      'meI-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2-[D-Chg]-dV-[Phe_ab-dehydro]-N-[D-Orn]-[D-aThr]--Phe_4Me',
      'Lys_Boc-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2-[D-Chg]-dV-Thr_PO3H2-N-[D-Orn]-[D-aThr]--Phe_4Me',
      'meI-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2-[D-Chg]-dV-Thr_PO3H2-N-[D-Orn]-[D-aThr]--Phe_4Me',
      'meI-hHis-Aca-N-T-dK-Thr_PO3H2-Aca-[D-Tyr_Et]-[Tyr_ab-dehydroMe]-dV-Chg-N-[D-Orn]-[D-aThr]--Phe_4Me'],
    helm: [
      'PEPTIDE1{[meI].[hHis].[Aca].N.T.[dE].[Thr_PO3H2].[Aca].[D-Tyr_Et].[Tyr_ab-dehydroMe].[dV].E.N.[D-Orn].[D-aThr].*.[Phe_4Me]}$$$$',
      'PEPTIDE1{[meI].[hHis].[Aca].[Cys_SEt].T.[dK].[Thr_PO3H2].[Aca].[Tyr_PO3H2].[D-Chg].[dV].[Phe_ab-dehydro].N.[D-Orn].[D-aThr].*.[Phe_4Me]}$$$$',
      'PEPTIDE1{[Lys_Boc].[hHis].[Aca].[Cys_SEt].T.[dK].[Thr_PO3H2].[Aca].[Tyr_PO3H2].[D-Chg].[dV].[Thr_PO3H2].N.[D-Orn].[D-aThr].*.[Phe_4Me]}$$$$',
      'PEPTIDE1{[meI].[hHis].[Aca].[Cys_SEt].T.[dK].[Thr_PO3H2].[Aca].[Tyr_PO3H2].[D-Chg].[dV].[Thr_PO3H2].N.[D-Orn].[D-aThr].*.[Phe_4Me]}$$$$',
      'PEPTIDE1{[meI].[hHis].[Aca].N.T.[dK].[Thr_PO3H2].[Aca].[D-Tyr_Et].[Tyr_ab-dehydroMe].[dV].[Chg].N.[D-Orn].[D-aThr].*.[Phe_4Me]}$$$$'
    ]
  }, {
    name: 'Cyclic',
    biln: [
      'C-C(1,3)-S-W-P-A-R-C(2,3)-L-H-Q-D-L-C(3,3)-NH2.[C](1,1)(2,2)(3,3)',
      'D-T-H-F-P-I-C(1,3)-I-F-C(2,3)-C(3,3)-G-C(2,3)-C(4,3)-H-R-S-K-C(3,3)-G-M-C(4,3)-C(1,3)-K-T',
      'L-C(1,3)-G-S-H-L-V-E-A-L-Y-L-V-C(2,3)-G.G-I-V-E-Q-C(3,3)-C(1,3)-T-S-I-C(3,3)-S-L-Y-Q-L-E-N-Y-C(2,3)-N',
      'H-Aib-E-G-T-F-T-S-D(2,3)-V-S-S-Y-L-E-G-Q-A-A-K(1,3)-E-F-I-A-W-L-V-R-G-R-G.C(2,3)-gGlu-G-G(1,2)',
      'F(4,2).dI(1,1)(3,3)-Trp_Ome-Asp_OMe-Cys_Bn-meG-Phe_3Cl-dD-T-dI(4,3)-T-dK-aG-3Pal-xiIle-meD-Ala_tBu(1,2).L(2,1)-Pro_4Me3OH-S-NMe2Abz-Q-3Pal-xiIle-D-Hyp-Ala_tBu-dI(3,3)-Trp_Ome-Asp_OMe-N-meG-Phe_34diCl-Phe_34diCl(2,2)'
    ],
    helm: [
      'PEPTIDE1{C.C.S.W.P.A.R.C.L.H.Q.D.L.C.[NH2]}|PEPTIDE2{C}$PEPTIDE1,PEPTIDE2,14:R3-1:R3|PEPTIDE1,PEPTIDE2,8:R3-1:R2|PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0',
      'PEPTIDE1{D.T.H.F.P.I.C.I.F.C.C.G.C.C.H.R.S.K.C.G.M.C.C.K.T}$PEPTIDE1,PEPTIDE1,10:R3-13:R3|PEPTIDE1,PEPTIDE1,11:R3-19:R3|PEPTIDE1,PEPTIDE1,14:R3-22:R3|PEPTIDE1,PEPTIDE1,7:R3-23:R3$$$V2.0',
      'PEPTIDE1{L.C.G.S.H.L.V.E.A.L.Y.L.V.C.G}|PEPTIDE2{G.I.V.E.Q.C.C.T.S.I.C.S.L.Y.Q.L.E.N.Y.C.N}$PEPTIDE1,PEPTIDE2,2:R3-7:R3|PEPTIDE2,PEPTIDE2,6:R3-11:R3|PEPTIDE1,PEPTIDE2,14:R3-20:R3$$$V2.0',
      'PEPTIDE1{H.[Aib].E.G.T.F.T.S.D.V.S.S.Y.L.E.G.Q.A.A.K.E.F.I.A.W.L.V.R.G.R.G}|PEPTIDE2{C.[gGlu].G.G}$PEPTIDE1,PEPTIDE2,9:R3-1:R3|PEPTIDE1,PEPTIDE2,20:R3-4:R2$$$V2.0',
      'PEPTIDE1{F}|PEPTIDE2{[dI].[Trp_Ome].[Asp_OMe].[Cys_Bn].[meG].[Phe_3Cl].[dD].T.[dI].T.[dK].[aG].[3Pal].[xiIle].[meD].[Ala_tBu]}|PEPTIDE3{L.[Pro_4Me3OH].S.[NMe2Abz].Q.[3Pal].[xiIle].D.[Hyp].[Ala_tBu].[dI].[Trp_Ome].[Asp_OMe].N.[meG].[Phe_34diCl].[Phe_34diCl]}$PEPTIDE1,PEPTIDE2,1:R2-9:R3|PEPTIDE2,PEPTIDE2,1:R1-16:R2|PEPTIDE2,PEPTIDE3,1:R3-11:R3|PEPTIDE3,PEPTIDE3,1:R1-17:R2$$$V2.0'
    ]
  }
];

export const helmToBilnTestsData: {name: string, helm: string[], biln: string[]}[] = [
  {
    name: 'Linear',
    biln: [
      'meI-hHis-Aca-N-T-dE-Thr_PO3H2-Aca-[D-Tyr_Et]-[Tyr_ab-dehydroMe]-dV-E-N-[D-Orn]-[D-aThr]--Phe_4Me',
      'meI-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2-[D-Chg]-dV-[Phe_ab-dehydro]-N-[D-Orn]-[D-aThr]--Phe_4Me',
      'Lys_Boc-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2-[D-Chg]-dV-Thr_PO3H2-N-[D-Orn]-[D-aThr]--Phe_4Me',
      'meI-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2-[D-Chg]-dV-Thr_PO3H2-N-[D-Orn]-[D-aThr]--Phe_4Me',
      'meI-hHis-Aca-N-T-dK-Thr_PO3H2-Aca-[D-Tyr_Et]-[Tyr_ab-dehydroMe]-dV-Chg-N-[D-Orn]-[D-aThr]--Phe_4Me'
    ],
    helm: [
      'PEPTIDE1{[meI].[hHis].[Aca].N.T.[dE].[Thr_PO3H2].[Aca].[D-Tyr_Et].[Tyr_ab-dehydroMe].[dV].E.N.[D-Orn].[D-aThr].*.[Phe_4Me]}$$$$',
      'PEPTIDE1{[meI].[hHis].[Aca].[Cys_SEt].T.[dK].[Thr_PO3H2].[Aca].[Tyr_PO3H2].[D-Chg].[dV].[Phe_ab-dehydro].N.[D-Orn].[D-aThr].*.[Phe_4Me]}$$$$',
      'PEPTIDE1{[Lys_Boc].[hHis].[Aca].[Cys_SEt].T.[dK].[Thr_PO3H2].[Aca].[Tyr_PO3H2].[D-Chg].[dV].[Thr_PO3H2].N.[D-Orn].[D-aThr].*.[Phe_4Me]}$$$$',
      'PEPTIDE1{[meI].[hHis].[Aca].[Cys_SEt].T.[dK].[Thr_PO3H2].[Aca].[Tyr_PO3H2].[D-Chg].[dV].[Thr_PO3H2].N.[D-Orn].[D-aThr].*.[Phe_4Me]}$$$$',
      'PEPTIDE1{[meI].[hHis].[Aca].N.T.[dK].[Thr_PO3H2].[Aca].[D-Tyr_Et].[Tyr_ab-dehydroMe].[dV].[Chg].N.[D-Orn].[D-aThr].*.[Phe_4Me]}$$$$'
    ]
  }, {
    name: 'Cyclic',
    biln: [
      'C-C(3,3)-S-W-P-A-R-C(2,3)-L-H-Q-D-L-C(1,3)-NH2.C(1,3)(2,2)(3,1)',
      'D-T-H-F-P-I-C(4,3)-I-F-C(1,3)-C(2,3)-G-C(1,3)-C(3,3)-H-R-S-K-C(2,3)-G-M-C(3,3)-C(4,3)-K-T',
      'L-C(1,3)-G-S-H-L-V-E-A-L-Y-L-V-C(3,3)-G.G-I-V-E-Q-C(2,3)-C(1,3)-T-S-I-C(2,3)-S-L-Y-Q-L-E-N-Y-C(3,3)-N',
      'H-Aib-E-G-T-F-T-S-D(1,3)-V-S-S-Y-L-E-G-Q-A-A-K(2,3)-E-F-I-A-W-L-V-R-G-R-G.C(1,3)-gGlu-G-G(2,2)',
      'F(1,2).dI(2,1)(3,3)-Trp_Ome-Asp_OMe-Cys_Bn-meG-Phe_3Cl-dD-T-dI(1,3)-T-dK-aG-3Pal-xiIle-meD-Ala_tBu(2,2).L(4,1)-Pro_4Me3OH-S-NMe2Abz-Q-3Pal-xiIle-D-Hyp-Ala_tBu-dI(3,3)-Trp_Ome-Asp_OMe-N-meG-Phe_34diCl-Phe_34diCl(4,2)'
    ],
    helm: [
      'PEPTIDE1{C.C.S.W.P.A.R.C.L.H.Q.D.L.C.[NH2]}|PEPTIDE2{C}$PEPTIDE1,PEPTIDE2,14:R3-1:R3|PEPTIDE1,PEPTIDE2,8:R3-1:R2|PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0',
      'PEPTIDE1{D.T.H.F.P.I.C.I.F.C.C.G.C.C.H.R.S.K.C.G.M.C.C.K.T}$PEPTIDE1,PEPTIDE1,10:R3-13:R3|PEPTIDE1,PEPTIDE1,11:R3-19:R3|PEPTIDE1,PEPTIDE1,14:R3-22:R3|PEPTIDE1,PEPTIDE1,7:R3-23:R3$$$V2.0',
      'PEPTIDE1{L.C.G.S.H.L.V.E.A.L.Y.L.V.C.G}|PEPTIDE2{G.I.V.E.Q.C.C.T.S.I.C.S.L.Y.Q.L.E.N.Y.C.N}$PEPTIDE1,PEPTIDE2,2:R3-7:R3|PEPTIDE2,PEPTIDE2,6:R3-11:R3|PEPTIDE1,PEPTIDE2,14:R3-20:R3$$$V2.0',
      'PEPTIDE1{H.[Aib].E.G.T.F.T.S.D.V.S.S.Y.L.E.G.Q.A.A.K.E.F.I.A.W.L.V.R.G.R.G}|PEPTIDE2{C.[gGlu].G.G}$PEPTIDE1,PEPTIDE2,9:R3-1:R3|PEPTIDE1,PEPTIDE2,20:R3-4:R2$$$V2.0',
      'PEPTIDE1{F}|PEPTIDE2{[dI].[Trp_Ome].[Asp_OMe].[Cys_Bn].[meG].[Phe_3Cl].[dD].T.[dI].T.[dK].[aG].[3Pal].[xiIle].[meD].[Ala_tBu]}|PEPTIDE3{L.[Pro_4Me3OH].S.[NMe2Abz].Q.[3Pal].[xiIle].D.[Hyp].[Ala_tBu].[dI].[Trp_Ome].[Asp_OMe].N.[meG].[Phe_34diCl].[Phe_34diCl]}$PEPTIDE1,PEPTIDE2,1:R2-9:R3|PEPTIDE2,PEPTIDE2,1:R1-16:R2|PEPTIDE2,PEPTIDE3,1:R3-11:R3|PEPTIDE3,PEPTIDE3,1:R1-17:R2$$$V2.0'
    ]
  }
];
