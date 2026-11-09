import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {ISeqHelper, getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {GapOriginals} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {ISeqGraphInfoHelm} from '@datagrok-libraries/bio/src/utils/macromolecule/types';

enum Tests {
  fasta = 'fasta',
  fastaMsa = 'fastaMsa',
  separator = 'separator',
  separatorMsa = 'separatorMsa',
  helm = 'helm',
}

category('SeqHandler: splitted', () => {
  let seqHelper: ISeqHelper;

  before(async () => {
    seqHelper = await getSeqHelper();
  });

  const fG = GapOriginals[NOTATION.FASTA];
  const hG = GapOriginals[NOTATION.HELM];
  const sG = GapOriginals[NOTATION.SEPARATOR];
  const data: {
    [testName: string]: {
      src: { csv: string },
      tgt: { notation: NOTATION, separator?: string, splitted: (string[] | string)[] }
    }
  } = {
    [Tests.fasta]: {
      src: {
        csv: `seq
ACGTCACGTC
CAGTGTCAGTGT
TTCAACTTCAAC`
      },
      tgt: {
        notation: NOTATION.FASTA,
        splitted: [
          ['A', 'C', 'G', 'T', 'C', 'A', 'C', 'G', 'T', 'C'],
          ['C', 'A', 'G', 'T', 'G', 'T', 'C', 'A', 'G', 'T', 'G', 'T'],
          ['T', 'T', 'C', 'A', 'A', 'C', 'T', 'T', 'C', 'A', 'A', 'C'],
        ]
      }
    },
    [Tests.fastaMsa]: {
      src: {
        csv: `seq
AC-GT-CTAC-GT-CT
CAC-T-GTCAC-T-GT
ACCGTACTACCGTACT`,
      },
      tgt: {
        notation: NOTATION.FASTA,
        splitted: [
          //@formatter:off
          ['A', 'C', '-', 'G', 'T', '-', 'C', 'T', 'A', 'C', '-', 'G', 'T', '-', 'C', 'T'],
          ['C', 'A', 'C', '-', 'T', '-', 'G', 'T', 'C', 'A', 'C', '-', 'T', '-', 'G', 'T'],
          ['A', 'C', 'C', 'G', 'T', 'A', 'C', 'T', 'A', 'C', 'C', 'G', 'T', 'A', 'C', 'T'],
          //@formatter:on
        ]
      }
    },
    [Tests.separator]: {
      src: {
        csv: `seq
abc-dfgg-abc1-cfr3-rty-wert-abc-dfgg-abc1-cfr3-rty-wert
rut12-her2-rty-wert-abc-abc1-dfgg-rut12-her2-rty-wert-abc
rut12-rty-her2-abc-cfr3-wert-rut12-rut12-rty-her2-abc-cfr3`,
      },
      tgt: {
        notation: NOTATION.SEPARATOR,
        separator: '-',
        splitted: [
          ['abc', 'dfgg', 'abc1', 'cfr3', 'rty', 'wert', 'abc', 'dfgg', 'abc1', 'cfr3', 'rty', 'wert'],
          ['rut12', 'her2', 'rty', 'wert', 'abc', 'abc1', 'dfgg', 'rut12', 'her2', 'rty', 'wert', 'abc'],
          ['rut12', 'rty', 'her2', 'abc', 'cfr3', 'wert', 'rut12', 'rut12', 'rty', 'her2', 'abc', 'cfr3']
        ]
      }
    },

    [Tests.separatorMsa]: {
      src: {
        csv: `seq
rut0-dfgg-abc1-cfr3-rty-wert-abc-dfgg-abc1-cfr3-rty-wert
rut1-her2-rty--abc1-dfgg-rut12-her2-rty--abc1-dfgg
rut2-rty-her2---wert-rut12-rty-her2---wert
\"rut3-rty-her2-\"\"-\"\"-\"\"-\"\"-wert-rut12-rty-her2-\"\"-\"\"-\"\"-\"\"-wert\"
\"\"\"-\"\"-rut4-her2-wert-rut12-rty-her2-wert\"
\"rut5-rty-her2-wert-rut12-rty-her2-wert-\"\"-\"\"\"`
      },
      tgt: {
        notation: NOTATION.SEPARATOR,
        separator: '-',
        splitted: [
          ['rut0', 'dfgg', 'abc1', 'cfr3', 'rty', 'wert', 'abc', 'dfgg', 'abc1', 'cfr3', 'rty', 'wert'],
          ['rut1', 'her2', 'rty', sG, 'abc1', 'dfgg', 'rut12', 'her2', 'rty', sG, 'abc1', 'dfgg'],
          ['rut2', 'rty', 'her2', sG, sG, 'wert', 'rut12', 'rty', 'her2', sG, sG, 'wert'],
          ['rut3', 'rty', 'her2', sG, sG, 'wert', 'rut12', 'rty', 'her2', sG, sG, 'wert'],
          [sG, 'rut4', 'her2', 'wert', 'rut12', 'rty', 'her2', 'wert'],
          ['rut5', 'rty', 'her2', 'wert', 'rut12', 'rty', 'her2', 'wert', sG],
        ]
      }
    },

    [Tests.helm]: {
      src: {
        csv: `seq
PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.Thr_PO3H2.Aca.D-Tyr_Et}$$$$
PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.dK.Thr_PO3H2.Aca}$$$$
PEPTIDE1{Lys_Boc.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.dK.Thr_PO3H2.Aca}$$$$
PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.T.dK.Thr_PO3H2}$$$$,
PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.T.dK}|PEPTIDE2{Thr_PO3H2}$$$$`,
      },
      tgt: {
        notation: NOTATION.HELM,
        splitted: [
          ['meI', 'hHis', 'Aca', 'N', 'T', 'dE', 'Thr_PO3H2', 'Aca', 'D-Tyr_Et', 'Thr_PO3H2', 'Aca', 'D-Tyr_Et'],
          ['meI', 'hHis', 'Aca', 'Cys_SEt', 'T', 'dK', 'Thr_PO3H2', 'Aca', 'dK', 'Thr_PO3H2', 'Aca'],
          ['Lys_Boc', 'hHis', 'Aca', 'Cys_SEt', 'T', 'dK', 'Thr_PO3H2', 'Aca', 'dK', 'Thr_PO3H2', 'Aca'],
          ['meI', 'hHis', 'Aca', 'Cys_SEt', 'T', 'dK', 'Thr_PO3H2', 'T', 'dK', 'Thr_PO3H2'],
          ['meI', 'hHis', 'Aca', 'Cys_SEt', 'T', 'dK', 'Thr_PO3H2', 'T', 'dK', 'Thr_PO3H2'],
        ]
      }
    }
  };

  for (const [testName, testData] of Object.entries(data)) {
    test(`${testName}`, async () => {
      const df: DG.DataFrame = DG.DataFrame.fromCsv(testData.src.csv);
      const col: DG.Column = df.getCol('seq');

      const semType = await grok.functions.call('Bio:detectMacromolecule', {col: col});
      if (semType) col.semType = semType;
      expect(col.semType, DG.SEMTYPE.MACROMOLECULE);

      const sh = seqHelper.getSeqHandler(col);
      expect(sh.notation, testData.tgt.notation);
      expect(sh.separator === testData.tgt.separator, true);

      const resSplitted: string[][] = wu.count(0).take(sh.length)
        .map((rowIdx) => {
          const seqSS = sh.getSplitted(rowIdx);
          return wu.count(0).take(seqSS.length).map((posIdx) => seqSS.getOriginal(posIdx)).toArray();
        }).toArray();
      expectArray(resSplitted, testData.tgt.splitted);
    });
  }

  // splitterAsHelm on the hwe parser (Bio on Crux, PLAN §5.2): mixed polymer types (B1), annotated chains (B2),
  // positions >= 100 (B3), inline SMILES with '.'
  test('helm-hwe-fixed', async () => {
    const big = `PEPTIDE1{${new Array(120).fill('A').join('.')}}$PEPTIDE1,PEPTIDE1,120:R2-1:R1$$$V2.0`;
    const col = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'seq', [
      'PEPTIDE1{A.C.D}|CHEM1{[PEG2]}$PEPTIDE1,CHEM1,3:R2-1:R1$$$V2.0',
      'PEPTIDE1{A.C.D}"HC"|PEPTIDE2{G.C}$$$$V2.0',
      big,
      'CHEM1{[C(=O)O.N[*:1]]}|PEPTIDE1{A.G.C}$PEPTIDE1,CHEM1,3:R3-1:R1$$$V2.0',
    ]) as DG.Column<string>;
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = NOTATION.HELM;
    const sh = seqHelper.getSeqHandler(col);
    const originals = (rowIdx: number): string[] => {
      const ss = sh.getSplitted(rowIdx);
      return wu.count(0).take(ss.length).map((posIdx) => ss.getOriginal(posIdx)).toArray();
    };
    const conns = (rowIdx: number): string[] => sh.getSplitted(rowIdx).graphInfo!.connections.map((c) =>
      [c.seqIndex1, c.monomerIndex1, c.rGroup1, c.seqIndex2, c.monomerIndex2, c.rGroup2].join(','));

    expectArray(originals(0), ['A', 'C', 'D', 'PEG2']);
    expectArray(conns(0), ['0,2,2,1,0,1']);
    expectArray(originals(1), ['A', 'C', 'D', 'G', 'C']);
    expectArray(sh.getSplitted(1).graphInfo!.disjointSeqStarts, [0, 3]);
    expect(sh.getSplitted(2).length, 120);
    expectArray(conns(2), ['0,119,2,0,0,1']);
    expectArray(originals(3), ['C(=O)O.N[*:1]', 'A', 'G', 'C']);
    expectArray(conns(3), ['1,2,3,0,0,1']);
  });

  // hwe 1.1.0: whitespace (B10), a BLOB type with a space, repeats listed once as written, an ambiguous connection
  // (left out of connections, an error for assembly)
  test('helm-hwe-1.1.0', async () => {
    const col = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'seq', [
      'PEPTIDE1{A . C}$$$$V2.0 ',
      'BLOB1{Gold Particle}$$$$V2.0',
      'RNA1{R(C)P.(R(U)P)\'3\'.R(A)P}$$$$V2.0',
      'PEPTIDE1{A.C.K.D}|CHEM1{[sDBL]}$PEPTIDE1,CHEM1,(1,3):R3-1:R1$$$V2.0',
    ]) as DG.Column<string>;
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = NOTATION.HELM;
    const sh = seqHelper.getSeqHandler(col);
    const originals = (rowIdx: number): string[] => {
      const ss = sh.getSplitted(rowIdx);
      return wu.count(0).take(ss.length).map((posIdx) => ss.getOriginal(posIdx)).toArray();
    };
    const giOf = (rowIdx: number) => sh.getSplitted(rowIdx).graphInfo as ISeqGraphInfoHelm;

    expectArray(originals(0), ['A', 'C']);
    expect(giOf(0).errors!.length, 0);
    expectArray(originals(1), ['Gold Particle']);
    expectArray(originals(2), ['R', 'C', 'P', 'R', 'U', 'P', 'R', 'A', 'P']);
    expectArray(giOf(2).repeats!.map((r) => [r.start, r.end, r.min, r.max].join(',')), ['3,6,3,3']);
    expect(giOf(3).connections.length, 0);
    expectArray(giOf(3).errors!, ['ambiguous connection \'PEPTIDE1,CHEM1,(1,3):R3-1:R1\' is not supported']);
  });
});
