import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import {GapSymbols, UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/types';

category('UnitsHandler', () => {
  const fG = GapSymbols[NOTATION.FASTA];
  const hG = GapSymbols[NOTATION.HELM];
  const sG = GapSymbols[NOTATION.SEPARATOR];
  const data: { [testName: string]: { src: { csv: string }, tgt: { splitted: (string[] | string)[] } } } = {
    fasta: {
      src: {
        csv: `seq
ACGTCACGTC
CAGTGTCAGTGT
TTCAACTTCAAC`
      },
      tgt: {
        splitted: [
          'ACGTCACGTC',
          'CAGTGTCAGTGT',
          'TTCAACTTCAAC',
        ]
      }
    },
    fastaMsa: {
      src: {
        csv: `seq
AC-GT-CTAC-GT-CT
CAC-T-GTCAC-T-GT
ACCGTACTACCGTACT`,
      },
      tgt: {
        splitted: [
          //@formatter:off
          'AC-GT-CTAC-GT-CT',
          'CAC-T-GTCAC-T-GT',
          'ACCGTACTACCGTACT',
          //@formatter:on
        ]
      }
    },
    separator: {
      src: {
        csv: `seq
abc-dfgg-abc1-cfr3-rty-wert-abc-dfgg-abc1-cfr3-rty-wert
rut12-her2-rty-wert-abc-abc1-dfgg-rut12-her2-rty-wert-abc
rut12-rty-her2-abc-cfr3-wert-rut12-rut12-rty-her2-abc-cfr3`,
      },
      tgt: {
        splitted: [
          ['abc', 'dfgg', 'abc1', 'cfr3', 'rty', 'wert', 'abc', 'dfgg', 'abc1', 'cfr3', 'rty', 'wert'],
          ['rut12', 'her2', 'rty', 'wert', 'abc', 'abc1', 'dfgg', 'rut12', 'her2', 'rty', 'wert', 'abc'],
          ['rut12', 'rty', 'her2', 'abc', 'cfr3', 'wert', 'rut12', 'rut12', 'rty', 'her2', 'abc', 'cfr3']
        ]
      }
    },

    separatorMsa: {
      src: {
        csv: `seq
abc-dfgg-abc1-cfr3-rty-wert-abc-dfgg-abc1-cfr3-rty-wert
rut12-her2-rty--abc1-dfgg-rut12-her2-rty--abc1-dfgg
rut12-rty-her2---wert-rut12-rty-her2---wert`
      },
      tgt: {
        splitted: [
          ['abc', 'dfgg', 'abc1', 'cfr3', 'rty', 'wert', 'abc', 'dfgg', 'abc1', 'cfr3', 'rty', 'wert'],
          ['rut12', 'her2', 'rty', sG, 'abc1', 'dfgg', 'rut12', 'her2', 'rty', sG, 'abc1', 'dfgg'],
          ['rut12', 'rty', 'her2', sG, sG, 'wert', 'rut12', 'rty', 'her2', sG, sG, 'wert'],
        ]
      }
    },
    helm: {
      src: {
        csv: `seq
PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.Thr_PO3H2.Aca.D-Tyr_Et}$$$$
PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.dK.Thr_PO3H2.Aca}$$$$
PEPTIDE1{Lys_Boc.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.dK.Thr_PO3H2.Aca}$$$$
PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.T.dK.Thr_PO3H2}$$$$`
      },
      tgt: {
        splitted: [
          ['meI', 'hHis', 'Aca', 'N', 'T', 'dE', 'Thr_PO3H2', 'Aca', 'D-Tyr_Et', 'Thr_PO3H2', 'Aca', 'D-Tyr_Et'],
          ['meI', 'hHis', 'Aca', 'Cys_SEt', 'T', 'dK', 'Thr_PO3H2', 'Aca', 'dK', 'Thr_PO3H2', 'Aca'],
          ['Lys_Boc', 'hHis', 'Aca', 'Cys_SEt', 'T', 'dK', 'Thr_PO3H2', 'Aca', 'dK', 'Thr_PO3H2', 'Aca'],
          ['meI', 'hHis', 'Aca', 'Cys_SEt', 'T', 'dK', 'Thr_PO3H2', 'T', 'dK', 'Thr_PO3H2'],
        ]
      }
    }
  };

  for (const [testName, testData] of Object.entries(data)) {
    test(`splitted-${testName}`, async () => {
      const df: DG.DataFrame = DG.DataFrame.fromCsv(testData.src.csv);
      const col: DG.Column = df.getCol('seq');

      const semType = await grok.functions.call('Bio:detectMacromolecule', {col: col});
      if (semType) col.semType = semType;
      expect(col.semType, DG.SEMTYPE.MACROMOLECULE);

      const uh = UnitsHandler.getOrCreate(col);
      const resSplitted: ISeqSplitted[] = uh.splitted;
      expectArray(resSplitted, testData.tgt.splitted);
    });
  }
});
