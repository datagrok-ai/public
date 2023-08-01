import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, expectArray} from '@datagrok-libraries/utils/src/test';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

category('UnitsHandler', () => {
  const data: { [testName: string]: { src: { csv: string }, tgt: { splitted: string[][] } } } = {
    fasta: {
      src: {
        csv: `seq
ACGTC
CAGTGT
TTCAAC`
      },
      tgt: {
        splitted: [
          ['A', 'C', 'G', 'T', 'C'],
          ['C', 'A', 'G', 'T', 'G', 'T'],
          ['T', 'T', 'C', 'A', 'A', 'C']
        ]
      }
    },
    fastaMsa: {
      src: {
        csv: `seq
AC-GT-CT
CAC-T-GT
ACCGTACT`,
      },
      tgt: {
        splitted: [
          //@formatter:off
          ['A', 'C', '' , 'G', 'T', '' , 'C', 'T'],
          ['C', 'A', 'C', '' , 'T', '' , 'G', 'T'],
          ['A', 'C', 'C', 'G', 'T', 'A', 'C', 'T'],
          //@formatter:on
        ]
      }
    },
    separator: {
      src: {
        csv: `seq
abc-dfgg-abc1-cfr3-rty-wert
rut12-her2-rty-wert-abc-abc1-dfgg
rut12-rty-her2-abc-cfr3-wert-rut12`,
      },
      tgt: {
        splitted: [
          ['abc', 'dfgg', 'abc1', 'cfr3', 'rty', 'wert'],
          ['rut12', 'her2', 'rty', 'wert', 'abc', 'abc1', 'dfgg'],
          ['rut12', 'rty', 'her2', 'abc', 'cfr3', 'wert', 'rut12']
        ]
      }
    },

    separatorMsa: {
      src: {
        csv: `seq
abc-dfgg-abc1-cfr3-rty-wert
rut12-her2-rty--abc1-dfgg
rut12-rty-her2---wert`
      },
      tgt: {
        splitted: [
          ['abc', 'dfgg', 'abc1', 'cfr3', 'rty', 'wert'],
          ['rut12', 'her2', 'rty', '', 'abc1', 'dfgg'],
          ['rut12', 'rty', 'her2', '', '', 'wert'],
        ]
      }
    },
    helm: {
      src: {
        csv: `seq
PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et}$$$$
PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca}$$$$
PEPTIDE1{Lys_Boc.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca}$$$$
PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2}$$$$`
      },
      tgt: {
        splitted: [
          ['meI', 'hHis', 'Aca', 'N', 'T', 'dE', 'Thr_PO3H2', 'Aca', 'D-Tyr_Et'],
          ['meI', 'hHis', 'Aca', 'Cys_SEt', 'T', 'dK', 'Thr_PO3H2', 'Aca'],
          ['Lys_Boc', 'hHis', 'Aca', 'Cys_SEt', 'T', 'dK', 'Thr_PO3H2', 'Aca'],
          ['meI', 'hHis', 'Aca', 'Cys_SEt', 'T', 'dK', 'Thr_PO3H2'],
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
      const splitted: string[][] = uh.splitted;
      expectArray(splitted, testData.tgt.splitted);
    });
  }
});
