import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {ALPHABET, NOTATION, TAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';

const seqDna = `seq
ACGTCACGTC
CAGTGTCAGTGT
TTCAACTTCAAC`;

const seqDnaMsa = `seq
AC-GT-CTAC-GT-CT
CAC-T-GTCAC-T-GT
ACCGTACTACCGTACT`;

const seqUn = `seq
abc-dfgg-abc1-cfr3-rty-wert-abc-dfgg-abc1-cfr3-rty-wert
rut12-her2-rty-wert-abc-abc1-dfgg-rut12-her2-rty-wert-abc-abc1-dfgg
rut12-rty-her2-abc-cfr3-wert-rut12-rut12-rty-her2-abc-cfr3-wert-rut12`;

const seqHelm = `seq
PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.D-Orn.D-aThr.Phe_4Me}$$$$
PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.Tyr_PO3H2.D-Chg.dV.Phe_ab-dehydro.N.D-Orn.D-aThr.Phe_4Me}$$$$
PEPTIDE1{Lys_Boc.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.Tyr_PO3H2.D-Chg.dV.Thr_PO3H2.N.D-Orn.D-aThr.Phe_4Me}$$$$
PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.Tyr_PO3H2.D-Chg.dV.Thr_PO3H2.N.D-Orn.D-aThr.Phe_4Me}$$$$`;

category('UnitsHandler', () => {
  test('Seq-Fasta', async () => {
    const [_df, uh] = await loadCsvWithDetection(seqDna);
    expect(uh.notation, NOTATION.FASTA);
    expect(uh.isMsa(), false);
  });

  test('Seq-Fasta-MSA', async () => {
    const [_df, uh] = await loadCsvWithDetection(seqDnaMsa);
    expect(uh.notation, NOTATION.FASTA);
    expect(uh.isMsa(), true);
  });

  test('Seq-Fasta-units', async () => {
    const [_df, uh] = await loadCsvWithDetection(seqDna);
    expect(uh.notation, NOTATION.FASTA);
    expect(uh.isMsa(), false);
  });

  test('Seq-Fasta-MSA-units', async () => {
    const [_df, uh] = await loadCsvWithDetection(seqDnaMsa);
    expect(uh.notation, NOTATION.FASTA);
    expect(uh.isMsa(), true);
  });

  test('Seq-Helm', async () => {
    const [_df, uh] = await loadCsvWithDetection(seqHelm);
    expect(uh.notation, NOTATION.HELM);
    expect(uh.isHelm(), true);
  });

  test('Seq-UN', async () => {
    const [_df, uh] = await loadCsvWithDetection(seqUn);
    expect(uh.notation, NOTATION.SEPARATOR);
    expect(uh.separator, '-');
    expect(uh.alphabet, ALPHABET.UN);
  });

  test('Seq-UN-auto', async () => {
    const [_df, uh] = await loadCsvWithDetection(seqUn);
    expect(uh.notation, NOTATION.SEPARATOR);
    expect(uh.separator, '-');
    expect(uh.alphabet, ALPHABET.UN);
  });

  async function loadCsvWithDetection(csv: string): Promise<[df: DG.DataFrame, uh: UnitsHandler]> {
    const df = DG.DataFrame.fromCsv(csv);
    await grok.data.detectSemanticTypes(df);
    const uh = UnitsHandler.getOrCreate(df.getCol('seq'));
    return [df, uh];
  }

  // async function loadCsvWithTag(csv: string, tag: string, value: string):
  //   Promise<[df: DG.DataFrame, uh: UnitsHandler]> {
  //   const df = DG.DataFrame.fromCsv(csv);
  //   const col = df.getCol('seq');
  //   col.setTag(tag, value);
  //   col.semType = DG.SEMTYPE.MACROMOLECULE;
  //   if (value === NOTATION.SEPARATOR)
  //     col.setTag(TAGS.separator, '-');
  //   const uh = UnitsHandler.getOrCreate(df.getCol('seq'));
  //   return [df, uh];
  // }
});
