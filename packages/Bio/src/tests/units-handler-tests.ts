import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {ALPHABET, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

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

category('SeqHandler', () => {
  test('Seq-Fasta', async () => {
    const [_df, sh] = await loadCsvWithDetection(seqDna);
    expect(sh.notation, NOTATION.FASTA);
    expect(sh.isMsa(), false);
  });

  test('Seq-Fasta-MSA', async () => {
    const [_df, sh] = await loadCsvWithDetection(seqDnaMsa);
    expect(sh.notation, NOTATION.FASTA);
    expect(sh.isMsa(), true);
  });

  test('Seq-Fasta-units', async () => {
    const [_df, sh] = await loadCsvWithDetection(seqDna);
    expect(sh.notation, NOTATION.FASTA);
    expect(sh.isMsa(), false);
  });

  test('Seq-Fasta-MSA-units', async () => {
    const [_df, sh] = await loadCsvWithDetection(seqDnaMsa);
    expect(sh.notation, NOTATION.FASTA);
    expect(sh.isMsa(), true);
  });

  test('Seq-Helm', async () => {
    const [_df, sh] = await loadCsvWithDetection(seqHelm);
    expect(sh.notation, NOTATION.HELM);
    expect(sh.isHelm(), true);
  });

  test('Seq-UN', async () => {
    const [_df, sh] = await loadCsvWithDetection(seqUn);
    expect(sh.notation, NOTATION.SEPARATOR);
    expect(sh.separator, '-');
    expect(sh.alphabet, ALPHABET.UN);
  });

  test('Seq-UN-auto', async () => {
    const [_df, sh] = await loadCsvWithDetection(seqUn);
    expect(sh.notation, NOTATION.SEPARATOR);
    expect(sh.separator, '-');
    expect(sh.alphabet, ALPHABET.UN);
  });

  async function loadCsvWithDetection(csv: string): Promise<[df: DG.DataFrame, sh: SeqHandler]> {
    const df = DG.DataFrame.fromCsv(csv);
    await grok.data.detectSemanticTypes(df);
    const sh = SeqHandler.forColumn(df.getCol('seq'));
    return [df, sh];
  }

  // async function loadCsvWithTag(csv: string, tag: string, value: string):
  //   Promise<[df: DG.DataFrame, uh: SeqHandler]> {
  //   const df = DG.DataFrame.fromCsv(csv);
  //   const col = df.getCol('seq');
  //   col.setTag(tag, value);
  //   col.semType = DG.SEMTYPE.MACROMOLECULE;
  //   if (value === NOTATION.SEPARATOR)
  //     col.setTag(TAGS.separator, '-');
  //   const sh = SeqHandler.forColumn(df.getCol('seq'));
  //   return [df, sh];
  // }
});
