import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {category, test, expectFloat, before, after, expect} from '@datagrok-libraries/test/src/test';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {IMonomerLibHelper, getMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';

import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

category('Scoring', () => {
  /* eslint-disable max-len */
  const sequence = 'sequence';
  const expectedSimilarity = 'expected_similarity';
  const expectedIdentity = 'expected_identity';
  /* eslint-disable max-len */
  const table = DG.DataFrame.fromCsv(`${sequence},${expectedSimilarity},${expectedIdentity}
PEPTIDE1{Aca.Orn.gGlu.Pqa.D-His_1Bn.dH.hHis.4Abz.D-Tic.D-Dap.Y.Iva.meS.F.P.F.D-1Nal}$$$$,1.0,1.0
PEPTIDE1{Iva.Gly_allyl.gGlu.Pqa.D-Dip.dH.hHis.4Abz.D-aHyp.D-Dap.Y.Iva.I.Tyr_26diMe.P.Asu.meC}$$$$,0.691,0.53
PEPTIDE1{[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal]}$$$$V2.0,0.37,0.0`
  );
  const seqCol: DG.Column<string> = table.getCol(sequence);
  seqCol.meta.units = NOTATION.HELM;
  seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
  const reference = seqCol.get(0)!;
  const shortReference = 'PEPTIDE1{Iva.Gly_allyl.gGlu.Pqa.D-Dip.dH.hHis.4Abz.D-aHyp.D-Dap.Y.Iva}$$$$';
  const longReference = 'PEPTIDE1{Iva.Gly_allyl.gGlu.Pqa.D-Dip.dH.hHis.4Abz.D-aHyp.D-Dap.Y.Iva.I.Tyr_26diMe.P.Asu.meC.I.Tyr_26diMe.P.Asu.meC}$$$$';
  /* eslint-enable max-len */

  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;

  before(async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    await monomerLibHelper.loadMonomerLibForTests(); // load default libraries
  });

  after(async () => {
    // UserDataStorage.put() replaces existing data
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true); // load user settings libraries
  });

  test('Identity', async () => {
    const scoresCol = await grok.functions.call('Bio:sequenceIdentityScoring',
      {table: table, macromolecule: seqCol, reference: reference}) as DG.Column<number>;
    for (let i = 0; i < scoresCol.length; i++) {
      const resScore = scoresCol.get(i)!;
      const tgtScore = table.get(expectedIdentity, i);
      expectFloat(resScore, tgtScore, 0.01,
        `Wrong identity score for sequence at position ${i}`);
    }
  });

  test('Identity-shortReference', async () => {
    const scoresCol = await grok.functions.call('Bio:sequenceIdentityScoring',
      {table: table, macromolecule: seqCol, reference: shortReference}) as DG.Column<number>;
    expect(wu.count(0).take(scoresCol.length).map((rowI) => scoresCol.get(rowI))
      .every((v) => v != null && !isNaN(v)), true);
  });

  test('Identity-longReference', async () => {
    const scoresCol = await grok.functions.call('Bio:sequenceIdentityScoring',
      {table: table, macromolecule: seqCol, reference: longReference}) as DG.Column<number>;
    expect(wu.count(0).take(scoresCol.length).map((rowI) => scoresCol.get(rowI))
      .every((v) => v != null && !isNaN(v)), true);
  });

  test('Similarity', async () => {
    const scoresCol = await grok.functions.call('Bio:sequenceSimilarityScoring',
      {table: table, macromolecule: seqCol, reference: reference}) as DG.Column<number>;
    for (let i = 0; i < scoresCol.length; i++) {
      const resScore = scoresCol.get(i)!;
      const tgtScore = table.get(expectedSimilarity, i);
      expectFloat(resScore, tgtScore, 0.01,
        `Wrong similarity score for sequence at position ${i}`);
    }
  });

  test('seqIdentity', async () => {
    const ref = 'MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW';
    expectFloat(await grok.functions.call('Bio:seqIdentity', {seq: ref, ref: ref}), 1, 0.001);
    const changed: number = await grok.functions.call('Bio:seqIdentity',
      {seq: 'MDYKETLLMPKTAAAAAAAANKEPQIQEKW', ref: ref});
    expect(changed > 0.1 && changed < 0.99, true, `identity of a changed sequence: ${changed}`);
  });

  test('seqIdentity-emptySeq', async () => {
    const res = await grok.functions.call('Bio:seqIdentity',
      {seq: '', ref: 'PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0'});
    expect(res == null, true, `identity of an empty sequence: ${res}`);
  });

  test('sequenceAlignment', async () => {
    const seq1 = 'MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW';
    const cases: [string, string, string, number][] = [
      ['Global alignment', 'BLOSUM62', 'MIEVFLFGIVLGLIPITLAGLFVTAYLQYRRGDQLDL', 37],
      ['Local alignment', 'BLOSUM45', 'AAAAKETLLMPKTDFPAAAA', 12],
    ];
    for (const [alignType, alignTable, seq2, minLength] of cases) {
      const res = await grok.functions.call('Bio:sequenceAlignment',
        {alignType: alignType, alignTable: alignTable, gap: -10, seq1: seq1, seq2: seq2});
      const length = Math.min(res.seq1.length, res.seq2.length);
      expect(length >= minLength, true,
        `${alignType} with ${alignTable}: ${length} aligned positions, expected at least ${minLength}`);
    }
  });
});
