import * as DG from 'datagrok-api/dg';

import {category, test, expectFloat, before} from '@datagrok-libraries/utils/src/test';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {sequenceIdentityScoring, sequenceSimilarityScoring} from '../package';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

category('Scoring', () => {
  const sequence = 'sequence';
  const expectedSimilarity = 'expected_similarity';
  const expectedIdentity = 'expected_identity';
  const table = DG.DataFrame.fromCsv(`${sequence},${expectedSimilarity},${expectedIdentity}
  PEPTIDE1{Aca.Orn.gGlu.Pqa.D-His_1Bn.dH.hHis.4Abz.D-Tic.D-Dap.Y.Iva.meS.F.P.F.D-1Nal}$$$$,1.0,1.0
  PEPTIDE1{Iva.Gly_allyl.gGlu.Pqa.D-Dip.dH.hHis.4Abz.D-aHyp.D-Dap.Y.Iva.I.Tyr_26diMe.P.Asu.meC}$$$$,0.68,0.53
  PEPTIDE1{[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal]}$$$$V2.0,0.34,0.0
  `);
  const seqCol: DG.Column<string> = table.getCol(sequence);
  seqCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
  seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
  const reference = seqCol.get(0)!;

  before(async () => {
    const monomerLibHelper = await getMonomerLibHelper();
    await monomerLibHelper.loadLibraries(true);
  });

  test('Identity', async () => {
    const scoresCol = await sequenceIdentityScoring(table, seqCol, reference);
    for (let i = 0; i < scoresCol.length; i++) {
      expectFloat(scoresCol.get(i)!, table.get(expectedIdentity, i), 0.01,
        `Wrong identity score for sequence at position ${i}`);
    }
  });

  test('Similarity', async () => {
    const scoresCol = await sequenceSimilarityScoring(table, seqCol, reference);
    for (let i = 0; i < scoresCol.length; i++) {
      expectFloat(scoresCol.get(i)!, table.get(expectedSimilarity, i), 0.01,
        `Wrong similarity score for sequence at position ${i}`);
    }
  });
});
