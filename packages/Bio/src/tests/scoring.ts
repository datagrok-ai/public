import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test, expectFloat, before, after} from '@datagrok-libraries/utils/src/test';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

import {getMonomerLibHelper, sequenceIdentityScoring, sequenceSimilarityScoring} from '../package';
import {
  getUserLibSettings, setUserLibSettings, setUserLibSettingsForTests
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

category('Scoring', () => {
  const sequence = 'sequence';
  const expectedSimilarity = 'expected_similarity';
  const expectedIdentity = 'expected_identity';
  /* eslint-disable max-len */
  const table = DG.DataFrame.fromCsv(`${sequence},${expectedSimilarity},${expectedIdentity}
PEPTIDE1{Aca.Orn.gGlu.Pqa.D-His_1Bn.dH.hHis.4Abz.D-Tic.D-Dap.Y.Iva.meS.F.P.F.D-1Nal}$$$$,1.0,1.0
PEPTIDE1{Iva.Gly_allyl.gGlu.Pqa.D-Dip.dH.hHis.4Abz.D-aHyp.D-Dap.Y.Iva.I.Tyr_26diMe.P.Asu.meC}$$$$,0.68,0.53
PEPTIDE1{[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal].[1Nal]}$$$$V2.0,0.34,0.0`
  );
  /* eslint-enable max-len */
  const seqCol: DG.Column<string> = table.getCol(sequence);
  seqCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
  seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
  const reference = seqCol.get(0)!;

  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;

  before(async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    await setUserLibSettingsForTests();
    await monomerLibHelper.loadLibraries(true); // load default libraries
  });

  after(async () => {
    // UserDataStorage.put() replaces existing data
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadLibraries(true); // load user settings libraries
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
