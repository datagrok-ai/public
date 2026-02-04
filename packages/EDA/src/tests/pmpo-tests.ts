// Tests for Probabilistic MPO (pMPO)
// Reference scores are pre-computed and stored in the 'drugs-props-train-scores.csv' file.
// This scores are computed using the library: https://github.com/Merck/pmpo

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {Pmpo} from '../probabilistic-scoring/prob-scoring';
import {P_VAL_TRES_DEFAULT, Q_CUTOFF_DEFAULT, R2_DEFAULT, SCORES_PATH,
  SOURCE_PATH} from '../probabilistic-scoring/pmpo-defs';

const TIMEOUT = 4000;
const MAD_THRESH = 1E-6;

const DESIRABILITY_COL_NAME = 'CNS';
const DESCRIPTOR_NAMES = ['TPSA', 'TPSA_S', 'HBA', 'HBD', 'MW', 'nAtoms',
  'cLogD_ACD_v15', 'mapKa', 'cLogP_Biobyte', 'mbpKa', 'cLogP_ACD_v15', 'ALogP98'];
const SCORES_NAME = 'Score';
const DRUG = 'Drug';

const SIGMOIDAL = 'Sigmoidal';
const GAUSSIAN = 'Gaussian';
const PMPO_MODES = [SIGMOIDAL, GAUSSIAN];

/** Computes the maximum absolute deviation between pMPO scores in two data frames */
function getScoreMaxDeviation(sourceDrugCol: DG.Column, sourceScores: DG.Column,
  referenceDrugCol: DG.Column, referenceScores: DG.Column): number {
  let mad = 0;

  const sourceDrugList = sourceDrugCol.toList();
  const referenceDrugList = referenceDrugCol.toList();

  const sourceScoresRaw = sourceScores.getRawData();
  const referenceScoresRaw = referenceScores.getRawData();

  sourceDrugList.forEach((name, idx) => {
    const refIdx = referenceDrugList.indexOf(name);

    if (refIdx < 0)
      throw new Error(`Failed to compare pMPO scores: the "${name}" drug is missing in the reference data.`);

    mad = Math.max(mad, Math.abs(sourceScoresRaw[idx] - referenceScoresRaw[refIdx]));
  });

  return mad;
} // getScoreMaxDeviation

category('Probabilistic MPO', () => {
  PMPO_MODES.forEach((refScoreName) => {
    const useSigmoid = (refScoreName == SIGMOIDAL);

    test(refScoreName, async () => {
      let sourceDf: DG.DataFrame | null = null;
      let referenceDf: DG.DataFrame | null = null;
      let desirability: DG.Column | null = null;
      let descriptors: DG.Column[] = [];
      let sourceDrugCol: DG.Column | null = null;
      let referenceDrugCol: DG.Column | null = null;
      let referencePrediction: DG.Column | null = null;
      let mad: number | null = null;

      try {
      // Load data
        sourceDf = await grok.dapi.files.readCsv(SOURCE_PATH);
        referenceDf = await grok.dapi.files.readCsv(SCORES_PATH);

        // Extract training items
        desirability = sourceDf.col(DESIRABILITY_COL_NAME);
        descriptors = sourceDf.columns.byNames(DESCRIPTOR_NAMES);

        if (desirability == null)
          throw new Error();

        // Train pMPO model
        const trainRes = Pmpo.fit(
          sourceDf,
          DG.DataFrame.fromColumns(descriptors).columns,
          desirability,
          P_VAL_TRES_DEFAULT,
          R2_DEFAULT,
          Q_CUTOFF_DEFAULT,
        );

        // Apply pMPO
        const prediction = Pmpo.predict(sourceDf, trainRes.params, useSigmoid, SCORES_NAME);

        // Compare with reference scores
        sourceDrugCol = sourceDf.col(DRUG);
        referenceDrugCol = referenceDf.col(DRUG);
        referencePrediction = referenceDf.col(refScoreName);

        mad = getScoreMaxDeviation(sourceDrugCol!, prediction, referenceDrugCol!, referencePrediction!);

        //console.log(refScoreName, ': max absolute deviation of pMPO scores:', mad);
      } catch (error) {
        grok.shell.error((error as Error).message);
      }

      expect(sourceDf !== null, true, 'Failed to load the source data: ' + SOURCE_PATH);
      expect(referenceDf !== null, true, 'Failed to load the scores data: ' + SCORES_PATH);
      expect(desirability !== null, true, 'Inconsistent source data: no column ' + DESIRABILITY_COL_NAME);
      expect(descriptors.length, DESCRIPTOR_NAMES.length, 'Inconsistent source data: no enough of columns');
      expect(sourceDrugCol !== null, true, 'Inconsistent source data: no column ' + DRUG);
      expect(referenceDrugCol !== null, true, 'Inconsistent reference data: no column ' + DRUG);
      expect(referencePrediction !== null, true, 'Inconsistent reference data: no column ' + SCORES_NAME);
      expect(mad !== null, true, 'Failed to compare pMPO scores with the reference data');
      expect(mad! < MAD_THRESH, true, `Max absolute deviation of pMPO scores exceeds the threshold (${MAD_THRESH})`);
    }, {timeout: TIMEOUT});
  });
});
