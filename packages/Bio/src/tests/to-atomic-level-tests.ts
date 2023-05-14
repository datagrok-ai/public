/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, category, test, expectArray} from '@datagrok-libraries/utils/src/test';

import {toAtomicLevel} from '../package';

const appPath = 'System:AppData/Bio';
const fileSource = new DG.FileSource(appPath);

const testNames: {[k: string]: string} = {
  PT: 'peptides fasta',
  DNA: 'dna fasta',
  MSA: 'msa separator',
};

const inputPath: {[k: string]: string} = {
  PT: 'tests/to-atomic-level-peptides-fasta-input.csv',
  DNA: 'tests/to-atomic-level-dna-fasta-input.csv',
  MSA: 'tests/to-atomic-level-msa-separator-input.csv',
}

const ouptputPath: {[k: string]: string} = {
  PT: 'tests/to-atomic-level-peptides-output.csv',
  DNA: 'tests/to-atomic-level-dna-output.csv',
  MSA: 'tests/to-atomic-level-msa-output.csv',
}

const inputColName = 'sequence';
const outputColName = 'molfile(sequence)';

category('to atomic level', async () => {
  const sourceDf: {[key: string]: DG.DataFrame} = {};
  const tragetDf: {[key: string]: DG.DataFrame} = {};

  before(async () => {
    for (const key in testNames) {
      sourceDf[key] = await fileSource.readCsv(inputPath[key]);
      await grok.data.detectSemanticTypes(sourceDf[key]);
      tragetDf[key] = await fileSource.readCsv(ouptputPath[key]);
    }
  });

  async function getTestResult(source: DG.DataFrame, target: DG.DataFrame): Promise<void> {
    const inputCol = source.getCol(inputColName);
    await toAtomicLevel(source, inputCol);
    const obtainedCol = source.getCol(outputColName);
    const expectedCol = target.getCol(outputColName);
    const obtainedArray = [...obtainedCol.values()];
    const expectedArray = [...expectedCol.values()];
    expectArray(obtainedArray, expectedArray);
  }

  for (const key in testNames) {
    test(`${testNames[key]}`, async () => {
      await getTestResult(sourceDf[key], tragetDf[key]);
    });
  }
});
