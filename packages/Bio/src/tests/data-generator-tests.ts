import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, expect, test, before, expectFloat} from '@datagrok-libraries/utils/src/test';
import {colNames, DistributionParams, generatePeptidesDataFrame} from '@datagrok-libraries/bio/src/utils/data-generator';
import {getSplitterForColumn, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

category('data-generator', () => {
  const naturalMonomers: string[] = Array.from(UnitsHandler.PeptideFastaAlphabet);
  const polymerType = 'PEPTIDE';
  const distParams: DistributionParams = {type: 'normal', mean: 0, std: 1};
  const clusterSequenceLength = [14];
  let monomerLibMonomers: string[];

  before(async () => {
    const bioLib = await grok.functions.call('Bio:getBioLib');
    monomerLibMonomers = bioLib.getMonomerNamesByType(polymerType);
  });

  test('Generate FASTA', async () => {
    return testDataGeneration(NOTATION.FASTA, 10, naturalMonomers, clusterSequenceLength, distParams);
  });

  test('Generate separator', async () => {
    return testDataGeneration(NOTATION.SEPARATOR, 10, monomerLibMonomers, clusterSequenceLength, distParams);
  });

  test('Generate HELM', async () => {
    return testDataGeneration(NOTATION.HELM, 10, monomerLibMonomers, clusterSequenceLength, distParams);
  });
});

async function testDataGeneration(notation: NOTATION, length: number, monomers: string[],
  clusterSequenceLength: number[], distParams: DistributionParams): Promise<void> {
  const df = generatePeptidesDataFrame(notation, length, monomers, clusterSequenceLength, distParams);
  expect(df.rowCount, length);

  grok.shell.addTableView(df);
  await grok.data.detectSemanticTypes(df);

  const sequencesCol = df.getCol(colNames.sequences);
  expect(sequencesCol.getTag(DG.TAGS.UNITS), notation);

  const splitter = getSplitterForColumn(sequencesCol);

  const clustersCol = df.getCol(colNames.clusters);
  for (let rowIdx = 0; rowIdx < length; ++rowIdx) {
    const sequence = splitter(sequencesCol.get(rowIdx));
    const cluster = clustersCol.get(rowIdx);
    expect(sequence.length, cluster);
  }

  const activityCol = df.getCol(colNames.activity);
  expectFloat(activityCol.stats.avg, distParams.mean, distParams.std);
}
