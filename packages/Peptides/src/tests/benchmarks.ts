import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test} from '@datagrok-libraries/utils/src/test';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

import {_package} from '../package-test';
import {findMutations} from '../utils/algorithms';
import * as type from '../utils/types';
import {scaleActivity} from '../utils/misc';
import {PeptidesModel} from '../model';
import {startAnalysis} from '../widgets/peptides';
import * as C from '../utils/constants';


const benchmarkDatasetSizes = [5, 50, 100, 200];

category('Benchmarks: Mutation Cliffs', () => {
  for (const size of benchmarkDatasetSizes)
    test(`${size}k sequences`, async () => await mutationCliffsBenchmark(size), {timeout: 60000});
});

category('Benchmarks: Cluster stats', () => {
  for (const size of benchmarkDatasetSizes) {
    test(`${size}k sequences`, async () => {
      if (!DG.Test.isInBenchmark)
        return null;

      const model = await getModel(size);
      DG.time(`Cluster stats benchmark - ${size}`, () => model!.calculateClusterStatistics());
    }, {timeout: 60000});
  }
});

category('Benchmarks: Monomer-Position stats', () => {
  for (const size of benchmarkDatasetSizes) {
    test(`${size}k sequences`, async () => {
      if (!DG.Test.isInBenchmark)
        return null;

      const model = await getModel(size);
      DG.time(`Monomer-Position stats benchmark - ${size}`, () => model!.calculateMonomerPositionStatistics());
    }, {timeout: 60000});
  }
});

category('Benchmarks: Analysis start', () => {
  for (const size of benchmarkDatasetSizes) {
    test(`${size}k sequences`, async () => {
      if (!DG.Test.isInBenchmark)
        return;

      const df = (await _package.files.readBinaryDataFrames(`tests/${size}k.d42`))[0];
      const activityCol = df.getCol('activity');
      const scaledActivityCol = scaleActivity(activityCol, C.SCALING_METHODS.NONE);
      const clustersCol = df.getCol('cluster');
      const sequenceCol = df.getCol('sequence');
      sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
      sequenceCol.setTag(DG.TAGS.UNITS, size === benchmarkDatasetSizes[0] ? NOTATION.HELM : NOTATION.FASTA);

      await DG.timeAsync('Analysis start', async () => {
        const model = await startAnalysis(activityCol, sequenceCol, clustersCol, df, scaledActivityCol, C.SCALING_METHODS.NONE);

        if (model)
          grok.shell.closeTable(model.df);
      });
    });
  }
});

async function mutationCliffsBenchmark(size: number): Promise<void> {
  if (!DG.Test.isInBenchmark)
    return;

  const df = (await _package.files.readBinaryDataFrames(`tests/${size}k.d42`))[0];
  const activityCol: type.RawData = df.getCol('activity').getRawData();
  const monomerCols: type.RawColumn[] = [];
  let i = 1;
  while (df.col(i.toString()) !== null) {
    const col = df.getCol(i.toString());
    monomerCols.push({name: col.name, rawData: col.getRawData(), cat: col.categories});
    ++i;
  }
  DG.time('Mutation Cliffs', () => findMutations(activityCol, monomerCols));
}

async function getModel(size: number): Promise<PeptidesModel | null> {
  if (!DG.Test.isInBenchmark)
    return null;

  const df = (await _package.files.readBinaryDataFrames(`tests/${size}k.d42`))[0];
  const activityCol = df.getCol('activity');
  const scaledActivityCol = scaleActivity(activityCol, C.SCALING_METHODS.NONE);
  const clustersCol = df.getCol('cluster');
  const sequenceCol = df.getCol('sequence');
  sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
  sequenceCol.setTag(DG.TAGS.UNITS, size === benchmarkDatasetSizes[0] ? NOTATION.HELM : NOTATION.FASTA);
  return await startAnalysis(activityCol, sequenceCol, clustersCol, df, scaledActivityCol, C.SCALING_METHODS.NONE);
}
