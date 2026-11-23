import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {readDataframe} from './utils';
import {_testSequenceSpaceOnNotation, _testSequenceSpaceReturnsResult} from './sequence-space-utils';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BYPASS_LARGE_DATA_WARNING} from '@datagrok-libraries/ml/src/functionEditors/consts';
import {getGPUDevice} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

category('sequenceSpace', async () => {
  let testFastaDf: DG.DataFrame;
  let testFastaTableView: DG.TableView;
  let testHelmWithEmptyRows: DG.DataFrame;
  let testHelmWithEmptyRowsTableView: DG.TableView;

  test('sequenceSpaceOpens', async () => {
    const testData = !DG.Test.isInBenchmark ?
      {fileName: 'tests/100_3_clustests.csv'} :
      {fileName: 'tests/peptides_with_random_motif_1600.csv'};
    testFastaDf = await readDataframe(testData.fileName);
    testFastaTableView = grok.shell.addTableView(testFastaDf);
    await _testSequenceSpaceReturnsResult(testFastaDf, DimReductionMethods.UMAP, 'sequence');
    //grok.shell.closeTable(testFastaDf);
    //testFastaTableView.close();
  }, {benchmark: true});

  test('sequenceSpaceWithEmptyRows', async () => {
    testHelmWithEmptyRows = await readDataframe('tests/100_3_clustests_empty_vals.csv');
    testHelmWithEmptyRowsTableView = grok.shell.addTableView(testHelmWithEmptyRows);
    await _testSequenceSpaceReturnsResult(testHelmWithEmptyRows, DimReductionMethods.UMAP, 'sequence');
    //grok.shell.closeTable(testHelmWithEmptyRows);
    //testHelmWithEmptyRowsTableView.close();
  });

  test('sequenceSpaceHelm', async () => {
    await _testSequenceSpaceOnNotation('samples/HELM.csv', 'HELM', NOTATION.HELM);
  });

  test('sequenceSpaceMsa', async () => {
    await _testSequenceSpaceOnNotation('samples/MSA.csv', 'MSA', NOTATION.SEPARATOR);
  });

  test('sequenceSpace: same seed, same coordinates, on the GPU and on the CPU', async () => {
    const df = await readDataframe('tests/100_3_clustests.csv');
    await grok.data.detectSemanticTypes(df);
    const tv = grok.shell.addTableView(df);
    try {
      const embed = async (randomSeed: string, useWebGPU: boolean) => {
        const before = df.columns.names();
        await grok.functions.call('Bio:sequenceSpaceTopMenu', {
          table: df, molecules: df.getCol('sequence'), methodName: DimReductionMethods.UMAP,
          similarityMetric: MmDistanceFunctionsNames.LEVENSHTEIN, plotEmbeddings: false, clusterEmbeddings: false,
          options: {randomSeed, useWebGPU, [BYPASS_LARGE_DATA_WARNING]: true},
        });
        return JSON.stringify(df.columns.names().filter((name) => !before.includes(name))
          .map((name) => df.col(name)!.toList()));
      };
      for (const useWebGPU of await getGPUDevice() ? [true, false] : [false]) {
        const device = useWebGPU ? 'GPU' : 'CPU';
        const first = await embed('42', useWebGPU);
        expect(await embed('42', useWebGPU), first, `${device}: seed 42 gave other coordinates on the second run`);
        expect(await embed('43', useWebGPU) !== first, true, `${device}: seeds 42 and 43 gave the same coordinates`);
      }
    } finally {tv.close();}
  }, {timeout: 120000});
});
