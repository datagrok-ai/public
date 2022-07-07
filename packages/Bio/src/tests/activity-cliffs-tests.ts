import {after, before, category, expect, expectFloat, test} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import {createTableView, readDataframe} from './utils';
import {_package} from '../package-test';
import {getEmbeddingColsNames, sequenceSpace} from '../utils/sequence-space';
import {drawTooltip, sequenceGetSimilarities} from '../utils/sequence-activity-cliffs';
import {getActivityCliffs} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';


category('activityCliffs', async () => {
  let actCliffsTableView: DG.TableView;
  let actCliffsDf: DG.DataFrame;

  before(async () => {
    actCliffsTableView = await createTableView('sample_MSA.csv');
    actCliffsDf = await readDataframe('sample_MSA.csv');
  });

  test('activityCliffsOpen', async () => {
    const axesNames = getEmbeddingColsNames(actCliffsDf);
    const units = actCliffsDf.col('MSA')!.tags[DG.TAGS.UNITS];
    const options = {
      'SPE': {cycles: 2000, lambda: 1.0, dlambda: 0.0005},
    };
    const scatterPlot = await getActivityCliffs(
      actCliffsDf,
      actCliffsDf.col('MSA')!,
      axesNames,
      'Activity cliffs',
      actCliffsDf.col('Activity')!,
      50,
      'Levenshtein',
      't-SNE',
      DG.SEMTYPE.MACROMOLECULE,
      units,
      sequenceSpace,
      sequenceGetSimilarities,
      drawTooltip);

    expect(scatterPlot != null, true);

    const cliffsLink = (Array.from(scatterPlot.root.children) as Element[])
      .filter((it) => it.className === 'ui-btn ui-btn-ok');
    expect((cliffsLink[0] as HTMLElement).innerText, '101 cliffs');
  });

  after(async () => {
    actCliffsTableView.close();
  });
});
