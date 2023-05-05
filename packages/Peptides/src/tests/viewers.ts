import * as DG from 'datagrok-api/dg';

import {before, category, expect, test, testViewer} from '@datagrok-libraries/utils/src/test';
import {aligned1} from './test-data';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import {_package} from '../package-test';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {scaleActivity} from '../utils/misc';
import {startAnalysis} from '../widgets/peptides';
import {MONOMER_POSITION_MODE, MonomerPosition} from '../viewers/sar-viewer';
import {SCALING_METHODS} from '../utils/constants';


category('Viewers: Basic', () => {
  const df = DG.DataFrame.fromCsv(aligned1);
  const viewers = DG.Func.find({package: 'Peptides', tags: ['viewer']}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      await testViewer(v, df.clone(), true);
    }, {skipReason: 'GROK-11534'});
  }
});

category('Viewers: Monomer-Position', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;

  before(async () => {
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol('activity');
    sequenceCol = df.getCol('sequence');
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, SCALING_METHODS.NONE);
    clusterCol = df.getCol('cluster');
    const tempModel = await startAnalysis(
      activityCol, sequenceCol, clusterCol, df, scaledActivityCol, SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
  });

  test('Tooltip', async () => {

  }, {skipReason: 'Not implemented yet'});

  test('Modes', async () => {
    const mpViewer = model.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition | null;
    if (mpViewer === null)
      throw new Error('Monomer-Position viewer doesn\'t exist');

    expect(mpViewer.mode, MONOMER_POSITION_MODE.MUTATION_CLIFFS,
      `Default Monomer-Position mode is not ${MONOMER_POSITION_MODE.MUTATION_CLIFFS}`);

    mpViewer.mode = MONOMER_POSITION_MODE.INVARIANT_MAP;
    expect(mpViewer.mode, MONOMER_POSITION_MODE.INVARIANT_MAP,
      `Monomer-Position mode is not ${MONOMER_POSITION_MODE.INVARIANT_MAP} after switching`);

    mpViewer.mode = MONOMER_POSITION_MODE.MUTATION_CLIFFS;
    expect(mpViewer.mode, MONOMER_POSITION_MODE.MUTATION_CLIFFS,
      `Monomer-Position mode is not ${MONOMER_POSITION_MODE.MUTATION_CLIFFS} after switching`);
  });
});

category('Viewers: Most Potent Residues', () => {
  test('Tooltip', async () => {

  }, {skipReason: 'Not implemented yet'});
});

category('Viewers: Logo Summary Table', () => {
  test('Properties', async () => {

  }, {skipReason: 'Not implemented yet'});

  test('Tooltip', async () => {

  }, {skipReason: 'Not implemented yet'});
});
