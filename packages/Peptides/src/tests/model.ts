import * as DG from 'datagrok-api/dg';

import {category, test, before, expect, expectArray, expectFloat} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {PeptidesModel} from '../model';
import {startAnalysis} from '../widgets/peptides';
import {scaleActivity} from '../utils/misc';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {COLUMNS_NAMES, SCALING_METHODS} from '../utils/constants';

category('Model: Settings', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;

  const firstMonomerPair = {monomer: 'N', position: '4', count: 7};
  const secondMonomerPair = {monomer: 'meI', position: '1', count: 10};
  const firstCluster = {name: '0', count: 3};
  const secondCluster = {name: '1', count: 3};

  before(async () => {
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol('activity');
    sequenceCol = df.getCol('sequence');
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, SCALING_METHODS.NONE);
    clusterCol = df.getCol('cluster');
    const tempModel = await startAnalysis(activityCol, sequenceCol, clusterCol, df, scaledActivityCol,
      SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
  });

  test('Activity scaling', async () => {
    const tolerance = 0.0001;
    const origActivityData =
      model.df.getCol(model.settings.activityColumnName!).getRawData();
    const scaledActivity = model.df.getCol(COLUMNS_NAMES.ACTIVITY_SCALED);
    const dfLen = model.df.rowCount;

    // Check 'none' scaling
    let scaledActivityData = scaledActivity.getRawData();
    for (let i = 0; i < dfLen; i++) {
      expectFloat(scaledActivityData[i], origActivityData[i], tolerance,
        `Activity mismatch at row ${i} for scaling method ` +
        `'${SCALING_METHODS.NONE}'`);
    }
    
    // Check 'lg' scaling
    model.settings = {scaling: SCALING_METHODS.LG};
    scaledActivityData = scaledActivity.getRawData();
    for (let i = 0; i < dfLen; i++) {
      expectFloat(scaledActivityData[i], Math.log10(origActivityData[i]),
        tolerance, `Activity mismatch at row ${i} for scaling method ` +
        `'${SCALING_METHODS.LG}'`);
    }

    // Check '-lg' scaling
    model.settings = {scaling: SCALING_METHODS.MINUS_LG};
    scaledActivityData = scaledActivity.getRawData();
    for (let i = 0; i < dfLen; i++) {
      expectFloat(scaledActivityData[i], -Math.log10(origActivityData[i]),
        tolerance, `Activity mismatch at row ${i} for scaling method ` +
        `'${SCALING_METHODS.MINUS_LG}'`);
    }
  });

  test('Bidirectional analysis', async () => {

  }, {skipReason: 'Not implemented yet'});

  test('Mutation Cliffs', async () => {

  }, {skipReason: 'Not implemented yet'});

  test('Include columns', async () => {

  }, {skipReason: 'Not implemented yet'});

  test('Dendrogram', async () => {

  }, {skipReason: 'Not implemented yet'});
});
