import * as DG from 'datagrok-api/dg';

import {category, test, before, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {PeptidesModel} from '../model';
import {scaleActivity} from '../utils/misc';
import {startAnalysis} from '../widgets/peptides';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SCALING_METHODS} from '../utils/constants';
import {PANES_INPUTS, SETTINGS_PANES, getSettingsDialog} from '../widgets/settings';
import {getDistributionWidget} from '../widgets/distribution';
import {mutationCliffsWidget} from '../widgets/mutation-cliffs';

category('Widgets: Settings', () => {
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
    const tempModel = await startAnalysis(activityCol, sequenceCol, clusterCol, df, scaledActivityCol,
      SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
  });

  test('UI', async () => {
    const settingsElements = getSettingsDialog(model);

    // Check number of panes
    const panes = settingsElements.accordion.panes.map((pane) => pane.name);
    expect(panes.length, 4, `Expected 4 panes, got ${settingsElements.accordion.panes.length}`);
    for (const paneName of Object.values(SETTINGS_PANES))
      expect(panes.includes(paneName), true, `Pane ${paneName} is missing`);

    // Check inputs in each pane
    for (const paneName of Object.values(SETTINGS_PANES)) {
      const paneInputs = settingsElements.inputs[paneName].map((input) => input.caption);
      for (const inputName of Object.values(PANES_INPUTS[paneName]))
        expect(paneInputs.includes(inputName), true, `Input ${inputName} is missing from ${paneName}`);
    }
  });
});

category('Widgets: Distribution panel', () => {
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
    const tempModel = await startAnalysis(activityCol, sequenceCol, clusterCol, df, scaledActivityCol,
      SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
  });

  test('UI', async () => {
    getDistributionWidget(model.df, model);
  });

  test('Split', async () => {

  }, {skipReason: 'Not implemented yet'});
});

category('Widgets: Mutation cliffs', () => {
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
    const tempModel = await startAnalysis(activityCol, sequenceCol, clusterCol, df, scaledActivityCol,
      SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
  });

  test('UI', async () => {
    mutationCliffsWidget(model.df, model);
  });

  test('General', async () => {

  }, {skipReason: 'Not implemented yet'});

  test('Filtering', async () => {

  }, {skipReason: 'Not implemented yet'});
});

category('Widgets: Actions', () => {
  test('New view', async () => {

  }, {skipReason: 'Not implemented yet'});

  test('New cluster', async () => {

  });

  test('Remove cluster', async () => {

  }, {skipReason: 'Not implemented yet'});
});
