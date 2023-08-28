import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as type from '../utils/types';
import * as C from '../utils/constants';
import {PeptidesModel, VIEWER_TYPE} from '../model';

import $ from 'cash-dom';
import wu from 'wu';
import {getTreeHelperInstance} from '../package';

type PaneInputs = {[paneName: string]: DG.InputBase[]};
type SettingsElements = {dialog: DG.Dialog, accordion: DG.Accordion, inputs: PaneInputs};

export enum SETTINGS_PANES {
  GENERAL = 'General',
  VIEWERS = 'Viewers',
  MUTATION_CLIFFS = 'Mutation Cliffs',
  COLUMNS = 'Columns',
};

export enum GENERAL_INPUTS {
  ACTIVITY = 'Activity',
  ACTIVITY_SCALING = 'Activity scaling',
  BIDIRECTIONAL_ANALYSIS = 'Bidirectional analysis',
}

export enum VIEWERS_INPUTS {
  DENDROGRAM = VIEWER_TYPE.DENDROGRAM,
}

export enum MUTATION_CLIFFS_INPUTS {
  MAX_MUTATIONS = 'Max mutations',
  MIN_ACTIVITY_DELTA = 'Min activity delta',
}

export enum COLUMNS_INPUTS {
  IS_INCLUDED = '',
  AGGREGATION = 'Aggregation',
}

export const PANES_INPUTS = {
  [SETTINGS_PANES.GENERAL]: GENERAL_INPUTS,
  [SETTINGS_PANES.VIEWERS]: VIEWERS_INPUTS,
  [SETTINGS_PANES.MUTATION_CLIFFS]: MUTATION_CLIFFS_INPUTS,
  [SETTINGS_PANES.COLUMNS]: COLUMNS_INPUTS,
};

//TODO: show sliderInput values
export function getSettingsDialog(model: PeptidesModel): SettingsElements {
  const accordion = ui.accordion();
  const settings = model.settings;
  const currentScaling = settings.scaling ?? C.SCALING_METHODS.NONE;
  const currentBidirectional = settings.isBidirectional ?? false;
  const currentMaxMutations = settings.maxMutations ?? 1;
  const currentMinActivityDelta = settings.minActivityDelta ?? 0;
  const currentColumns = settings.columns ?? {};

  const result: type.PeptidesSettings = {};
  const inputs: PaneInputs = {};

  // General pane options
  const activityCol = ui.columnInput(GENERAL_INPUTS.ACTIVITY, model.df,
    model.df.getCol(model.settings.activityColumnName!), () => result.activityColumnName = activityCol.value!.name,
    {filter: (col: DG.Column) => (col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) && col.name !== C.COLUMNS_NAMES.ACTIVITY_SCALED});
  activityCol.setTooltip('Numeric activity column');
  const activityScaling =
    ui.choiceInput(GENERAL_INPUTS.ACTIVITY_SCALING, currentScaling, Object.values(C.SCALING_METHODS),
      () => result.scaling = activityScaling.value as C.SCALING_METHODS) as DG.InputBase<C.SCALING_METHODS>;
  activityScaling.setTooltip('Activity column transformation method');
  const bidirectionalAnalysis = ui.boolInput(GENERAL_INPUTS.BIDIRECTIONAL_ANALYSIS, currentBidirectional,
    () => result.isBidirectional = bidirectionalAnalysis.value) as DG.InputBase<boolean>;
  bidirectionalAnalysis.setTooltip('Distinguish between positive and negative mean activity difference in ' +
    'Monomer-Position and Most Potent Residues viewers');

  accordion.addPane(SETTINGS_PANES.GENERAL, () => ui.inputs([activityCol, activityScaling, bidirectionalAnalysis]), true);
  inputs[SETTINGS_PANES.GENERAL] = [activityCol, activityScaling, bidirectionalAnalysis];

  // Viewers pane options
  /* FIXME: combinations of adding and deleting viewers are not working properly
  const isMPEnabled = wu(model.analysisView.viewers).some((v) => v.type === VIEWER_TYPE.MONOMER_POSITION);
  const monomerPosition = ui.boolInput(VIEWER_TYPE.MONOMER_POSITION, isMPEnabled ?? false,
    () => result.showMostPotentResidues = monomerPosition.value!);
  const isMPREnabled = wu(model.analysisView.viewers).some((v) => v.type === VIEWER_TYPE.MOST_POTENT_RESIDUES);
  const mostPotentResidues = ui.boolInput(VIEWER_TYPE.MOST_POTENT_RESIDUES, isMPREnabled ?? false,
    () => result.showMonomerPosition = mostPotentResidues.value!);
  const isLSTEnabled = wu(model.analysisView.viewers).some((v) => v.type === VIEWER_TYPE.LOGO_SUMMARY_TABLE);
  const logoSummaryTable = ui.boolInput(VIEWER_TYPE.LOGO_SUMMARY_TABLE, isLSTEnabled ?? false,
    () => result.showLogoSummaryTable = logoSummaryTable.value!);
  logoSummaryTable.enabled = typeof settings.clustersColumnName !== 'undefined';
  */
  const isDendrogramEnabled = wu(model.analysisView.viewers).some((v) => v.type === VIEWER_TYPE.DENDROGRAM);
  const dendrogram = ui.boolInput(VIEWER_TYPE.DENDROGRAM, isDendrogramEnabled ?? false,
    () => result.showDendrogram = dendrogram.value) as DG.InputBase<boolean>;
  dendrogram.setTooltip('Show dendrogram viewer');
  dendrogram.enabled = getTreeHelperInstance() !== null;

  accordion.addPane(SETTINGS_PANES.VIEWERS, () => ui.inputs([dendrogram]), true);
  inputs[SETTINGS_PANES.VIEWERS] = [dendrogram];

  // Mutation Cliffs pane options
  const maxMutations = ui.sliderInput(MUTATION_CLIFFS_INPUTS.MAX_MUTATIONS, currentMaxMutations, 1, 50, () => {
    const val = Math.round(maxMutations.value);
    $(maxMutations.root).find('label.ui-input-description').remove();
    result.maxMutations = val;
    maxMutations.addPostfix(val.toString());
  }) as DG.InputBase<number>;
  maxMutations.setTooltip('Maximum number of mutations between reference and mutated sequences');
  maxMutations.addPostfix((settings.maxMutations ?? 1).toString());
  const minActivityDelta = ui.sliderInput(MUTATION_CLIFFS_INPUTS.MIN_ACTIVITY_DELTA, currentMinActivityDelta, 0,
    100, () => {
      const val = minActivityDelta.value.toFixed(3);
      result.minActivityDelta = parseFloat(val);
      $(minActivityDelta.root).find('label.ui-input-description').remove();
      minActivityDelta.addPostfix(val);
    }) as DG.InputBase<number>;
  minActivityDelta.setTooltip('Minimum activity difference between reference and mutated sequences');
  minActivityDelta.addPostfix((settings.minActivityDelta ?? 0).toString());
  accordion.addPane(SETTINGS_PANES.MUTATION_CLIFFS, () => ui.inputs([maxMutations, minActivityDelta]), true);
  inputs[SETTINGS_PANES.MUTATION_CLIFFS] = [maxMutations, minActivityDelta];

  // Columns to include pane options
  const inputsRows: HTMLElement[] = [];
  const includedColumnsInputs: DG.InputBase[] = [];
  for (const col of model.df.columns.numerical) {
    const colName = col.name;
    if (colName === settings.activityColumnName || colName === C.COLUMNS_NAMES.ACTIVITY_SCALED)
      continue;

    const isIncludedInput = ui.boolInput(COLUMNS_INPUTS.IS_INCLUDED, typeof (currentColumns)[colName] !== 'undefined',
      () => {
        result.columns ??= {};
        if (isIncludedInput.value)
          result.columns[colName] = aggregationInput.value;
        else {
          delete result.columns[colName];
          if (Object.keys(result.columns).length === Object.keys(currentColumns).length)
            delete result.columns;
        }
      }) as DG.InputBase<boolean>;
    isIncludedInput.setTooltip('Include aggregated column value in tooltips, Logo Summary Table and Distribution panel');

    const aggregationInput = ui.choiceInput(COLUMNS_INPUTS.AGGREGATION, (currentColumns)[colName] ?? DG.AGG.AVG,
      Object.values(DG.STATS), () => {
        result.columns ??= {};
        if (isIncludedInput.value)
          result.columns[colName] = aggregationInput.value;
        else {
          delete result.columns[col.name];
          if (Object.keys(result.columns).length === Object.keys(currentColumns).length)
            delete result.columns;
        }
      }) as DG.InputBase<DG.AggregationType>;
    aggregationInput.setTooltip('Aggregation method');
    $(aggregationInput.root).find('label').css('width', 'auto');
    const inputsRow = ui.inputsRow(col.name, [isIncludedInput, aggregationInput]);
    includedColumnsInputs.push(...[isIncludedInput, aggregationInput]);
    $(inputsRow).find('div.ui-div').css('display', 'inline-flex');
    inputsRows.push(inputsRow);
  }
  if (inputsRows.length !== 0) {
    accordion.addPane(SETTINGS_PANES.COLUMNS, () => ui.divV(inputsRows), false);
    inputs[SETTINGS_PANES.COLUMNS] = includedColumnsInputs;
  }

  const dialog = ui.dialog('Peptides settings').add(accordion);
  dialog.root.style.width = '400px';
  dialog.onOK(() => model.settings = result);
  dialog.show();

  return {dialog, accordion, inputs};
}
