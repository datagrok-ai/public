import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as type from '../utils/types';
import * as C from '../utils/constants';
import {PeptidesModel, VIEWER_TYPE} from '../model';

import $ from 'cash-dom';
import wu from 'wu';

//TODO: show sliderInput values
export function getSettingsDialog(model: PeptidesModel): DG.Dialog {
  const accordion = ui.accordion();
  const settings = model.settings;
  const result: type.PeptidesSettings = {columns: {}};

  // General pane options
  const activityScaling = ui.choiceInput('Activity scaling', settings.scaling ?? 'none', ['none', 'lg', '-lg'],
    () => result.scaling = activityScaling.value! as type.ScalingMethods);
  const bidirectionalAnalysis = ui.boolInput('Bidirectional analysis', settings.isBidirectional ?? false,
    () => result.isBidirectional = bidirectionalAnalysis.value!);

  accordion.addPane('General', () => ui.inputs([activityScaling, bidirectionalAnalysis]), true);

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
    () => result.showDendrogram = dendrogram.value!);

  accordion.addPane('Viewers',
    () => ui.inputs([dendrogram]), true);

  // Mutation Cliffs pane options
  const maxMutations = ui.sliderInput('Max mutations', settings.maxMutations ?? 1, 0, 50, () => {
    const val = Math.round(maxMutations.value!);
    $(maxMutations.root).find('label.ui-input-description').remove();
    result.maxMutations = val;
    maxMutations.addPostfix(val.toString());
  });
  maxMutations.addPostfix((settings.maxMutations ?? 1).toString());
  const minActivityDelta = ui.sliderInput('Min activity delta', settings.minActivityDelta ?? 0, 0, 100, () => {
    const val = minActivityDelta.value!.toFixed(3);
    result.minActivityDelta = parseFloat(val);
    $(minActivityDelta.root).find('label.ui-input-description').remove();
    minActivityDelta.addPostfix(val);
  });
  minActivityDelta.addPostfix((settings.minActivityDelta ?? 0).toString());
  accordion.addPane('Mutation Cliffs', () => ui.inputs([maxMutations, minActivityDelta]), true);

  // Columns to include pane options
  const inputsRows: HTMLElement[] = [];
  for (const col of model.df.columns.numerical) {
    const colName = col.name;
    if (colName == settings.activityColumnName || colName == C.COLUMNS_NAMES.ACTIVITY_SCALED)
      continue;

    const isIncludedInput = ui.boolInput('', typeof (settings.columns ?? {})[colName] !== 'undefined',
      () => {
        if (isIncludedInput.value)
          result.columns![colName] = aggregationInput.value;
        else
          delete result.columns![colName];
      }) as DG.InputBase<boolean>;
    const aggregationInput = ui.choiceInput('Aggregation', (settings.columns ?? {})[colName] ?? 'avg',
      Object.values(DG.STATS), () => {
        if (isIncludedInput.value)
        result.columns![colName] = aggregationInput.value;
        else
          delete result.columns![col.name];
      }) as DG.InputBase<DG.AggregationType>;
    $(aggregationInput.root).find('label').css('width', 'auto');
    const inputsRow = ui.inputsRow(col.name, [isIncludedInput, aggregationInput]);
    $(inputsRow).find('div.ui-div').css('display', 'inline-flex');
    inputsRows.push(inputsRow);
  }
  if (inputsRows.length != 0)
    accordion.addPane('Columns to include', () => ui.divV(inputsRows), false);

  const dialog = ui.dialog('Peptides settings').add(accordion);
  dialog.root.style.width = '400px';
  dialog.onOK(() => model.settings = result);
  dialog.show();

  return dialog.show();
}
