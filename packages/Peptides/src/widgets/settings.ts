import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import * as type from '../utils/types';
import * as C from '../utils/constants';
import {PeptidesModel, VIEWER_TYPE} from '../model';

import $ from 'cash-dom';
import wu from 'wu';
import {getTreeHelperInstance} from '../package';
import {MmDistanceFunctionsNames as distFNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

type PaneInputs = { [paneName: string]: DG.InputBase[] };
type SettingsElements = { dialog: DG.Dialog, accordion: DG.Accordion, inputs: PaneInputs };

export enum SETTINGS_PANES {
  GENERAL = 'General',
  VIEWERS = 'Viewers',
  COLUMNS = 'Columns',
  SEQUENCE_SPACE = 'Sequence space',
  MCL = 'MCL',
}

export enum GENERAL_INPUTS {
  ACTIVITY = 'Activity',
  ACTIVITY_SCALING = 'Activity scaling',
}

export enum VIEWERS_INPUTS {
  DENDROGRAM = VIEWER_TYPE.DENDROGRAM,
}

export enum COLUMNS_INPUTS {
  IS_INCLUDED = '',
  AGGREGATION = 'Aggregation',
}
export enum SEQUENCE_SPACE_INPUTS {
  DISTANCE_FUNCTION = 'Distance function',
  GAP_OPEN = 'Gap open penalty',
  GAP_EXTEND = 'Gap extend penalty',
  CLUSTER_EMBEDDINGS = 'Cluster embeddings',
  EPSILON = 'Epsilon',
  MIN_PTS = 'Minimum points',
  FINGERPRINT_TYPE = 'Fingerprint type',
}

export enum MCL_INPUTS {
  DISTANCE_FUNCTION = 'Distance function',
  GAP_OPEN = 'Gap open penalty',
  GAP_EXTEND = 'Gap extend penalty',
  FINGERPRINT_TYPE = 'Fingerprint type',
  THRESHOLD = 'Similarity threshold',
  MAX_ITERATIONS = 'Max iterations',
}


export const PANES_INPUTS = {
  [SETTINGS_PANES.GENERAL]: GENERAL_INPUTS,
  [SETTINGS_PANES.VIEWERS]: VIEWERS_INPUTS,
  [SETTINGS_PANES.COLUMNS]: COLUMNS_INPUTS,
  [SETTINGS_PANES.SEQUENCE_SPACE]: SEQUENCE_SPACE_INPUTS,
  [SETTINGS_PANES.MCL]: MCL_INPUTS,
};

/**
 * Creates settings dialog for peptides analysis.
 * @param model - Peptides analysis model.
 * @return - Settings dialog elements.
 */
export function getSettingsDialog(model: PeptidesModel): SettingsElements {
  if (model.settings == null)
    grok.log.error('PeptidesError: Settings are not initialized');

  const accordion = ui.accordion();
  const settings = model.settings;
  const currentScaling = settings?.activityScaling ?? C.SCALING_METHODS.NONE;
  const currentColumns = settings?.columns ?? {};

  const result: type.PartialPeptidesSettings = {};
  const inputs: PaneInputs = {};
  const seqSpaceParams = settings?.sequenceSpaceParams ?? new type.SequenceSpaceParams();
  const mclParams = settings?.mclSettings ?? new type.MCLSettings();
  // General pane options
  const activityCol = ui.columnInput(GENERAL_INPUTS.ACTIVITY, model.df,
    model.df.getCol(model.settings!.activityColumnName!), () => result.activityColumnName = activityCol.value!.name,
    {
      filter: (col: DG.Column) => (col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) &&
        col.name !== C.COLUMNS_NAMES.ACTIVITY && col.stats.missingValueCount === 0,
    });
  activityCol.setTooltip('Numeric activity column');
  const activityScaling =
    ui.choiceInput(GENERAL_INPUTS.ACTIVITY_SCALING, currentScaling, Object.values(C.SCALING_METHODS),
      () => result.activityScaling = activityScaling.value as C.SCALING_METHODS) as DG.InputBase<C.SCALING_METHODS>;
  activityScaling.setTooltip('Activity column transformation method');

  accordion.addPane(SETTINGS_PANES.GENERAL, () => ui.inputs([activityCol, activityScaling]), true);
  inputs[SETTINGS_PANES.GENERAL] = [activityCol, activityScaling];

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
  const showSeqSpace = ui.boolInput('Sequence space', !!settings?.showSequenceSpace, () => {
    result.showSequenceSpace = showSeqSpace.value ?? undefined;
    if (showSeqSpace.value) {
      seqSpacePane.root.style.display = 'flex';
      if (!settings?.showSequenceSpace)
        result.sequenceSpaceParams = seqSpaceParams;
    } else {
      seqSpacePane.root.style.display = 'none';
      delete result.sequenceSpaceParams;
    }
    if (result.showSequenceSpace === settings?.showSequenceSpace)
      delete result.showSequenceSpace;
  });
  dendrogram.setTooltip('Show dendrogram viewer');
  dendrogram.enabled = getTreeHelperInstance() !== null;

  accordion.addPane(SETTINGS_PANES.VIEWERS, () => ui.inputs([dendrogram, showSeqSpace]), true);
  inputs[SETTINGS_PANES.VIEWERS] = [dendrogram];

  // Columns to include pane options
  const inputsRows: HTMLElement[] = [];
  const includedColumnsInputs: DG.InputBase[] = [];
  for (const col of model.df.columns.numerical) {
    const colName = col.name;
    if (colName === settings!.activityColumnName || colName === C.COLUMNS_NAMES.ACTIVITY)
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
    isIncludedInput.setTooltip('Include aggregated column value in tooltips, Logo Summary Table and ' +
      'Distribution panel');

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

  // Sequence space pane options
  const modifiedSeqSpaceParams: Partial<type.SequenceSpaceParams> = {};
  function onSeqSpaceParamsChange(fieldName: keyof type.SequenceSpaceParams, value: any): void {
    correctSeqSpaceInputs();
    if (value === null || value === undefined || value === '')
      return;
    modifiedSeqSpaceParams[fieldName] = value;
    let isAllSame = true;
    for (const [key, val] of Object.entries(modifiedSeqSpaceParams)) {
      if (val !== seqSpaceParams[key as keyof type.SequenceSpaceParams]) {
        isAllSame = false;
        break;
      }
    }
    if (isAllSame)
      delete result.sequenceSpaceParams;
    else
      result.sequenceSpaceParams = {...seqSpaceParams, ...modifiedSeqSpaceParams};
  }

  function toggleInputs(nwInputs: DG.InputBase[], condition: boolean): void {
    nwInputs.forEach((input) => {
      if (condition)
        input.root.style.display = 'flex';
      else
        input.root.style.display = 'none';
    });
  }
  // SEQ SPACE INPUTS
  const distanceFunctionInput = ui.choiceInput(SEQUENCE_SPACE_INPUTS.DISTANCE_FUNCTION, seqSpaceParams.distanceF,
    [distFNames.NEEDLEMANN_WUNSCH, distFNames.HAMMING, distFNames.LEVENSHTEIN, distFNames.MONOMER_CHEMICAL_DISTANCE],
    () => onSeqSpaceParamsChange('distanceF', distanceFunctionInput.value));
  distanceFunctionInput.setTooltip('Distance function');
  const gapOpenInput = ui.floatInput(SEQUENCE_SPACE_INPUTS.GAP_OPEN, seqSpaceParams.gapOpen,
    () => onSeqSpaceParamsChange('gapOpen', gapOpenInput.value));
  const gapExtendInput = ui.floatInput(SEQUENCE_SPACE_INPUTS.GAP_EXTEND, seqSpaceParams.gapExtend,
    () => onSeqSpaceParamsChange('gapExtend', gapExtendInput.value));
  const clusterEmbeddingsInput =
    ui.boolInput(SEQUENCE_SPACE_INPUTS.CLUSTER_EMBEDDINGS, seqSpaceParams.clusterEmbeddings ?? false,
      () => onSeqSpaceParamsChange('clusterEmbeddings', clusterEmbeddingsInput.value));
  clusterEmbeddingsInput.setTooltip('Cluster embeddings using DBSCAN algorithm');
  const epsilonInput = ui.floatInput(SEQUENCE_SPACE_INPUTS.EPSILON, seqSpaceParams.epsilon,
    () => onSeqSpaceParamsChange('epsilon', epsilonInput.value));
  epsilonInput.setTooltip(
    'Epsilon parameter for DBSCAN. Minimum distance between two points to be considered as a cluster');
  const minPtsInput = ui.intInput(SEQUENCE_SPACE_INPUTS.MIN_PTS, seqSpaceParams.minPts,
    () => onSeqSpaceParamsChange('minPts', minPtsInput.value));
  minPtsInput.setTooltip('Minimum number of points in a cluster');
  const fingerprintTypesInput = ui.choiceInput('Fingerprint type', seqSpaceParams.fingerprintType,
    ['Morgan', 'RDKit', 'Pattern', 'AtomPair', 'MACCS', 'TopologicalTorsion'],
    () => onSeqSpaceParamsChange('fingerprintType', fingerprintTypesInput.value));
  function correctSeqSpaceInputs(): void {
    toggleInputs([gapOpenInput, gapExtendInput], distanceFunctionInput.value === distFNames.NEEDLEMANN_WUNSCH);
    toggleInputs([epsilonInput, minPtsInput], clusterEmbeddingsInput.value === true);
    toggleInputs([fingerprintTypesInput],
      distanceFunctionInput.value === distFNames.MONOMER_CHEMICAL_DISTANCE ||
      distanceFunctionInput.value === distFNames.NEEDLEMANN_WUNSCH);
  }
  correctSeqSpaceInputs();
  // END OF SEQ SPACE INPUTS

  //MCL INPUTS

  const modifiedMCLParams: Partial<type.MCLSettings> = {};
  function onMCLParamsChange(fieldName: keyof type.MCLSettings, value: any): void {
    correctMCLInputs();
    //correctSeqSpaceInputs();
    if (value === null || value === undefined || value === '')
      return;
    modifiedMCLParams[fieldName] = value;
    let isAllSame = true;
    for (const [key, val] of Object.entries(modifiedMCLParams)) {
      if (val !== mclParams[key as keyof type.MCLSettings]) {
        isAllSame = false;
        break;
      }
    }
    if (isAllSame)
      delete result.mclSettings;
    else
      result.mclSettings = {...mclParams, ...modifiedMCLParams};
  }

  function correctMCLInputs(): void {
    toggleInputs([mclGapOpenInput, mclGapExtendInput], mclDistanceFunctionInput.value === distFNames.NEEDLEMANN_WUNSCH);
    toggleInputs([mclFingerprintTypesInput],
      mclDistanceFunctionInput.value === distFNames.MONOMER_CHEMICAL_DISTANCE ||
      mclDistanceFunctionInput.value === distFNames.NEEDLEMANN_WUNSCH);
  }

  const mclDistanceFunctionInput = ui.choiceInput(SEQUENCE_SPACE_INPUTS.DISTANCE_FUNCTION, mclParams.distanceF,
    [distFNames.NEEDLEMANN_WUNSCH, distFNames.MONOMER_CHEMICAL_DISTANCE, distFNames.HAMMING, distFNames.LEVENSHTEIN],
    () => onMCLParamsChange('distanceF', mclDistanceFunctionInput.value));
  const mclGapOpenInput = ui.floatInput(SEQUENCE_SPACE_INPUTS.GAP_OPEN, mclParams.gapOpen,
    () => onMCLParamsChange('gapOpen', mclGapOpenInput.value));
  const mclGapExtendInput = ui.floatInput(SEQUENCE_SPACE_INPUTS.GAP_EXTEND, mclParams.gapExtend,
    () => onMCLParamsChange('gapExtend', mclGapExtendInput.value));
  const mclFingerprintTypesInput = ui.choiceInput('Fingerprint type', mclParams.fingerprintType,
    ['Morgan', 'RDKit', 'Pattern', 'AtomPair', 'MACCS', 'TopologicalTorsion'],
    () => onMCLParamsChange('fingerprintType', mclFingerprintTypesInput.value));
  const mclThresholdInput = ui.intInput('Similarity threshold', mclParams.threshold ?? 80,
    () => onMCLParamsChange('threshold', mclThresholdInput.value));
  const mclMaxIterationsInput = ui.intInput('Max iterations', mclParams.maxIterations ?? 5,
    () => onMCLParamsChange('maxIterations', mclMaxIterationsInput.value));
  correctMCLInputs();

  const mclInputs = [mclThresholdInput, mclDistanceFunctionInput, mclFingerprintTypesInput,
    mclGapOpenInput, mclGapExtendInput, mclMaxIterationsInput];

  accordion.addPane(SETTINGS_PANES.MCL, () => ui.inputs(mclInputs), true);
  inputs[SETTINGS_PANES.MCL] = mclInputs;
  // END OF MCL INPUTS

  const seqSpaceInputs = [distanceFunctionInput, fingerprintTypesInput, gapOpenInput,
    gapExtendInput, clusterEmbeddingsInput, epsilonInput, minPtsInput];
  const seqSpacePane = accordion.addPane(SETTINGS_PANES.SEQUENCE_SPACE, () => ui.inputs(seqSpaceInputs), true);
  inputs[SETTINGS_PANES.SEQUENCE_SPACE] = seqSpaceInputs;
  showSeqSpace.fireChanged();
  const dialog = ui.dialog('Peptides settings').add(accordion);
  dialog.root.style.width = '400px';
  dialog.onOK(() => {
    model.settings = result;
  });
  dialog.show();

  return {dialog, accordion, inputs};
}
