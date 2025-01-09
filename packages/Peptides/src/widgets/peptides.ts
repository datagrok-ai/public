/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../styles.css';
import * as C from '../utils/constants';
import * as type from '../utils/types';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import $ from 'cash-dom';
import {scaleActivity} from '../utils/misc';
import {ALIGNMENT, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ILogoSummaryTable, LogoSummaryTable} from '../viewers/logo-summary';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {MCL_INPUTS} from './settings';

export type DialogParameters = { host: HTMLElement, callback: () => Promise<boolean> };

/**
 * Peptides analysis parameters UI
 * @param df - Dataframe with peptides
 * @param [col] - Peptides column
 * @return - UI host and analysis start callback
 */
export function analyzePeptidesUI(df: DG.DataFrame, col?: DG.Column<string>): DialogParameters {
  const mclOptions = new type.MCLSettings();
  const logoHost = ui.div();
  let seqColInput: DG.InputBase | null = null;
  if (typeof col === 'undefined') {
    // Building UI for starting analysis from dialog (top menu)
    const potentialCol = DG.Utils.firstOrNull(
      df.columns.toList().filter((dfCol) => dfCol.semType === DG.SEMTYPE.MACROMOLECULE));
    if (potentialCol === null)
      throw new Error('Peptides Error: table doesn\'t contain sequence columns');
    else if (potentialCol.stats.missingValueCount !== 0)
      grok.shell.info('Sequences column contains missing values. They will be ignored during analysis');


    seqColInput = ui.input.column('Sequence', {table: df, value: potentialCol, onValueChanged: (value) => {
      $(logoHost).empty().append(ui.wait(async () => {
        const viewer = await df.plot.fromType('WebLogo', {sequenceColumnName: value.name});
        viewer.root.style.setProperty('height', '130px');
        return viewer.root;
      }));
      if (value.stats.missingValueCount !== 0)
        grok.shell.info('Sequences column contains missing values. They will be ignored during analysis');
    }, filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE});
    seqColInput.setTooltip('Macromolecule column in FASTA, HELM or separated format');
  } else if (!(col.getTag(bioTAGS.aligned) === ALIGNMENT.SEQ_MSA) &&
    col.getTag(DG.TAGS.UNITS) !== NOTATION.HELM) {
    return {
      host: ui.label('Peptides analysis only works with aligned sequences'),
      callback: async (): Promise<boolean> => false,
    };
  }

  let funcs = DG.Func.find({package: 'Bio', name: 'webLogoViewer'});
  if (funcs.length === 0) {
    return {
      host: ui.label('Bio package is missing or out of date. Please install the latest version.'),
      callback: async (): Promise<boolean> => false,
    };
  }

  funcs = DG.Func.find({package: 'Bio', name: 'getBioLib'});
  if (funcs.length === 0) {
    return {
      host: ui.label('Bio package is missing or out of date. Please install the latest version.'),
      callback: async (): Promise<boolean> => false,
    };
  }

  // Activity column properties
  let scaledCol: DG.Column<number>;
  const defaultActivityColumn: DG.Column<number> | null = df.col('activity') || df.col('IC50') ||
    DG.Utils.firstOrNull(df.columns.numerical);
  const histogramHost = ui.div([], {id: 'pep-hist-host'});

  const activityScalingMethod = ui.input.choice(
    'Scaling', {value: C.SCALING_METHODS.NONE, items: Object.values(C.SCALING_METHODS),
      onValueChanged: async (value): Promise<void> => {
        scaledCol = scaleActivity(activityColumnChoice.value!, value);

        const hist = DG.DataFrame.fromColumns([scaledCol]).plot.histogram({
          filteringEnabled: false, valueColumnName: C.COLUMNS_NAMES.ACTIVITY, legendVisibility: 'Never', showXAxis: true,
          showColumnSelector: false, showRangeSlider: false, showBinSelector: false,
        });
        histogramHost.lastChild?.remove();
        histogramHost.appendChild(hist.root);
      }}) as DG.InputBase<C.SCALING_METHODS | null>;
  activityScalingMethod.setTooltip('Activity column transformation method');

  const activityScalingMethodState = (): void => {
    activityScalingMethod.enabled = (activityColumnChoice.value ?? false) && activityColumnChoice.value!.stats.min > 0;
    activityScalingMethod.value = C.SCALING_METHODS.NONE;
    if (activityColumnChoice.value!.stats.missingValueCount !== 0)
      grok.shell.info('Activity column contains missing values. They will be ignored during analysis');
  };
  const activityColumnChoice = ui.input.column('Activity', {table: df, value: defaultActivityColumn!,
    onValueChanged: activityScalingMethodState, filter: (col: DG.Column) => col.type === DG.TYPE.INT || col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.QNUM});
  activityColumnChoice.setTooltip('Numerical activity column');
  const clustersColumnChoice = ui.input.column('Clusters', {table: df, onValueChanged: (value) => {
    if (value) {
      generateClustersInput.value = false;
      generateClustersInput.fireChanged();
    }
  }});
  clustersColumnChoice.setTooltip('Optional. Clusters column is used to create Logo Summary Table');
  clustersColumnChoice.nullable = true;
  // clustering input
  const generateClustersInput = ui.input.bool('Generate clusters', {value: true, onValueChanged: (value) => {
    if (value) {
      //@ts-ignore
      clustersColumnChoice.value = null;
      clustersColumnChoice.fireChanged();
    }
  }});
  generateClustersInput
    .setTooltip('Generate clusters column based on MCL embeddings for Logo Summary Table');
  activityColumnChoice.fireChanged();
  activityScalingMethod.fireChanged();
  generateClustersInput.fireChanged();


  const inputsList = [activityColumnChoice, activityScalingMethod, clustersColumnChoice, generateClustersInput];
  if (seqColInput !== null)
    inputsList.splice(0, 0, seqColInput);

  // ### MCL INPUTS ###
  const similarityThresholdInput = ui.input.float(MCL_INPUTS.THRESHOLD, {
    value: mclOptions.threshold, nullable: false, onValueChanged: (value) => mclOptions.threshold = value ?? mclOptions.threshold,
  });
  similarityThresholdInput.setTooltip('Threshold for similarity between two sequences to create edges');

  const inflationInput = ui.input.float(MCL_INPUTS.INFLATION, {
    value: mclOptions.inflation, nullable: false, onValueChanged: (value) => mclOptions.inflation = value ?? mclOptions.inflation,
  });
  inflationInput.setTooltip('Inflation value for MCL algorithm');

  const maxIterationsInput = ui.input.int(MCL_INPUTS.MAX_ITERATIONS, {
    value: mclOptions.maxIterations, nullable: false, onValueChanged: (value) => mclOptions.maxIterations = value ?? mclOptions.maxIterations,
  });
  maxIterationsInput.setTooltip('Maximum iterations for MCL algorithm');

  // disable input if there is no gpu
  const useWebGPUInput = ui.input.bool(MCL_INPUTS.USE_WEBGPU, {
    value: mclOptions.useWebGPU, onValueChanged: (value) => mclOptions.useWebGPU = value,
  });
  useWebGPUInput.enabled = false;
  mclOptions.webGPUDescriptionPromise.then(() => {
    if (mclOptions.webGPUDescription !== type.webGPUNotSupported) {
      useWebGPUInput.setTooltip(`Use WebGPU for MCL algorithm (${mclOptions.webGPUDescription})`);
      useWebGPUInput.enabled = true;
    } else {
      useWebGPUInput.setTooltip(type.webGPUNotSupported);
      useWebGPUInput.enabled = false;
      useWebGPUInput.value = false;
    }
  });

  const minClusterSizeInput = ui.input.int(MCL_INPUTS.MIN_CLUSTER_SIZE, {
    value: mclOptions.minClusterSize, nullable: false, onValueChanged: (value) => mclOptions.minClusterSize = value ?? mclOptions.minClusterSize,
  });
  minClusterSizeInput.setTooltip('Minimum cluster size for MCL algorithm');

  const mclDistanceFunctionInput= ui.input.choice(MCL_INPUTS.DISTANCE_FUNCTION,
    {value: mclOptions.distanceF, items: [MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH, MmDistanceFunctionsNames.MONOMER_CHEMICAL_DISTANCE,
      MmDistanceFunctionsNames.HAMMING, MmDistanceFunctionsNames.LEVENSHTEIN], nullable: false,
    onValueChanged: (value) => mclOptions.distanceF = value}) as DG.ChoiceInput<MmDistanceFunctionsNames>;
  const mclGapOpenInput = ui.input.float(MCL_INPUTS.GAP_OPEN, {value: mclOptions.gapOpen,
    onValueChanged: (value) => mclOptions.gapOpen = value ?? mclOptions.gapOpen});
  const mclGapExtendInput = ui.input.float(MCL_INPUTS.GAP_EXTEND, {value: mclOptions.gapExtend,
    onValueChanged: (value) => mclOptions.gapExtend = value ?? mclOptions.gapExtend});
  const mclFingerprintTypesInput: DG.ChoiceInput<string> = ui.input.choice(MCL_INPUTS.FINGERPRINT_TYPE, {value: mclOptions.fingerprintType,
    items: ['Morgan', 'RDKit', 'Pattern', 'AtomPair', 'MACCS', 'TopologicalTorsion'], nullable: false,
    onValueChanged: (value) => mclOptions.fingerprintType = value}) as DG.ChoiceInput<string>;


  const mclInputs = [similarityThresholdInput, inflationInput, maxIterationsInput, minClusterSizeInput,
    mclDistanceFunctionInput, mclFingerprintTypesInput, mclGapOpenInput, mclGapExtendInput, useWebGPUInput];

  const mclInputsHost = ui.form(mclInputs);
  mclInputsHost.style.display = 'none';

  const settingsIcon = ui.icons.settings(() => {
    mclInputsHost.style.display = mclInputsHost.style.display === 'none' ? 'flex' : 'none';
    mclInputsHost.classList.remove('ui-form-condensed');
  }, 'Adjust clustering parameters');
  settingsIcon.style.fontSize = '16px';
  generateClustersInput.root.appendChild(settingsIcon);
  // ### END MCL INPUTS ###


  const bitsetChanged = df.filter.onChanged.subscribe(() => activityScalingMethodState());

  const startAnalysisCallback = async (): Promise<boolean> => {
    const sequencesCol = col ?? seqColInput!.value;
    bitsetChanged.unsubscribe();
    if (sequencesCol) {
      const model = await startAnalysis(activityColumnChoice.value!, sequencesCol, clustersColumnChoice.value, df,
        scaledCol, activityScalingMethod.value ?? C.SCALING_METHODS.NONE, {addSequenceSpace: false, addMCL: true,
          useEmbeddingsClusters: generateClustersInput.value ?? false, mclSettings: mclOptions});
      return model !== null;
    }
    return false;
  };

  let bottomHeight = 'auto';
  const inputElements: HTMLElement[] = [ui.divV(inputsList)];
  $(inputElements[0]).find('label').css('width', 'unset');
  if (typeof col !== 'undefined') {
    const startBtn = ui.button('Launch SAR', startAnalysisCallback, '');
    startBtn.style.alignSelf = 'center';
    inputElements.push(startBtn);
    bottomHeight = '215px';
  }

  $(logoHost).empty().append(ui.wait(async () => {
    const viewer = await df.plot.fromType('WebLogo', {sequenceColumnName: col?.name ?? seqColInput!.value!.name});
    viewer.root.style.setProperty('height', '130px');
    return viewer.root;
  }));

  const mainHost = ui.divV([
    logoHost,
    ui.splitH([
      ui.splitV(inputElements),
      histogramHost,
    ], {style: {height: bottomHeight, minWidth: '500px', maxWidth: '600px'}}),
    mclInputsHost,
  ]);
  return {host: mainHost, callback: startAnalysisCallback};
}

type AnalysisOptions = {
  addSequenceSpace?: boolean,
  useEmbeddingsClusters?: boolean,
  addMCL?: boolean,
  mclSettings?: type.MCLSettings,
};

/**
 * Creates dataframe to use in analysis, model instance and adds viewers
 * @param activityColumn - Activity column
 * @param peptidesCol - Peptides column
 * @param clustersColumn - Clusters column or null
 * @param sourceDf - Source dataframe
 * @param scaledCol - Scaled activity column
 * @param scaling - Activity scaling method
 * @param options - Additional options
 * @return - Peptides model instance or null
 */
export async function startAnalysis(activityColumn: DG.Column<number>, peptidesCol: DG.Column<string>,
  clustersColumn: DG.Column | null, sourceDf: DG.DataFrame, scaledCol: DG.Column<number>, scaling: C.SCALING_METHODS,
  options: AnalysisOptions = {}): Promise<PeptidesModel | null> {
  let model: PeptidesModel | null = null;
  if (activityColumn.type !== DG.COLUMN_TYPE.FLOAT && activityColumn.type !== DG.COLUMN_TYPE.INT &&
    activityColumn.type !== DG.COLUMN_TYPE.QNUM
  ) {
    grok.shell.error('The activity column must be of numeric type!');
    return model;
  }
  const progress = DG.TaskBarProgressIndicator.create('Loading SAR...');

  // Prepare new DF
  // const newDf = DG.DataFrame.create(sourceDf.rowCount);
  // newDf.name = 'Peptides analysis';
  // const newDfCols = newDf.columns;
  // newDfCols.add(scaledCol);
  // for (const col of sourceDf.columns) {
  //   if (col.getTag(C.TAGS.ANALYSIS_COL) !== `${true}`) {
  //     if (col.name.toLowerCase() === scaledCol.name.toLowerCase())
  //       col.name = sourceDf.columns.getUnusedName(col.name);
  //     newDfCols.add(col);
  //   }
  // }

  //make sure the data sync is turned off for the dataframe:
  sourceDf.tags.delete && sourceDf.tags.delete('.script');

  const sourceCols = sourceDf.columns;
  const oldActivityCol = sourceDf.col(scaledCol.name);
  if (oldActivityCol)
    oldActivityCol.name = sourceCols.getUnusedName(oldActivityCol.name);
  const scaleColRawData = scaledCol.getRawData();
  sourceDf.columns.addNew(scaledCol.name, scaledCol.type).init((i) => scaleColRawData[i]);
  sourceCols.setOrder([scaledCol.name, peptidesCol.name, ...sourceCols.names().filter((name) => name !== peptidesCol.name && name !== scaledCol.name)]);
  const settings: type.PeptidesSettings = {
    sequenceColumnName: peptidesCol.name, activityColumnName: activityColumn.name, activityScaling: scaling,
    columns: {}, showDendrogram: false, showSequenceSpace: false,
    sequenceSpaceParams: new type.SequenceSpaceParams(!!options.useEmbeddingsClusters && !clustersColumn),
    mclSettings: options.mclSettings ?? new type.MCLSettings(),
  };

  if (clustersColumn) {
    const clusterCol = sourceDf.getCol(clustersColumn.name);
    if (clusterCol.type !== DG.COLUMN_TYPE.STRING)
      sourceCols.replace(clusterCol, clusterCol.convertTo(DG.COLUMN_TYPE.STRING));
  }
  sourceDf.setTag(C.TAGS.SETTINGS, JSON.stringify(settings));

  // const bitset = DG.BitSet.create(sourceDf.rowCount,
  //   (i) => !activityColumn.isNone(i) && !peptidesCol.isNone(i) && sourceDf.filter.get(i));

  // Cloning dataframe with applied filter. If filter is not applied, cloning is
  // needed anyway to allow filtering on the original dataframe
  model = PeptidesModel.getInstance(sourceDf);
  model.init(settings);
  if (clustersColumn) {
    const lstProps: ILogoSummaryTable = {
      clustersColumnName: clustersColumn.name, sequenceColumnName: peptidesCol.name, activityScaling: scaling,
      activityColumnName: activityColumn.name,
    };
    await model.addLogoSummaryTable(lstProps);
  }
  await model.addMonomerPosition();
  await model.addMostPotentResidues();

  // FIXME: enable by default for tests
  if (options.addSequenceSpace ?? false) {
    await model.addSequenceSpace({clusterCol: clustersColumn, clusterEmbeddings: options.useEmbeddingsClusters});
    if (!clustersColumn && (options.useEmbeddingsClusters ?? false)) {
      const clusterCol = model._sequenceSpaceCols
        .find((col) => model!.df.col(col) && model!.df.col(col)?.type === DG.COLUMN_TYPE.STRING);
      if (clusterCol) {
        const lstProps: ILogoSummaryTable = {
          clustersColumnName: clusterCol, sequenceColumnName: peptidesCol.name, activityScaling: scaling,
          activityColumnName: activityColumn.name,
        };
        await model.addLogoSummaryTable(lstProps);
        setTimeout(() => {
          model && (model?.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable)?.render &&
          (model?.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable)?.render();
        }, 100);
      }
    }
  } else if (options.addMCL ?? false) {
    await model.addMCLClusters();
    if (!clustersColumn && (options.useEmbeddingsClusters ?? false)) {
      const mclClusterCol = model._mclCols
        .find(
          (col) => model?.df.col(col) && col.toLowerCase().startsWith('cluster') && !col.toLowerCase().includes('size'),
        );
      if (mclClusterCol) {
        const lstProps: ILogoSummaryTable = {
          clustersColumnName: mclClusterCol, sequenceColumnName: peptidesCol.name, activityScaling: scaling,
          activityColumnName: activityColumn.name,
        };
        await model.addLogoSummaryTable(lstProps);
        setTimeout(() => {
          model && (model?.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable)?.render &&
          (model?.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable)?.render();
        }, 100);
      }
    }
  }


  progress.close();
  return model;
}
