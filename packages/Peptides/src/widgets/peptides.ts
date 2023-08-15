import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as uuid from 'uuid';

import '../styles.css';
import * as C from '../utils/constants';
import * as type from '../utils/types';
import {PeptidesModel} from '../model';
import $ from 'cash-dom';
import {scaleActivity} from '../utils/misc';
import {ALIGNMENT, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';

/** Peptide analysis widget.
 * @param {DG.DataFrame} df Working table
 * @param {DG.Column} col Aligned sequence column
 * @return {Promise<DG.Widget>} Widget containing peptide analysis */
export function analyzePeptidesUI(df: DG.DataFrame, col?: DG.Column<string>):
  { host: HTMLElement, callback: () => Promise<boolean> } {
  const logoHost = ui.div();
  // logoHost.style.alignContent = 'center';
  let seqColInput: DG.InputBase | null = null;
  if (typeof col === 'undefined') {
    const sequenceColumns = df.columns.toList().filter((dfCol) => dfCol.semType === DG.SEMTYPE.MACROMOLECULE);
    const potentialCol = DG.Utils.firstOrNull(sequenceColumns);
    if (potentialCol === null)
      throw new Error('Peptides Error: table doesn\'t contain sequence columns');

    seqColInput = ui.columnInput('Sequence', df, potentialCol, () => {
      const seqCol = seqColInput!.value;
      if (!(seqCol.getTag(DG.TAGS.SEMTYPE) === DG.SEMTYPE.MACROMOLECULE)) {
        grok.shell.warning('Peptides analysis only works with macromolecules');
        seqColInput!.value = potentialCol;
      }
      $(logoHost).empty().append(ui.wait(async () => {
        const viewer = await df.plot.fromType('WebLogo', {sequenceColumnName: seqCol.name});
        viewer.root.style.setProperty('height', '130px');
        return viewer.root;
      }));
      //TODO: add when new version of datagrok-api is available
    }, {filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE});
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

  let scaledCol: DG.Column<number>;

  const defaultActivityColumn: DG.Column<number> | null =
    df.col('activity') || df.col('IC50') || DG.Utils.firstOrNull(df.columns.numerical);
  ;
  const histogramHost = ui.div([], {id: 'pep-hist-host'});

  const activityScalingMethod = ui.choiceInput(
    'Scaling', C.SCALING_METHODS.NONE, Object.values(C.SCALING_METHODS),
    async (currentMethod: C.SCALING_METHODS): Promise<void> => {
      scaledCol = scaleActivity(activityColumnChoice.value!, currentMethod);

      const hist = DG.DataFrame.fromColumns([scaledCol]).plot.histogram({
        filteringEnabled: false,
        valueColumnName: C.COLUMNS_NAMES.ACTIVITY_SCALED,
        legendVisibility: 'Never',
        showXAxis: true,
        showColumnSelector: false,
        showRangeSlider: false,
        showBinSelector: false,
      });
      histogramHost.lastChild?.remove();
      histogramHost.appendChild(hist.root);
    }) as DG.InputBase<C.SCALING_METHODS | null>;
  activityScalingMethod.setTooltip('Function to apply for each value in activity column');

  const activityScalingMethodState = (): void => {
    activityScalingMethod.enabled = (activityColumnChoice.value ?? false) &&
      DG.Stats.fromColumn(activityColumnChoice.value!).min > 0;
    activityScalingMethod.fireChanged();
  };
  //TODO: add when new version of datagrok-api is available
  const activityColumnChoice = ui.columnInput('Activity', df, defaultActivityColumn, activityScalingMethodState,
    {filter: (col: DG.Column) => col.type === DG.TYPE.INT || col.type === DG.TYPE.FLOAT});
  const clustersColumnChoice = ui.columnInput('Clusters', df, null);
  clustersColumnChoice.nullable = true;
  activityColumnChoice.fireChanged();
  activityScalingMethod.fireChanged();

  const targetColumnChoice = ui.columnInput('Target', df, null, null,
    {filter: (col: DG.Column) => col.type === DG.TYPE.STRING});
  targetColumnChoice.nullable = true;

  const inputsList = [activityColumnChoice, activityScalingMethod, clustersColumnChoice, targetColumnChoice];
  if (seqColInput !== null)
    inputsList.splice(0, 0, seqColInput);

  const bitsetChanged = df.filter.onChanged.subscribe(() => activityScalingMethodState());

  const startAnalysisCallback = async (): Promise<boolean> => {
    const sequencesCol = col ?? seqColInput!.value;
    bitsetChanged.unsubscribe();
    if (sequencesCol) {
      const model = await startAnalysis(activityColumnChoice.value!, sequencesCol, clustersColumnChoice.value, df,
        scaledCol, activityScalingMethod.value ?? C.SCALING_METHODS.NONE, targetColumnChoice.value);
      return model !== null;
    }
    return false;
  };

  let bottomHeight = 'auto';
  const inputElements: HTMLElement[] = [ui.inputs(inputsList)];
  $(inputElements[0]).find('label').css('width', 'unset');
  if (typeof col !== 'undefined') {
    const startBtn = ui.button('Launch SAR', startAnalysisCallback);
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
    ], {style: {height: bottomHeight}}),
  ]);
  mainHost.style.maxWidth = '400px';
  return {host: mainHost, callback: startAnalysisCallback};
}

export async function startAnalysis(activityColumn: DG.Column<number>, peptidesCol: DG.Column<string>,
  clustersColumn: DG.Column | null, currentDf: DG.DataFrame, scaledCol: DG.Column<number>, scaling: C.SCALING_METHODS,
  targetColumn: DG.Column<string> | null = null): Promise<PeptidesModel | null> {
  const progress = DG.TaskBarProgressIndicator.create('Loading SAR...');
  let model = null;
  if (activityColumn.type === DG.COLUMN_TYPE.FLOAT || activityColumn.type === DG.COLUMN_TYPE.INT) {
    //prepare new DF
    const newDf = DG.DataFrame.create(currentDf.rowCount);
    const newDfCols = newDf.columns;
    newDfCols.add(scaledCol);
    for (const col of currentDf.columns)
      newDfCols.add(col);

    newDf.name = 'Peptides analysis';
    const settings: type.PeptidesSettings = {
      sequenceColumnName: peptidesCol.name,
      activityColumnName: activityColumn.name,
      scaling: scaling,
      isBidirectional: false,
      columns: {},
      maxMutations: 1,
      minActivityDelta: 0,
      showDendrogram: false,
    };
    if (targetColumn !== null)
      settings.targetColumnName = targetColumn.name;

    if (clustersColumn) {
      const clusterCol = newDf.getCol(clustersColumn.name);
      if (clusterCol.type !== DG.COLUMN_TYPE.STRING)
        newDfCols.replace(clusterCol, clusterCol.convertTo(DG.COLUMN_TYPE.STRING));
      settings.clustersColumnName = clustersColumn.name;
    }
    newDf.setTag(C.TAGS.SETTINGS, JSON.stringify(settings));

    let monomerType = 'HELM_AA';
    if (peptidesCol.getTag(DG.TAGS.UNITS) === NOTATION.HELM) {
      const sampleSeq = peptidesCol.get(0)!;
      monomerType = sampleSeq.startsWith('PEPTIDE') ? 'HELM_AA' : 'HELM_BASE';
    } else {
      const alphabet = peptidesCol.tags[C.TAGS.ALPHABET];
      monomerType = alphabet === 'DNA' || alphabet === 'RNA' ? 'HELM_BASE' : 'HELM_AA';
    }
    const dfUuid = uuid.v4();
    newDf.setTag(C.TAGS.UUID, dfUuid);
    newDf.setTag('monomerType', monomerType);

    // Cloning dataframe with applied filter. If filter is not applied, cloning is
    // needed anyway to allow filtering on the original dataframe
    model = PeptidesModel.getInstance(newDf.clone(currentDf.filter));
    if (clustersColumn) await model.addLogoSummaryTable();
    await model.addMonomerPosition();
    await model.addMostPotentResidues();
  } else
    grok.shell.error('The activity column must be of numeric type!');
  progress.close();
  return model;
}
