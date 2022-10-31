import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import '../styles.css';
import * as C from '../utils/constants';
import {PeptidesModel} from '../model';
import $ from 'cash-dom';
import {scaleActivity} from '../utils/misc';

/** Peptide analysis widget.
 *
 * @param {DG.DataFrame} df Working table
 * @param {DG.Column} col Aligned sequence column
 * @return {Promise<DG.Widget>} Widget containing peptide analysis */
export async function analyzePeptidesWidget(df: DG.DataFrame, col?: DG.Column<string>): Promise<DG.Widget> {
  let seqColInput: DG.InputBase | null = null;
  if (typeof col === 'undefined') {
    const sequenceColumns = df.columns.toList().filter((dfCol) => dfCol.semType === DG.SEMTYPE.MACROMOLECULE);
    let potentialCol = DG.Utils.firstOrNull(sequenceColumns);
    if (potentialCol === null)
      throw new Error('Peptides Error: table doesn\'t contain sequence columns');
    seqColInput = ui.columnInput('Sequence', df, potentialCol, () => {
      const seqCol = seqColInput!.value;
      if (!seqCol.tags['aligned']?.includes('MSA') && seqCol.tags[DG.TAGS.UNITS].toLowerCase() !== 'helm')
        grok.shell.warning('Peptides analysis only works with aligned sequences');
    });
  } else if (!col.tags['aligned']?.includes('MSA') && col.tags[DG.TAGS.UNITS].toLowerCase() !== 'helm')
    return new DG.Widget(ui.divText('Peptides analysis only works with aligned sequences'));

  let funcs = DG.Func.find({package: 'Bio', name: 'webLogoViewer'});
  if (funcs.length == 0)
    return new DG.Widget(ui.label('Bio package is missing or out of date. Please install the latest version.'));

  funcs = DG.Func.find({package: 'Helm', name: 'getMonomerLib'});
  if (funcs.length == 0)
    return new DG.Widget(ui.label('Helm package is missing or out of date. Please install the latest version.'));

  let scaledCol: DG.Column<number>;

  const defaultActivityColumn: DG.Column<number> | null =
    df.col('activity') || df.col('IC50') || DG.Utils.firstOrNull(df.columns.numerical); ;
  const histogramHost = ui.div([], {id: 'pep-hist-host'});

  const activityScalingMethod = ui.choiceInput(
    'Scaling', 'none', ['none', 'lg', '-lg'],
    async (currentMethod: string): Promise<void> => {
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
    });
  activityScalingMethod.setTooltip('Function to apply for each value in activity column');

  const activityScalingMethodState = (): void => {
    activityScalingMethod.enabled = (activityColumnChoice.value ?? false) &&
      DG.Stats.fromColumn(activityColumnChoice.value!).min > 0;
    activityScalingMethod.fireChanged();
  };
  const activityColumnChoice = ui.columnInput('Activity', df, defaultActivityColumn, activityScalingMethodState);
  const clustersColumnChoice = ui.columnInput('Clusters', df, null);
  activityColumnChoice.fireChanged();
  activityScalingMethod.fireChanged();

  const inputsList = [activityColumnChoice, activityScalingMethod, clustersColumnChoice];
  if (seqColInput !== null)
    inputsList.splice(0, 0, seqColInput);

  const bitsetChanged = df.filter.onChanged.subscribe(() => {
    activityScalingMethodState();
  });

  const startBtn = ui.button('Launch SAR', async () => {
    const sequencesCol = col ?? seqColInput!.value;
    if (sequencesCol)
    await startAnalysis(activityColumnChoice.value!, sequencesCol, clustersColumnChoice.value, df, scaledCol,
      activityScalingMethod.value ?? 'none');
    bitsetChanged.unsubscribe();
  });
  startBtn.style.alignSelf = 'center';

  const viewer = await df.plot.fromType('WebLogo') as bio.WebLogoViewer;
  viewer.root.style.setProperty('height', '130px');
  const logoHost = ui.div();
  $(logoHost).empty().append(viewer.root);

  const mainHost = ui.divV([
    logoHost,
    ui.splitH([
      ui.splitV([ui.inputs(inputsList), startBtn]),
      histogramHost,
    ], {style: {height: '215px'}}),
  ]);
  mainHost.style.maxWidth = '400px';
  return new DG.Widget(mainHost);
}

export async function startAnalysis(activityColumn: DG.Column<number>, peptidesCol: DG.Column<string>,
  clustersColumn: DG.Column | null, currentDf: DG.DataFrame, scaledCol: DG.Column<number>, scaling: string,
): Promise<PeptidesModel | null> {
  const progress = DG.TaskBarProgressIndicator.create('Loading SAR...');
  let model = null;
  if (activityColumn.type === DG.TYPE.FLOAT || activityColumn.type === DG.TYPE.INT) {
    //prepare new DF
    const newDf = DG.DataFrame.create(currentDf.rowCount);
    const newDfCols = newDf.columns;
    for (const col of currentDf.columns.toList()) {
      const currentCol = newDfCols.add(col);
      if (col === activityColumn)
        currentCol.name = C.COLUMNS_NAMES.ACTIVITY;
      else if (col === peptidesCol)
        currentCol.name = C.COLUMNS_NAMES.MACROMOLECULE;
      col.setTag(C.TAGS.VISIBLE, '0');
    }
    activityColumn.semType = C.SEM_TYPES.ACTIVITY;
    newDfCols.add(scaledCol);
    newDf.name = 'Peptides analysis';
    if (clustersColumn) {
      newDf.getCol(clustersColumn.name).name = C.COLUMNS_NAMES.CLUSTERS;
      newDf.tags[C.TAGS.CLUSTERS] = C.COLUMNS_NAMES.CLUSTERS;
    }
    // newDf.tags['scaling'] = scaling;
    newDf.setTag('settings', JSON.stringify({scaling: scaling}));

    let monomerType = 'HELM_AA';
    if (peptidesCol.getTag(DG.TAGS.UNITS).toLowerCase() == 'helm') {
      const sampleSeq = peptidesCol.get(0)!;
      monomerType = sampleSeq.startsWith('PEPTIDE') ? 'HELM_AA' : 'HELM_BASE';
    } else {
      const alphabet = peptidesCol.tags[C.TAGS.ALPHABET];
      monomerType = alphabet == 'DNA' || alphabet == 'RNA' ? 'HELM_BASE' : 'HELM_AA';
    }

    newDf.setTag('monomerType', monomerType);
    model = await PeptidesModel.getInstance(newDf);
  } else
    grok.shell.error('The activity column must be of numeric type!');
  progress.close();
  return model;
}
