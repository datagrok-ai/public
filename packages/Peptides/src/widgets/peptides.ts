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
export async function analyzePeptidesWidget(df: DG.DataFrame, col: DG.Column): Promise<DG.Widget> {
  if (!col.tags['aligned']?.includes('MSA') && col.tags[DG.TAGS.UNITS].toLowerCase() != 'helm')
    return new DG.Widget(ui.divText('Peptides analysis only works with aligned sequences'));

  let funcs = DG.Func.find({package: 'Bio', name: 'webLogoViewer'});
  if (funcs.length == 0)
    return new DG.Widget(ui.label('Bio package is missing or out of date. Please install the latest version.'));

  funcs = DG.Func.find({package: 'Helm', name: 'getMonomerLib'});
  if (funcs.length == 0)
    return new DG.Widget(ui.label('Helm package is missing or out of date. Please install the latest version.'));

  let tempCol = null;
  let scaledDf: DG.DataFrame;
  let newScaledColName: string;
  let scalingFormula: (x: number) => number;

  for (const column of df.columns.numerical)
    tempCol = column.type === DG.TYPE.FLOAT ? column : null;

  const defaultActivityColumn: DG.Column<number> | null = df.col('activity') || df.col('IC50') || tempCol;
  const histogramHost = ui.div([], {id: 'pep-hist-host'});

  const indexes: number[] = [];
  const f = df.filter;
  df.onFilterChanged.subscribe(() => {
    for (let i = 0; i < f.length; ++i) {
      if (f.get(i))
        indexes.push(i);
    }
  });
  const activityScalingMethod = ui.choiceInput(
    'Scaling', 'none', ['none', 'lg', '-lg'],
    async (currentMethod: string): Promise<void> => {
      [scaledDf, scalingFormula, newScaledColName] =
        scaleActivity(currentMethod, activityColumnChoice.value!, indexes.length !== 0 ? indexes : undefined);

      const hist = scaledDf.plot.histogram({
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

  const activityScalingMethodState = (_: any): void => {
    activityScalingMethod.enabled = (activityColumnChoice.value ?? false) &&
      DG.Stats.fromColumn(activityColumnChoice.value!, df.filter).min > 0;
    activityScalingMethod.fireChanged();
  };
  const activityColumnChoice = ui.columnInput('Activity', df, defaultActivityColumn, activityScalingMethodState);
  const clustersColumnChoice = ui.columnInput('Clusters', df, null);
  activityColumnChoice.fireChanged();
  activityScalingMethod.fireChanged();

  const inputsList = [activityColumnChoice, activityScalingMethod, clustersColumnChoice];

  const startBtn = ui.button('Launch SAR', async () => {
    await startAnalysis(activityColumnChoice.value, col, clustersColumnChoice.value, df, scalingFormula,
      newScaledColName, activityScalingMethod.value ?? 'none', indexes);
  });
  startBtn.style.alignSelf = 'center';

  const viewer = await df.plot.fromType('WebLogo') as bio.WebLogo;
  viewer.root.style.setProperty('height', '130px');
  const logoHost = ui.div();
  $(logoHost).empty().append(viewer.root);

  return new DG.Widget(
    ui.divV([
      logoHost,
      ui.splitH([
        ui.splitV([ui.inputs(inputsList), startBtn]),
        histogramHost,
      ], {style: {height: '215px'}}),
    ]),
  );
}

export async function startAnalysis(activityColumn: DG.Column<number> | null, peptidesCol: DG.Column<string>,
  clustersColumn: DG.Column | null, currentDf: DG.DataFrame, scaleNum: (x: number) => number, newScaledColName: string,
  scaling: string, indexes: number[]): Promise<PeptidesModel | null> {
  const progress = DG.TaskBarProgressIndicator.create('Loading SAR...');
  let model = null;
  if (activityColumn?.type === DG.TYPE.FLOAT) {
    const f = currentDf.filter;
    //prepare new DF
    const newDf = DG.DataFrame.create(f.trueCount);
    const getIndex = indexes.length !== 0 ? (i: number): number => indexes[i] : (i: number): number => i;
    let activityCol: DG.Column<number> | null = null;
    for (const col of currentDf.columns.toList()) {
      let virtualCol: DG.Column<any>;
      if (col === activityColumn) {
        virtualCol = newDf.columns.addNewVirtual(
          C.COLUMNS_NAMES.ACTIVITY, (i) => activityColumn.get(getIndex(i)!), DG.TYPE.FLOAT);
        activityCol = virtualCol;
      } else if (col === peptidesCol) {
        virtualCol = newDf.columns.addNewVirtual(
          C.COLUMNS_NAMES.MACROMOLECULE, (i) => peptidesCol.get(getIndex(i)!), DG.TYPE.STRING);
      } else
        virtualCol = newDf.columns.addNewVirtual(col.name, (i) => col.get(getIndex(i)!), col.type as DG.TYPE);
      virtualCol.setTag(C.TAGS.VISIBLE, '0');
    }
    activityCol!.semType = C.SEM_TYPES.ACTIVITY;
    const activityScaledCol = newDf.columns.addNewVirtual(C.COLUMNS_NAMES.ACTIVITY_SCALED, (i) => {
      const val = activityCol!.get(getIndex(i)!);
      return val ? scaleNum(val) : val;
    }, DG.TYPE.FLOAT);
    activityScaledCol.semType = C.SEM_TYPES.ACTIVITY_SCALED;
    newDf.name = 'Peptides analysis';
    newDf.tags[C.COLUMNS_NAMES.ACTIVITY_SCALED] = newScaledColName;
    if (clustersColumn) {
      newDf.getCol(clustersColumn.name).name = C.COLUMNS_NAMES.CLUSTERS;
      newDf.tags[C.TAGS.CLUSTERS] = C.COLUMNS_NAMES.CLUSTERS;
    }
    // newDf.tags[C.PEPTIDES_ANALYSIS] = 'true';
    newDf.tags['scaling'] = scaling;

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
    grok.shell.error('The activity column must be of floating point number type!');
  progress.close();
  return model;
}
