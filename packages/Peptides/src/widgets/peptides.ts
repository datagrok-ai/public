import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../styles.css';
import * as C from '../utils/constants';
import * as type from '../utils/types';
import {PeptidesModel} from '../model';
import $ from 'cash-dom';
import {scaleActivity} from '../utils/misc';
import * as bio from '@datagrok-libraries/bio';

/** Peptide analysis widget.
 *
 * @param {DG.DataFrame} df Working table
 * @param {DG.Column} col Aligned sequence column
 * @return {Promise<DG.Widget>} Widget containing peptide analysis */
export function analyzePeptidesUI(df: DG.DataFrame, col?: DG.Column<string>):
{host: HTMLElement, callback: () => Promise<boolean>} {
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
      if (!(seqCol.getTag(DG.TAGS.SEMTYPE) == DG.SEMTYPE.MACROMOLECULE)) {
        grok.shell.warning('Peptides analysis only works with macromolecules');
        seqColInput!.value = potentialCol;
      }
      $(logoHost).empty().append(ui.wait(async () => {
        const viewer = await df.plot.fromType('WebLogo', {sequenceColumnName: seqCol.name});
        viewer.root.style.setProperty('height', '130px');
        return viewer.root;
      }));
    });
  } else if (!(col.getTag(bio.TAGS.aligned) == bio.ALIGNMENT.SEQ_MSA) && col.getTag(DG.TAGS.UNITS) !== bio.NOTATION.HELM) {
    return {
      host: ui.label('Peptides analysis only works with aligned sequences'),
      callback: async (): Promise<boolean> => false,
    };
  }

  let funcs = DG.Func.find({package: 'Bio', name: 'webLogoViewer'});
  if (funcs.length == 0) {
    return {
      host: ui.label('Bio package is missing or out of date. Please install the latest version.'),
      callback: async (): Promise<boolean> => false,
    };
  }

  funcs = DG.Func.find({package: 'Bio', name: 'getBioLib'});
  if (funcs.length == 0) {
    return {
      host: ui.label('Bio package is missing or out of date. Please install the latest version.'),
      callback: async (): Promise<boolean> => false,
    };
  }

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

  const startAnalysisCallback = async (): Promise<boolean> => {
    const sequencesCol = col ?? seqColInput!.value;
    if (sequencesCol) {
      const model = await startAnalysis(activityColumnChoice.value!, sequencesCol, clustersColumnChoice.value, df,
        scaledCol, activityScalingMethod.stringValue as type.ScalingMethods ?? 'none');
      if (model != null) {
        bitsetChanged.unsubscribe();
        return true;
      }
    }
    return false;
  };

  const inputElements: HTMLElement[] = [ui.inputs(inputsList)];
  $(inputElements[0]).find('label').css('width', 'unset');
  if (typeof col !== 'undefined') {
    const startBtn = ui.button('Launch SAR', startAnalysisCallback);
    startBtn.style.alignSelf = 'center';
    inputElements.push(startBtn);
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
    ], {style: {height: '215px'}}),
  ]);
  mainHost.style.maxWidth = '400px';
  return {host: mainHost, callback: startAnalysisCallback};
}

export async function startAnalysis(activityColumn: DG.Column<number>, peptidesCol: DG.Column<string>,
  clustersColumn: DG.Column | null, currentDf: DG.DataFrame, scaledCol: DG.Column<number>, scaling: type.ScalingMethods,
): Promise<PeptidesModel | null> {
  const progress = DG.TaskBarProgressIndicator.create('Loading SAR...');
  let model = null;
  if (activityColumn.type === DG.TYPE.FLOAT || activityColumn.type === DG.TYPE.INT) {
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
    };
    if (clustersColumn) {
      const clusterCol = newDf.getCol(clustersColumn.name);
      newDf.columns.replace(clusterCol, clusterCol.convertTo(DG.COLUMN_TYPE.STRING));
      settings.clustersColumnName = clustersColumn.name;
    }
    newDf.setTag(C.TAGS.SETTINGS, JSON.stringify(settings));

    let monomerType = 'HELM_AA';
    if (peptidesCol.getTag(DG.TAGS.UNITS) == bio.NOTATION.HELM) {
      const sampleSeq = peptidesCol.get(0)!;
      monomerType = sampleSeq.startsWith('PEPTIDE') ? 'HELM_AA' : 'HELM_BASE';
    } else {
      const alphabet = peptidesCol.tags[C.TAGS.ALPHABET];
      monomerType = alphabet == 'DNA' || alphabet == 'RNA' ? 'HELM_BASE' : 'HELM_AA';
    }

    newDf.setTag('monomerType', monomerType);
    newDf.setTag('newAnalysis', '1');
    model = PeptidesModel.getInstance(newDf);
    await model.addViewers();
  } else
    grok.shell.error('The activity column must be of numeric type!');
  progress.close();
  return model;
}
