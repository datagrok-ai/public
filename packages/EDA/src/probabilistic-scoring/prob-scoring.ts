/* eslint-disable max-len */
// Probabilistic scoring (pMPO) features
// Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC4716604/

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

import '../../css/pmpo.css';

import {getDesiredTables, getDescriptorStatistics, gaussDesirabilityFunc, sigmoidS,
  getBoolPredictionColumn, getPmpoEvaluation} from './stat-tools';
import {MIN_SAMPLES_COUNT, PMPO_NON_APPLICABLE, DescriptorStatistics, P_VAL_TRES_MIN, DESCR_TITLE,
  R2_MIN, Q_CUTOFF_MIN, PmpoParams, SCORES_TITLE, DESCR_TABLE_TITLE, PMPO_COMPUTE_FAILED, SELECTED_TITLE,
  P_VAL, DESIRABILITY_COL_NAME, STAT_GRID_HEIGHT, DESIRABILITY_COLUMN_WIDTH, WEIGHT_TITLE,
  P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT, USE_SIGMOID_DEFAULT, ROC_TRESHOLDS,
  FPR_TITLE, TPR_TITLE, COLORS, THRESHOLD, AUTO_TUNE_MAX_APPLICABLE_ROWS, AUTO_TUNE_WARNING_MIN_ROWS,
  DEFAULT_OPTIMIZATION_SETTINGS, P_VAL_TRES_MAX, R2_MAX, Q_CUTOFF_MAX, OptimalPoint, LOW_PARAMS_BOUNDS,
  HIGH_PARAMS_BOUNDS,
  FORMAT} from './pmpo-defs';
import {addSelectedDescriptorsCol, getDescriptorStatisticsTable, getFilteredByPvalue, getFilteredByCorrelations,
  getModelParams, getDescrTooltip, saveModel, getScoreTooltip, getDesirabilityProfileJson, getCorrelationTriples,
  addCorrelationColumns, setPvalColumnColorCoding, setCorrColumnColorCoding, PmpoError} from './pmpo-utils';
import {getOutputPalette} from '../pareto-optimization/utils';
import {OPT_TYPE} from '../pareto-optimization/defs';
import {optimizeNM} from './nelder-mead';

export type PmpoTrainingResult = {
  params: Map<string, PmpoParams>,
  descrStatsTable: DG.DataFrame,
  selectedByPvalue: string[],
  selectedByCorr: string[],
};

/** Type for pMPO training controls */
export type Controls = {form: HTMLElement, saveBtn: HTMLButtonElement};

/** Type for pMPO elements */
export type PmpoAppItems = {
  statsGrid: DG.Viewer;
  rocCurve: DG.Viewer;
  confusionMatrix: DG.Viewer;
  controls: Controls;
};

/** Class implementing probabilistic MPO (pMPO) model training and prediction */
export class Pmpo {
  /** Checks if pMPO model can be applied to the given descriptors and desirability column */
  static isApplicable(descriptors: DG.ColumnList, desirability: DG.Column, pValThresh: number,
    r2Tresh: number, qCutoff: number, toShowWarning: boolean = false): boolean {
    const rows = desirability.length;

    const showWarning = (msg: string) => {
      if (toShowWarning)
        grok.shell.warning(PMPO_NON_APPLICABLE + msg);
    };

    // Check p-value threshold
    if (pValThresh < P_VAL_TRES_MIN) {
      showWarning(`: too small p-value threshold - ${pValThresh}, minimum - ${P_VAL_TRES_MIN}`);
      return false;
    }

    // Check R2 threshold
    if (r2Tresh < R2_MIN) {
      showWarning(`: too small R² threshold - ${r2Tresh}, minimum - ${R2_MIN}`);
      return false;
    }

    // Check q-cutoff
    if (qCutoff < Q_CUTOFF_MIN) {
      showWarning(`: too small q-cutoff - ${qCutoff}, minimum - ${Q_CUTOFF_MIN}`);
      return false;
    }

    // Check samples count
    if (rows < MIN_SAMPLES_COUNT) {
      showWarning(`: not enough of samples - ${rows}, minimum - ${MIN_SAMPLES_COUNT}`);
      return false;
    }

    // Check desirability
    if (desirability.type !== DG.COLUMN_TYPE.BOOL) {
      showWarning(`: "${desirability.name}" must be boolean column.`);
      return false;
    }

    if (desirability.stats.stdev === 0) { // TRUE & FALSE
      showWarning(`: "${desirability.name}" has a single category.`);
      return false;
    }

    // Check descriptors
    let nonConstantCols = 0;

    for (const col of descriptors) {
      if (!col.isNumerical) {
        showWarning(`: "${col.name}" is not numerical.`);
        return false;
      }

      if (col.stats.missingValueCount > 0) {
        showWarning(`: "${col.name}" contains missing values.`);
        return false;
      }

      if (col.stats.stdev > 0)
        ++nonConstantCols;
    }

    if (nonConstantCols < 1) {
      showWarning(`: not enough of non-constant descriptors.`);
      return false;
    }

    return true;
  } // isApplicable

  /** Validates the input data frame for pMPO applicability */
  static isTableValid(df: DG.DataFrame, toShowMsg: boolean = true): boolean {
    // Check row count
    if (df.rowCount < 2) {
      if (toShowMsg)
        grok.shell.warning(PMPO_NON_APPLICABLE + `. Not enough of samples: ${df.rowCount}, minimum: 2.`);
      return false;
    }

    let boolColsCount = 0;
    let validNumericColsCount = 0;

    // Check numeric columns and boolean columns
    for (const col of df.columns) {
      if (col.isNumerical) {
        if ((col.stats.missingValueCount < 1) && (col.stats.stdev > 0))
          ++validNumericColsCount;
      } else if (col.type == DG.COLUMN_TYPE.BOOL)
        ++boolColsCount;
    }

    // Check boolean columns count
    if (boolColsCount < 1) {
      if (toShowMsg)
        grok.shell.warning(PMPO_NON_APPLICABLE + ': no boolean columns.');
      return false;
    }

    // Check valid numeric columns count
    if (validNumericColsCount < 1) {
      if (toShowMsg)
        grok.shell.warning(PMPO_NON_APPLICABLE + ': no numeric columns without missing values and non-zero variance.');
      return false;
    }

    return true;
  } // isTableValid

  /** Fits the pMPO model to the given data and returns training results */
  static fit(df: DG.DataFrame, descriptors: DG.ColumnList, desirability: DG.Column,
    pValTresh: number, r2Tresh: number, qCutoff: number, toCheckApplicability: boolean = true): PmpoTrainingResult {
    if (toCheckApplicability) {
      if (!Pmpo.isApplicable(descriptors, desirability, pValTresh, r2Tresh, qCutoff))
        throw new Error('Failed to train pMPO model: the method is not applicable to the inputs');
    }

    const descriptorNames = descriptors.names();
    const {desired, nonDesired} = getDesiredTables(df, desirability);

    // Compute descriptors' statistics
    const descrStats = new Map<string, DescriptorStatistics>();
    descriptorNames.forEach((name) => {
      descrStats.set(name, getDescriptorStatistics(desired.col(name)!, nonDesired.col(name)!));
    });
    const descrStatsTable = getDescriptorStatisticsTable(descrStats);

    // Set p-value column color coding
    setPvalColumnColorCoding(descrStatsTable, pValTresh);

    // Filter by p-value
    const selectedByPvalue = getFilteredByPvalue(descrStatsTable, pValTresh);

    if (selectedByPvalue.length < 1)
      throw new PmpoError('Cannot train pMPO model: all descriptors have high p-values (not significant).');

    // Compute correlation triples
    const correlationTriples = getCorrelationTriples(descriptors, selectedByPvalue);

    // Filter by correlations
    const selectedByCorr = getFilteredByCorrelations(descriptors, selectedByPvalue, descrStats, r2Tresh, correlationTriples);

    // Add the Selected column
    addSelectedDescriptorsCol(descrStatsTable, selectedByCorr);

    // Add correlation columns
    addCorrelationColumns(descrStatsTable, descriptorNames, correlationTriples, selectedByCorr);

    // Set correlation columns color coding
    setCorrColumnColorCoding(descrStatsTable, descriptorNames, r2Tresh);

    // Compute pMPO parameters - training
    const params = getModelParams(desired, nonDesired, selectedByCorr, qCutoff);

    return {
      params: params,
      descrStatsTable: descrStatsTable,
      selectedByPvalue: selectedByPvalue,
      selectedByCorr: selectedByCorr,
    };
  } // fitModelParams

  /** Predicts pMPO scores for the given data frame using provided pMPO parameters */
  static predict(df: DG.DataFrame, params: Map<string, PmpoParams>, useSigmoid: boolean, predictionName: string): DG.Column {
    const count = df.rowCount;
    const scores = new Float64Array(count).fill(0);

    // Compute pMPO scores (see https://pmc.ncbi.nlm.nih.gov/articles/PMC4716604/
    params.forEach((param, name) => {
      const col = df.col(name);
      const b = param.b;
      const c = param.c;
      const x0 = param.cutoff;
      let weight = param.weight;
      const avg = param.desAvg;
      const std = param.desStd;
      const frac = 1.0 / (2 * std**2);

      if (col == null)
        throw new Error(`Failed to apply pMPO: inconsistent data, no column "${name}" in the table "${df.name}"`);

      const vals = col.getRawData();

      if (useSigmoid) {
        if (c > 0) {
          for (let i = 0; i < count; ++i)
            scores[i] += weight * Math.exp(-((vals[i] - avg)**2) * frac) / (1.0 + b * (c ** (-(vals[i] - x0))));
        } else {
          weight = weight / (1.0 + b);

          for (let i = 0; i < count; ++i)
            scores[i] += weight * Math.exp(-((vals[i] - avg)**2) * frac);
        }
      } else {
        for (let i = 0; i < count; ++i)
          scores[i] += weight * Math.exp(-((vals[i] - avg)**2) * frac);
      }
    });

    return DG.Column.fromFloat64Array(predictionName, scores);
  } // predict

  private params: Map<string, PmpoParams> | null = null;

  private table: DG.DataFrame;
  private view: DG.TableView;
  private boolCols: DG.Column[];
  private numericCols: DG.Column[];

  private initTable = DG.DataFrame.create();

  private statGrid = DG.Viewer.grid(this.initTable, {showTitle: true, title: DESCR_TABLE_TITLE});

  private predictionName = SCORES_TITLE;
  private boolPredictionName = '';

  private desirabilityProfileRoots = new Map<string, HTMLElement>();

  private rocCurve = DG.Viewer.scatterPlot(this.initTable, {
    showTitle: true,
    showSizeSelector: false,
    showColorSelector: false,
  });

  private confusionMatrix = DG.Viewer.fromType('Confusion matrix', this.initTable, {
    xColumnName: 'control',
    yColumnName: 'control',
    showTitle: true,
    title: 'Confusion Matrix',
    descriptionPosition: 'Bottom',
    description: 'Confusion matrix for the predicted vs actual desirability labels.',
    descriptionVisibilityMode: 'Always',
  });

  constructor(df: DG.DataFrame, view?: DG.TableView) {
    this.table = df;
    this.view = view ?? (grok.shell.tableView(df.name) ?? grok.shell.addTableView(df));
    this.boolCols = this.getBoolCols();
    this.numericCols = this.getValidNumericCols();
    this.predictionName = df.columns.getUnusedName(SCORES_TITLE);
  };

  /** Sets the ribbon panels in the table view (removes the first panel) */
  private setRibbons(): void {
    const ribPanel = this.view.getRibbonPanels();

    if (ribPanel.length < 1)
      return;

    this.view.setRibbonPanels(ribPanel.slice(1));
  }

  /** Updates the statistics grid viewer with the given statistics table and selected descriptors */
  private updateStatisticsGrid(table: DG.DataFrame, descriptorNames: string[], selectedByPvalue: string[], selectedByCorr: string[]): void {
    const grid = this.statGrid;
    grid.dataFrame = table;
    grid.setOptions({
      showTitle: true,
      title: table.name,
    });

    grid.sort([SELECTED_TITLE], [false]);
    grid.col(P_VAL)!.format = 'scientific';

    // set color coding
    const descrCol = grid.col(DESCR_TITLE)!;
    descrCol.isTextColorCoded = true;

    const pValCol = grid.col(P_VAL)!;
    pValCol.isTextColorCoded = true;

    descriptorNames.forEach((name) => {
      const col = grid.col(name);
      if (col == null)
        return;

      col.isTextColorCoded = true;
      col.format = '0.000';
    });

    // set tooltips
    grid.onCellTooltip((cell, x, y) =>{
      if (cell.isColHeader) {
        const cellCol = cell.tableColumn;

        if (cellCol == null)
          return false;

        const colName = cellCol.name;

        switch (colName) {
        case DESCR_TITLE:
          ui.tooltip.show(getDescrTooltip(
            DESCR_TITLE,
            'Use of descriptors in model construction:',
            'selected',
            'excluded',
          ), x, y);
          return true;

        case DESIRABILITY_COL_NAME:
          ui.tooltip.show(ui.divV([
            ui.h2(DESIRABILITY_COL_NAME),
            ui.divText('Desirability profile charts for each descriptor. Only profiles for selected descriptors are shown.'),
          ]), x, y);
          return true;

        case WEIGHT_TITLE:
          ui.tooltip.show(ui.divV([
            ui.h2(WEIGHT_TITLE),
            ui.divText('Weights of selected descriptors.'),
          ]), x, y);
          return true;

        case P_VAL:
          ui.tooltip.show(getDescrTooltip(
            P_VAL,
            'Filtering descriptors by p-value:',
            'selected',
            'excluded',
          ), x, y);
          return true;

        default:
          if (descriptorNames.includes(colName)) {
            ui.tooltip.show(
              getDescrTooltip(
                colName,
                `Correlation of ${colName} with other descriptors, measured by R²:`,
                'weakly correlated',
                'highly correlated',
              ), x, y);

            return true;
          }

          return false;
        }
      } else {
        if (cell.isTableCell) {
          const cellCol = cell.tableColumn;

          if (cellCol == null)
            return false;

          const colName = cellCol.name;
          const value = cell.value;

          if (colName === DESCR_TITLE) {
            if (selectedByCorr.includes(value))
              ui.tooltip.show('Selected for model construction.', x, y);
            else if (selectedByPvalue.includes(value))
              ui.tooltip.show('Excluded due to a high correlation with other descriptors.', x, y);
            else
              ui.tooltip.show('Excluded due to a high p-value.', x, y);

            return true;
          } else {
            const descriptor = grid.cell(DESCR_TITLE, cell.gridRow).value;

            if (colName === WEIGHT_TITLE) {
              if (!this.desirabilityProfileRoots.has(descriptor)) {
                if (selectedByPvalue.includes(descriptor))
                  ui.tooltip.show(`No weight: <b>${descriptor}</b> is excluded due to a high correlation with other descriptors.`, x, y);
                else
                  ui.tooltip.show(`No weight: <b>${descriptor}</b> is excluded due to a high p-value.`, x, y);

                return true;
              }

              return false;
            } else {
              if (descriptorNames.includes(colName) && (!selectedByPvalue.includes(descriptor))) {
                ui.tooltip.show(`<b>${descriptor}</b> is excluded due to a high p-value; so correlation with <b>${colName}</b> is not needed.`, x, y);
                return true;
              }
            }
          }

          return false;
        }
      }
    }); // grid.onCellTooltip

    const desirabilityCol = grid.col(DESIRABILITY_COL_NAME);
    grid.setOptions({'rowHeight': STAT_GRID_HEIGHT});
    desirabilityCol!.width = DESIRABILITY_COLUMN_WIDTH;
    desirabilityCol!.cellType = 'html';

    // show desirability profile
    grid.onCellPrepare((cell) => {
      const cellCol = cell.tableColumn;
      if (cellCol == null)
        return;

      if (cell.tableColumn == null)
        return;

      if (!cell.isTableCell)
        return;

      if (cell.tableColumn.name !== DESIRABILITY_COL_NAME)
        return;

      const descriptor = grid.cell(DESCR_TITLE, cell.gridRow).value;
      const element = this.desirabilityProfileRoots.get(descriptor);

      if (element != null)
        cell.element = element;
      else {
        const selected = selectedByPvalue.includes(descriptor);
        const text = selected ? 'highly correlated with other descriptors' : 'statistically insignificant';
        const tooltipMsg = selected ?
          `No chart shown: <b>${descriptor}</b> is excluded due to a high correlation with other descriptors.` :
          `No chart shown: <b>${descriptor}</b> is excluded due to a high p-value.`;

        const divWithDescription = ui.divText(text);
        divWithDescription.style.color = COLORS.SKIPPED;
        divWithDescription.classList.add('eda-pmpo-centered-text');
        ui.tooltip.bind(divWithDescription, tooltipMsg);
        cell.element = divWithDescription;
      }
    }); // grid.onCellPrepare
  } // updateGrid

  /** Updates the main grid viewer with the pMPO scores column */
  private updateGrid(): void {
    const grid = this.view.grid;
    const name = this.predictionName;

    grid.sort([this.predictionName], [false]);

    grid.col(name)!.format = '0.0000';

    // set tooltips
    grid.onCellTooltip((cell, x, y) => {
      if (cell.isColHeader) {
        const cellCol = cell.tableColumn;
        if (cellCol) {
          if (cell.tableColumn.name === name) {
            ui.tooltip.show(getScoreTooltip(), x, y);

            return true;
          }

          return false;
        }
      }
    });
  } // updateGrid

  /** Updates the desirability profile data */
  private updateDesirabilityProfileData(descrStatsTable: DG.DataFrame, useSigmoidalCorrection: boolean): void {
    if (this.params == null)
      return;

    // Clear existing roots
    this.desirabilityProfileRoots.forEach((root) => root.remove());
    this.desirabilityProfileRoots.clear();

    const desirabilityProfile = getDesirabilityProfileJson(this.params, useSigmoidalCorrection, '', '', true);

    // Set weights
    const descrNames = descrStatsTable.col(DESCR_TITLE)!.toList();
    const weightsRaw = descrStatsTable.col(WEIGHT_TITLE)!.getRawData();
    const props = desirabilityProfile.properties;

    for (const name of Object.keys(props))
      weightsRaw[descrNames.indexOf(name)] = props[name].weight;

    // Set HTML elements
    const mpoEditor = new MpoProfileEditor();
    mpoEditor.setProfile(desirabilityProfile);
    const container = mpoEditor.root;
    const rootsCol = container.querySelector('div.d4-flex-col.ui-div');

    if (rootsCol == null)
      return;

    const rows = rootsCol.querySelectorAll('div.d4-flex-row.ui-div');

    rows.forEach((row) => {
      const children = row.children;
      if (children.length < 2) // expecting descriptor name, weight & profile
        return;

      const descrDivChildren = (children[0] as HTMLElement).children;
      if (descrDivChildren.length < 1) // expecting 1 div with descriptor name
        return;

      const descrName = (descrDivChildren[0] as HTMLElement).innerText;

      this.desirabilityProfileRoots.set(descrName, children[2] as HTMLElement);
    });
  } // updateDesirabilityProfileData

  /** Updates the ROC curve viewer with the given desirability (labels) and prediction columns
   * @return Best threshold according to Youden's J statistic
   */
  private updateRocCurve(desirability: DG.Column, prediction: DG.Column): number {
    const evaluation = getPmpoEvaluation(desirability, prediction);

    const rocDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat32Array(THRESHOLD, ROC_TRESHOLDS),
      DG.Column.fromFloat32Array(FPR_TITLE, evaluation.fpr),
      DG.Column.fromFloat32Array(TPR_TITLE, evaluation.tpr),
    ]);

    // Add baseline
    rocDf.meta.formulaLines.addLine({
      title: 'Non-informative baseline',
      formula: `\${${TPR_TITLE}} = \${${FPR_TITLE}}`,
      width: 1,
      style: 'dashed',
      min: 0,
      max: 1,
    });

    this.rocCurve.dataFrame = rocDf;
    this.rocCurve.setOptions({
      xColumnName: FPR_TITLE,
      yColumnName: TPR_TITLE,
      linesOrderColumnName: FPR_TITLE,
      linesWidth: 5,
      markerType: 'dot',
      title: `ROC Curve (AUC = ${evaluation.auc.toFixed(3)})`,
    });

    return evaluation.threshold;
  } // updateRocCurve

  /** Updates the confusion matrix viewer with the given data frame, desirability column name, and best threshold */
  private updateConfusionMatrix(df: DG.DataFrame, desColName: string, bestThreshold: number): void {
    this.confusionMatrix.dataFrame = df;
    this.confusionMatrix.setOptions({
      xColumnName: desColName,
      yColumnName: this.boolPredictionName,
      description: `Threshold: ${bestThreshold.toFixed(3)} (optimized via Youden's J)`,
      title: desColName + ' Confusion Matrix',
    });
  } // updateConfusionMatrix

  /** Fits the pMPO model to the given data and updates the viewers accordingly */
  private fitAndUpdateViewers(df: DG.DataFrame, descriptors: DG.ColumnList, desirability: DG.Column,
    pValTresh: number, r2Tresh: number, qCutoff: number, useSigmoid: boolean): void {
    const trainResult = Pmpo.fit(df, descriptors, desirability, pValTresh, r2Tresh, qCutoff);
    this.params = trainResult.params;
    const descrStatsTable = trainResult.descrStatsTable;
    const selectedByPvalue = trainResult.selectedByPvalue;
    const selectedByCorr = trainResult.selectedByCorr;

    const descriptorNames = descriptors.names();

    const prediction = Pmpo.predict(df, this.params, useSigmoid, this.predictionName);

    // Mark predictions with a color
    prediction.colors.setLinear(getOutputPalette(OPT_TYPE.MAX), {min: prediction.stats.min, max: prediction.stats.max});

    // Remove existing prediction column and add the new one
    df.columns.remove(this.predictionName);
    df.columns.add(prediction);

    // Update viewers
    this.updateGrid();

    // Update desirability profile roots map
    this.updateDesirabilityProfileData(descrStatsTable, useSigmoid);

    // Update statistics grid
    this.updateStatisticsGrid(descrStatsTable, descriptorNames, selectedByPvalue, selectedByCorr);

    // Update ROC curve
    const bestThreshold = this.updateRocCurve(desirability, prediction);

    // Update desirability prediction column
    const desColName = desirability.name;
    df.columns.remove(this.boolPredictionName);
    this.boolPredictionName = df.columns.getUnusedName(desColName + '(predicted)');
    const boolPrediction = getBoolPredictionColumn(prediction, bestThreshold, this.boolPredictionName);
    df.columns.add(boolPrediction);

    // Update confusion matrix
    this.updateConfusionMatrix(df, desColName, bestThreshold);

    this.view.dataFrame.selection.setAll(false, true);
  } // fitAndUpdateViewers

  /** Runs the pMPO model training application */
  public runTrainingApp(): void {
    const dockMng = this.view.dockManager;

    // Inputs form
    dockMng.dock(this.getInputForm(true).form, DG.DOCK_TYPE.LEFT, null, undefined, 0.1);

    // Dock viewers
    const gridNode = dockMng.findNode(this.view.grid.root);

    if (gridNode == null)
      throw new Error('Failed to train pMPO: missing a grid in the table view.');

    // Dock statistics grid
    const statGridNode = dockMng.dock(this.statGrid, DG.DOCK_TYPE.DOWN, gridNode, undefined, 0.5);

    // Dock ROC curve
    const rocNode = dockMng.dock(this.rocCurve, DG.DOCK_TYPE.RIGHT, statGridNode, undefined, 0.3);

    // Dock confusion matrix
    dockMng.dock(this.confusionMatrix, DG.DOCK_TYPE.RIGHT, rocNode, undefined, 0.2);

    this.setRibbons();
  } // runTrainingApp

  /** Runs the pMPO model training application */
  public getPmpoAppItems(): PmpoAppItems {
    return {
      statsGrid: this.statGrid,
      rocCurve: this.rocCurve,
      confusionMatrix: this.confusionMatrix,
      controls: this.getInputForm(false),
    };
  } // getViewers

  /** Creates and returns the input form for pMPO model training */
  private getInputForm(addBtn: boolean): Controls {
    const form = ui.form([]);
    form.append(ui.h2('Training data'));
    const numericColNames = this.numericCols.map((col) => col.name);

    // Function to run computations on input changes
    const runComputations = () => {
      try {
        //grok.shell.info('Running...');

        this.fitAndUpdateViewers(
          this.table,
          DG.DataFrame.fromColumns(descrInput.value).columns,
          this.table.col(desInput.value!)!,
          pInput.value!,
          rInput.value!,
          qInput.value!,
          useSigmoidInput.value,
        );
      } catch (err) {
        err instanceof PmpoError ?
          grok.shell.warning(err.message) :
          grok.shell.error(err instanceof Error ? err.message : PMPO_COMPUTE_FAILED + ': the platform issue.');
      }
    };

    // Descriptor columns input
    const descrInput = ui.input.columns('Descriptors', {
      table: this.table,
      nullable: false,
      available: numericColNames,
      checked: numericColNames,
      tooltipText: 'Descriptor columns used for model construction.',
      onValueChanged: (value) => {
        if (value != null) {
          areTunedSettingsUsed = false;
          checkAutoTuneAndRun();
        }
      },
    });
    form.append(descrInput.root);

    // Desirability column input
    const desInput = ui.input.choice('Desirability', {
      nullable: false,
      value: this.boolCols[0].name,
      items: this.boolCols.map((col) => col.name),
      tooltipText: 'Desirability column.',
      onValueChanged: (value) => {
        if (value != null) {
          areTunedSettingsUsed = false;
          checkAutoTuneAndRun();
        }
      },
    });
    form.append(desInput.root);

    const header = ui.h2('Settings');
    form.append(header);
    ui.tooltip.bind(header, 'Settings of the pMPO model.');

    // use sigmoid correction
    const useSigmoidInput = ui.input.bool('\u03C3 correction', {
      value: USE_SIGMOID_DEFAULT,
      tooltipText: 'Use the sigmoidal correction to the weighted Gaussian scores.',
      onValueChanged: (_value) => {
        areTunedSettingsUsed = false;
        checkAutoTuneAndRun();
      },
    });
    form.append(useSigmoidInput.root);

    const toUseAutoTune = (this.table.rowCount <= AUTO_TUNE_MAX_APPLICABLE_ROWS);
    const toShowAutoTuneWarning = (this.table.rowCount > AUTO_TUNE_WARNING_MIN_ROWS);

    // Flag indicating whether optimal parameters from auto-tuning are currently used
    let areTunedSettingsUsed = false;

    const setOptimalParametersAndRun = async () => {
      if (!areTunedSettingsUsed) {
        const optimalSettings = await this.getOptimalSettings(
          DG.DataFrame.fromColumns(descrInput.value).columns,
          this.table.col(desInput.value!)!,
          useSigmoidInput.value,
        );

        if (optimalSettings.success) {
          pInput.value = Math.max(optimalSettings.pValTresh, P_VAL_TRES_MIN);
          rInput.value = Math.max(optimalSettings.r2Tresh, R2_MIN);
          qInput.value = Math.max(optimalSettings.qCutoff, Q_CUTOFF_MIN);
          areTunedSettingsUsed = true;
        } else
          autoTuneInput.value = false; // revert to manual mode if optimization failed
      }

      runComputations();
    };

    const checkAutoTuneAndRun = () => {
      if (autoTuneInput.value)
        setOptimalParametersAndRun();
      else
        runComputations();
    };

    // autotuning input
    const autoTuneInput = ui.input.bool('Auto-tuning', {
      value: toUseAutoTune,
      tooltipText: 'Automatically select optimal p-value, R², and q-cutoff by maximizing AUC.',
      onValueChanged: (value) => {
        setEnability(!value);

        if (areTunedSettingsUsed)
          return;

        // If auto-tuning is turned on, set optimal parameters and run computations
        if (value) {
          if (toShowAutoTuneWarning) {
            const dlg = ui.dialog('⚠️ Long Computation')
              .add(ui.divText('Auto-tuning is time-consuming for large datasets.'))
              .add(ui.divText('Do you want to continue?'))
              .onCancel(() => autoTuneInput.value = false)
              .addButton('Run Anyway', async () => {
                dlg.close();
                setTimeout(async () => await setOptimalParametersAndRun(), 10);
              })
              .show();
          } else
            setOptimalParametersAndRun();
        }
      },
    });
    form.append(autoTuneInput.root);

    // p-value threshold input
    const pInput = ui.input.float('p-value', {
      nullable: false,
      min: P_VAL_TRES_MIN,
      max: P_VAL_TRES_MAX,
      step: 0.001,
      value: P_VAL_TRES_DEFAULT,
      // @ts-ignore
      format: FORMAT,
      tooltipText: 'P-value threshold. Descriptors with p-values above this threshold are excluded.',
      onValueChanged: (value) => {
        // Prevent running computations when auto-tuning is on, since parameters will be set automatically
        if (autoTuneInput.value)
          return;

        areTunedSettingsUsed = false;
        if ((value != null) && (value >= P_VAL_TRES_MIN) && (value <= P_VAL_TRES_MAX))
          runComputations();
      },
    });
    form.append(pInput.root);

    // R² threshold input
    const rInput = ui.input.float('R²', {
      // @ts-ignore
      format: FORMAT,
      nullable: false,
      min: R2_MIN,
      value: R2_DEFAULT,
      max: R2_MAX,
      step: 0.01,
      // eslint-disable-next-line max-len
      tooltipText: 'Squared correlation threshold. Descriptors with squared correlation above this threshold are considered highly correlated. Among them, the descriptor with the lower p-value is retained.',
      onValueChanged: (value) => {
        // Prevent running computations when auto-tuning is on, since parameters will be set automatically
        if (autoTuneInput.value)
          return;

        areTunedSettingsUsed = false;

        if ((value != null) && (value >= R2_MIN) && (value <= R2_MAX))
          runComputations();
      },
    });
    form.append(rInput.root);

    // q-cutoff input
    const qInput = ui.input.float('q-cutoff', {
      // @ts-ignore
      format: FORMAT,
      nullable: false,
      min: Q_CUTOFF_MIN,
      value: Q_CUTOFF_DEFAULT,
      max: Q_CUTOFF_MAX,
      step: 0.01,
      tooltipText: 'Q-cutoff for the pMPO model computation.',
      onValueChanged: (value) => {
        // Prevent running computations when auto-tuning is on, since parameters will be set automatically
        if (autoTuneInput.value)
          return;

        areTunedSettingsUsed = false;

        if ((value != null) && (value >= Q_CUTOFF_MIN) && (value <= Q_CUTOFF_MAX))
          runComputations();
      },
    });
    form.append(qInput.root);

    const setEnability = (toEnable: boolean) => {
      pInput.enabled = toEnable;
      rInput.enabled = toEnable;
      qInput.enabled = toEnable;
    };

    setTimeout(() => {
      if (toUseAutoTune)
        autoTuneInput.value = true; // this will trigger setting optimal parameters and running computations
      else
        runComputations();
    }, 10);

    // Save model button
    const saveBtn = ui.button('Save', async () => {
      if (this.params == null) {
        grok.shell.warning('Failed to save pMPO model: null parameters.');
        return;
      }

      saveModel(this.params, this.table.name, useSigmoidInput.value);
    }, 'Save model as platform file.');

    if (addBtn)
      form.append(saveBtn);

    const div = ui.div([form]);
    div.classList.add('eda-pmpo-input-form');

    return {
      form: div,
      saveBtn: saveBtn,
    };
  } // getInputForm

  /** Retrieves boolean columns from the data frame */
  private getBoolCols(): DG.Column[] {
    const res: DG.Column[] = [];

    for (const col of this.table.columns) {
      if ((col.type === DG.COLUMN_TYPE.BOOL) && (col.stats.stdev > 0))
        res.push(col);
    }

    return res;
  } // getBoolCols

  /** Retrieves valid (numerical, no missing values, non-zero standard deviation) numeric columns from the data frame */
  private getValidNumericCols(): DG.Column[] {
    const res: DG.Column[] = [];

    for (const col of this.table.columns) {
      if ((col.isNumerical) && (col.stats.missingValueCount < 1) && (col.stats.stdev > 0))
        res.push(col);
    }

    return res;
  } // getValidNumericCols

  /** Fits the pMPO model to the given data and updates the viewers accordingly */
  private async getOptimalSettings(descriptors: DG.ColumnList, desirability: DG.Column, useSigmoid: boolean): Promise<OptimalPoint> {
    const failedResult: OptimalPoint = {
      pValTresh: 0,
      r2Tresh: 0,
      qCutoff: 0,
      success: false,
    };

    const descriptorNames = descriptors.names();
    const {desired, nonDesired} = getDesiredTables(this.table, desirability);

    // Compute descriptors' statistics
    const descrStats = new Map<string, DescriptorStatistics>();
    descriptorNames.forEach((name) => {
      descrStats.set(name, getDescriptorStatistics(desired.col(name)!, nonDesired.col(name)!));
    });
    const descrStatsTable = getDescriptorStatisticsTable(descrStats);

    // Filter by p-value
    const selectedByPvalue = getFilteredByPvalue(descrStatsTable, P_VAL_TRES_DEFAULT);
    if (selectedByPvalue.length < 1)
      return failedResult;

    const correlationTriples = getCorrelationTriples(descriptors, selectedByPvalue);

    const funcToBeMinimized = (point: Float32Array) => {
      // Filter by correlations
      const selectedByCorr = getFilteredByCorrelations(descriptors, selectedByPvalue, descrStats, point[0], correlationTriples);

      // Compute pMPO parameters - training
      const params = getModelParams(desired, nonDesired, selectedByCorr, point[1]);

      // Get predictions
      const prediction = Pmpo.predict(this.table, params, useSigmoid, this.predictionName);

      // Evaluate predictions and return 1 - AUC (since optimization minimizes the function, but we want to maximize AUC)
      return 1 - getPmpoEvaluation(desirability, prediction).auc;
    }; // funcToBeMinimized

    const pi = DG.TaskBarProgressIndicator.create('Optimizing... ', {cancelable: true});

    try {
      const optimalResult = await optimizeNM(
        pi,
        funcToBeMinimized,
        new Float32Array([R2_DEFAULT, Q_CUTOFF_DEFAULT]),
        DEFAULT_OPTIMIZATION_SETTINGS,
        LOW_PARAMS_BOUNDS,
        HIGH_PARAMS_BOUNDS,
      );

      const success = !pi.canceled;
      pi.close();

      if (success) {
        return {
          pValTresh: P_VAL_TRES_DEFAULT,
          r2Tresh: optimalResult.optimalPoint[0],
          qCutoff: optimalResult.optimalPoint[1],
          success: true,
        };
      } else
        return failedResult;
    } catch (err) {
      pi.close();

      return failedResult;
    }
  } // getOptimalSettings
}; // Pmpo
