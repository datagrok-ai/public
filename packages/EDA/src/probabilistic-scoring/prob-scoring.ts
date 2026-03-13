/* eslint-disable max-len */
// Probabilistic scoring (pMPO) features
// Source paper https://pmc.ncbi.nlm.nih.gov/articles/PMC4716604/

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

import '../../css/pmpo.css';

import {getDesiredTables, getDescriptorStatistics, getBoolPredictionColumn, getPmpoEvaluation} from './stat-tools';
import {MIN_SAMPLES_COUNT, PMPO_NON_APPLICABLE, DescriptorStatistics, P_VAL_TRES_MIN, DESCR_TITLE,
  R2_MIN, Q_CUTOFF_MIN, PmpoParams, SCORES_TITLE, DESCR_TABLE_TITLE, SELECTED_TITLE,
  P_VAL, DESIRABILITY_COL_NAME, STAT_GRID_HEIGHT, DESIRABILITY_COLUMN_WIDTH, WEIGHT_TITLE,
  P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT, USE_SIGMOID_DEFAULT, ROC_TRESHOLDS,
  FPR_TITLE, TPR_TITLE, COLORS, THRESHOLD, AUTO_TUNE_MAX_APPLICABLE_ROWS, DEFAULT_OPTIMIZATION_SETTINGS,
  P_VAL_TRES_MAX, R2_MAX, Q_CUTOFF_MAX, OptimalPoint, LOW_PARAMS_BOUNDS, HIGH_PARAMS_BOUNDS, FORMAT,
  EQUALITY_SIGN, SIGN_OPTIONS, THRESHOLDED_DESIRABILITY_COL_NAME, PMPO_COMPUTE_FAILED,
  PmpoInputId, TooltipContent, PmpoValidationResult} from './pmpo-defs';
import {addSelectedDescriptorsCol, getDescriptorStatisticsTable, getFilteredByPvalue, getFilteredByCorrelations,
  getModelParams, getDescrTooltip, saveModel, getScoreTooltip, getDesirabilityProfileJson, getCorrelationTriples,
  addCorrelationColumns, setPvalColumnColorCoding, setCorrColumnColorCoding, PmpoError, getInitCol,
  getBoolDesirabilityColData, isDesirabilityValid,
  getDesirabilityColumnFromCategories,
  getSelectedCategories} from './pmpo-utils';
import {getOutputPalette} from '../pareto-optimization/utils';
import {OPT_TYPE} from '../pareto-optimization/defs';
import {optimizeNM} from './nelder-mead';
import {getMissingValsIndices} from '../missing-values-imputation/knn-imputer';
import {DesirabilityProfile} from '@datagrok-libraries/statistics/src/mpo/mpo';

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
  profile: DesirabilityProfile | null;
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

      if (col.stats.missingValueCount === col.length) {
        showWarning(`: "${col.name}" contains only missing values.`);
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

    let validColsCount = 0;

    // Check numeric columns and boolean columns
    for (const col of df.columns) {
      if (col.isNumerical || (col.type === DG.TYPE.BOOL)) {
        if (col.stats.stdev > 0)
          ++validColsCount;
      }
    }

    // Check valid numeric columns count
    if (validColsCount < 2) {
      if (toShowMsg)
        grok.shell.warning(PMPO_NON_APPLICABLE + ': not enough of non-constant columns.');
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
    const colsWithMissingVals: DG.Column[] = [];

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

      if (col.stats.missingValueCount > 0)
        colsWithMissingVals.push(col);

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
  private desirabilityProfile: DesirabilityProfile | null = null;

  private table: DG.DataFrame;
  private view: DG.TableView;
  private desirabilityColumns: DG.Column[];
  private numericCols: DG.Column[];
  private missingValsIndeces: Map<string, number[]>;

  private initTable = DG.DataFrame.create();

  private statGrid = DG.Viewer.grid(this.initTable, {showTitle: true, title: DESCR_TABLE_TITLE});

  private predictionName = SCORES_TITLE;
  private boolPredictionName = '';

  private desirabilityProfileRoots = new Map<string, HTMLElement>();

  private tresholdedColumn: DG.Column | null = null;
  private threshColTooltip: string | null = null;

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
    this.desirabilityColumns = this.getDesirabilityColumns();
    this.numericCols = this.getValidNumericCols();
    this.predictionName = df.columns.getUnusedName(SCORES_TITLE);
    this.missingValsIndeces = getMissingValsIndices(this.numericCols);
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

    const scoresCol = grid.col(name);
    scoresCol!.format = '0.0000';
    scoresCol!.isTextColorCoded = true;

    // set tooltips
    grid.onCellTooltip((cell, x, y) => {
      if (cell.isColHeader) {
        const cellCol = cell.tableColumn;
        if (cellCol) {
          if (cell.tableColumn.name === name) {
            ui.tooltip.show(getScoreTooltip(), x, y);

            return true;
          } else {
            if (this.tresholdedColumn != null && cell.tableColumn.name === this.tresholdedColumn.name) {
              ui.tooltip.show(ui.markdown(this.threshColTooltip ?? ''), x, y);

              return true;
            }
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
    this.desirabilityProfile = getDesirabilityProfileJson(this.params, useSigmoidalCorrection, '', '', false);

    // Set weights
    const descrNames = descrStatsTable.col(DESCR_TITLE)!.toList();
    const weightsRaw = descrStatsTable.col(WEIGHT_TITLE)!.getRawData();
    const props = desirabilityProfile.properties;
    const names: string[] = Object.keys(props);

    for (const name of names)
      weightsRaw[descrNames.indexOf(name)] = props[name].weight;

    // Set HTML elements
    const mpoEditor = new MpoProfileEditor();
    mpoEditor.setProfile(desirabilityProfile);
    const container = mpoEditor.root;
    const rootsCol = container.querySelector('div.d4-flex-col.ui-div');

    if (rootsCol == null)
      return;

    const rows = rootsCol.querySelectorAll('div.d4-flex-row.ui-div.statistics-mpo-row');

    rows.forEach((row, idx) => {
      const children = row.children;
      if (children.length < 2) // expecting descriptor name, weight & profile
        return;

      const profileRoot = children[2] as HTMLElement;
      profileRoot.style.width = '100%';
      this.desirabilityProfileRoots.set(names[idx], profileRoot);
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

  /** Sets null values for the predicted scores in rows with missing values in any of the descriptors */
  private getIndecesOfMissingValues(colNames: string[]): number[] {
    const indeces: number[] = [];

    colNames.forEach((name) => {
      const inds = this.missingValsIndeces.get(name);

      if (inds != null)
        indeces.push(...inds);
    });

    return indeces;
  }

  /** Sets null values for the predicted scores in rows with missing values in any of the descriptors */
  private setNulls(scores: DG.Column, indeces: number[]): void {
    const raw = scores.getRawData();
    indeces.forEach((ind) => raw[ind] = DG.FLOAT_NULL);
  }

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

    // Set nulls for rows with missing values in any of the selected descriptors
    const indecesOfMissingVals = this.getIndecesOfMissingValues(selectedByCorr);
    this.setNulls(prediction, indecesOfMissingVals);

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
      profile: this.desirabilityProfile,
    };
  } // getViewers

  /** Creates and returns the input form for pMPO model training */
  private getInputForm(addBtn: boolean): Controls {
    const form = ui.form([]);
    form.append(ui.h2('Training data'));
    const initDesirability = getInitCol(this.desirabilityColumns);

    // returns the desirability column to be used for computations, based on the input desirability column and threshold settings
    const getDesirabilityColumn = (): DG.Column => {
      // remove existing thresholded column if exists
      if (this.tresholdedColumn != null) {
        this.table.columns.remove(this.tresholdedColumn.name);
        this.tresholdedColumn = null;
      }

      if (desInput.value!.type === DG.COLUMN_TYPE.BOOL)
        return desInput.value!;

      const boolDesirabilityData = (desInput.value!.type === DG.COLUMN_TYPE.STRING) ?
        getDesirabilityColumnFromCategories(desInput.value!, desirableCategoriesInput!.value!) :
        getBoolDesirabilityColData(
          desInput.value!,
          desirabilityThresholdInput.value!,
          signInput.value as EQUALITY_SIGN,
        );

      this.tresholdedColumn = boolDesirabilityData.column;
      this.threshColTooltip = boolDesirabilityData.tooltip;

      this.tresholdedColumn.name = this.table.columns.getUnusedName(THRESHOLDED_DESIRABILITY_COL_NAME);
      this.table.columns.add(this.tresholdedColumn);

      return this.tresholdedColumn;
    }; // getDesirabilityColumn

    // Function to run computations on input changes
    const runComputations = () => {
      if (!areInputsValid())
        return;

      try {
        this.fitAndUpdateViewers(
          this.table,
          DG.DataFrame.fromColumns(descrInput.value).columns,
          getDesirabilityColumn(),
          pInput.value!,
          rInput.value!,
          qInput.value!,
          useSigmoidInput.value,
        );
      } catch (err) {
        if (err instanceof PmpoError) {
          grok.shell.warning(err.message);
          ui.tooltip.bind(desInput.input, err.message);
          ui.tooltip.bind(descrInput.input, err.message);
        } else {
          const msg = err instanceof Error ? err.message : PMPO_COMPUTE_FAILED + ': the platform issue.';
          grok.shell.error(msg);
          ui.tooltip.bind(desInput.input, msg);
          ui.tooltip.bind(descrInput.input, msg);
        };

        desInput.input.classList.add('d4-invalid');
        descrInput.input.classList.add('d4-invalid');
      }
    }; // runComputations

    // Descriptor columns input
    const descrInput = ui.input.columns('Descriptors', {
      table: this.table,
      nullable: false,
      available: this.numericCols.map((col) => col.name),
      checked: this.numericCols.filter((col) => {
        return (col.name !== initDesirability.name) && (col.stats.stdev > 0) && (col.stats.missingValueCount < col.length);
      }).map((col) => col.name),
      tooltipText: 'Descriptor columns used for model construction.',
      onValueChanged: (value) => {
        if (value != null) {
          areTunedSettingsUsed = false;
          checkAutoTuneAndRun();
        }
      },
    });
    form.append(descrInput.root);

    descrInput.addValidator(() => {
      if (descrInput.value == null || descrInput.value.length < 1)
        return 'Select at least one descriptor column.';
      if (desInput.value != null && descrInput.value.includes(desInput.value))
        return 'Desirability column cannot be used as a descriptor.';
      const zeroStdevCols = descrInput.value.filter((col) => col.stats.stdev === 0).map((col) => col.name);
      if (zeroStdevCols.length > 0)
        return `Descriptor columns with zero variance: ${zeroStdevCols.join(', ')}`;
      const nullCols = descrInput.value.filter((col) => col.stats.missingValueCount === col.length).map((col) => col.name);
      if (nullCols.length > 0)
        return `Descriptor columns with only missing values: ${nullCols.join(', ')}`;
      return null;
    });

    // Desirability column input and related controls
    const setVisibilityOfDesirabilityAuxInputs = (value: DG.Column) => {
      if (value.type === DG.COLUMN_TYPE.BOOL)
        desOptionsInputDiv.hidden = true;
      else {
        desOptionsInputDiv.hidden = false;
        const isString = (value.type === DG.COLUMN_TYPE.STRING);
        desirabilityThresholdInput.root.hidden = isString;
        signInput.root.hidden = isString;
      }
    }; // setVisibilityOfDesirabilityAuxInputs

    const desInput = ui.input.column('Desirability', {
      nullable: false,
      value: initDesirability,
      table: this.table,
      filter: (col) => this.desirabilityColumns.includes(col),
      tooltipText: 'Desirability column.',
      onValueChanged: (value) => {
        if (value != null) {
          updateDesirableCategoriesInput();
          setVisibilityOfDesirabilityAuxInputs(value);
          areComputationsBlocked = true;
          desirabilityThresholdInput.value = Math.round(value.stats.avg * 100) / 100;
          areComputationsBlocked = false;
          areTunedSettingsUsed = false;
          checkAutoTuneAndRun();
        }
      }, // onValueChanged
    });
    form.append(desInput.root);

    desInput.addValidator(() => {
      if (desInput.value == null)
        return 'Select a desirability column.';
      if (descrInput.value != null && descrInput.value.includes(desInput.value))
        return 'Desirability column cannot be used as a descriptor.';
      if (desInput.value.type === DG.COLUMN_TYPE.BOOL) {
        if (desInput.value.stats.stdev === 0)
          return 'All desirability values are the same - scoring is not feasible.';
      } else if (desInput.value.type === DG.COLUMN_TYPE.STRING) {
        if (desInput.value.categories.length < 2)
          return 'String desirability column must have at least 2 categories.';
      } else {
        if (desInput.value.stats.stdev === 0) {
          return desInput.value.stats.missingValueCount < desInput.value.length ?
            'All desirability values are the same - scoring is not feasible.' :
            'Empty column cannot be used as desirability column.';
        }
        if (desirabilityThresholdInput.value == null)
          return 'Specify non-null desirability threshold.';
        if (!isDesirabilityValid(desInput.value, desirabilityThresholdInput.value, signInput.value as EQUALITY_SIGN)) {
          return `All compounds are either desired or non-desired for ${desInput.value.name} ` +
            `${signInput.value} ${desirabilityThresholdInput.value}. Adjust the threshold or condition.`;
        }
      }
      return null;
    });

    let areComputationsBlocked = false;

    const signInput = ui.input.choice('Condition', {
      value: EQUALITY_SIGN.DEFAULT,
      items: SIGN_OPTIONS,
      nullable: false,
      tooltipText: 'How to compare numeric Desirability column values against the threshold.',
      onValueChanged: (_value) => {
        areTunedSettingsUsed = false;
        checkAutoTuneAndRun();
      },
    });

    const desirabilityThresholdInput = ui.input.float('Threshold', {
      value: Math.round(initDesirability.stats.avg * 100) / 100,
      nullable: false,
      tooltipText: 'Boundary value that separates desired from non-desired compounds.',
      format: '0.00',
      onValueChanged: (value) => {
        if (value != null) {
          if (areComputationsBlocked)
            return;
          areTunedSettingsUsed = false;
          checkAutoTuneAndRun();
        }
      },
    });

    desirabilityThresholdInput.addValidator(() => {
      if (desInput.value == null || desInput.value.type === DG.COLUMN_TYPE.BOOL ||
        desInput.value.type === DG.COLUMN_TYPE.STRING)
        return null;
      if (desirabilityThresholdInput.value == null)
        return 'Specify non-null desirability threshold.';
      if (!isDesirabilityValid(desInput.value, desirabilityThresholdInput.value, signInput.value as EQUALITY_SIGN))
        return 'Adjust the threshold to get both desired and non-desired groups.';
      return null;
    });

    const desOptionsInputDiv = ui.divV([signInput.root, desirabilityThresholdInput.root]);

    form.append(desOptionsInputDiv);

    let desirableCategoriesInput: DG.InputBase<string[] | null> | null = null;

    // For string columns - input for selecting which categories are considered desirable
    const updateDesirableCategoriesInput = () => {
      if (desirableCategoriesInput != null) {
        desirableCategoriesInput.root.remove();
        desirableCategoriesInput = null;
      }

      if (desInput.value?.type === DG.COLUMN_TYPE.STRING) {
        desirableCategoriesInput = ui.input.multiChoice('Preferred', {
          value: getSelectedCategories(desInput.value!.categories),
          items: desInput.value!.categories,
          nullable: false,
          tooltipText: 'Select which categories should be treated as desirable.',
          onValueChanged: (value) => {
            if (value != null) {
              if (areComputationsBlocked)
                return;
              areTunedSettingsUsed = false;
              checkAutoTuneAndRun();
            }
          },
        });

        desirableCategoriesInput.addValidator(() => {
          if (desirableCategoriesInput!.value == null || desirableCategoriesInput!.value.length === 0)
            return 'Select at least one preferable category.';
          if (desInput.value != null && desirableCategoriesInput!.value.length === desInput.value.categories.length)
            return 'At least one category must be non-preferable.';
          return null;
        });

        desOptionsInputDiv.append(desirableCategoriesInput.root);
      }
    }; // updateDesirableCategoriesInput

    setVisibilityOfDesirabilityAuxInputs(desInput.value!);

    // Settings inputs

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

    // Flag indicating whether optimal parameters from auto-tuning are currently used
    let areTunedSettingsUsed = false;

    // Auto-tune parameters and run computations; if auto-tune is not applicable, just run computations with current settings
    const setOptimalParametersAndRun = async () => {
      await new Promise((resolve) => setTimeout(resolve, 50));

      if (!areInputsValid())
        return;

      if (!areTunedSettingsUsed) {
        const optimalSettings = await this.getOptimalSettings(
          DG.DataFrame.fromColumns(descrInput.value).columns,
          getDesirabilityColumn(),
          useSigmoidInput.value,
        );

        if (optimalSettings.state === 'success') {
          pInput.value = Math.max(optimalSettings.pValTresh, P_VAL_TRES_MIN);
          rInput.value = Math.max(optimalSettings.r2Tresh, R2_MIN);
          qInput.value = Math.max(optimalSettings.qCutoff, Q_CUTOFF_MIN);
          areTunedSettingsUsed = true;
          runComputations();
        } else
          grok.shell.warning(optimalSettings.msg);
          /*descrInput.input.classList.add('d4-invalid');
          desInput.input.classList.add('d4-invalid');
          ui.tooltip.bind(descrInput.input, optimalSettings.msg);
          ui.tooltip.bind(desInput.input, optimalSettings.msg);*/
      } else
        runComputations();
    }; // setOptimalParametersAndRun

    // Validates all inputs before running computations using registered validators
    const areInputsValid = (): boolean => {
      const results = [
        descrInput.validate(),
        desInput.validate(),
        desirabilityThresholdInput.validate(),
        pInput.validate(),
        rInput.validate(),
        qInput.validate(),
      ];

      if (desirableCategoriesInput != null)
        results.push(desirableCategoriesInput.validate());

      return results.every((r) => r);
    }; // areInputsValid

    const checkAutoTuneAndRun = () => {
      if (autoTuneInput.value)
        setOptimalParametersAndRun();
      else
        runComputations();
    };

    // autotuning input
    const autoTuneInput = ui.input.bool('Auto-tuning', {
      value: false,
      tooltipText: 'Automatically select optimal p-value, R², and q-cutoff by maximizing AUC.',
      onValueChanged: async (value) => {
        setEnability(!value);

        if (areTunedSettingsUsed)
          return;

        // If auto-tuning is turned on, set optimal parameters and run computations
        if (value)
          await setOptimalParametersAndRun();
        else
          runComputations();
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

    pInput.addValidator(() => {
      if (pInput.value == null)
        return 'P-value is required.';
      if (pInput.value < P_VAL_TRES_MIN || pInput.value > P_VAL_TRES_MAX)
        return `P-value must be between ${P_VAL_TRES_MIN} and ${P_VAL_TRES_MAX}.`;
      return null;
    });

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

    rInput.addValidator(() => {
      if (rInput.value == null)
        return 'R² is required.';
      if (rInput.value < R2_MIN || rInput.value > R2_MAX)
        return `R² must be between ${R2_MIN} and ${R2_MAX}.`;
      return null;
    });

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

    qInput.addValidator(() => {
      if (qInput.value == null)
        return 'Q-cutoff is required.';
      if (qInput.value < Q_CUTOFF_MIN || qInput.value > Q_CUTOFF_MAX)
        return `Q-cutoff must be between ${Q_CUTOFF_MIN} and ${Q_CUTOFF_MAX}.`;
      return null;
    });

    const setEnability = (toEnable: boolean) => {
      pInput.enabled = toEnable;
      rInput.enabled = toEnable;
      qInput.enabled = toEnable;
    };

    setTimeout(() => {
      runComputations();

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

  /** Validates all pMPO inputs and returns structured errors without mutating the DOM */
  static validateInputs(params: {
    descriptors: DG.Column[] | null,
    desirability: DG.Column | null,
    threshold: number | null,
    sign: EQUALITY_SIGN,
    desirableCategories: string[] | null,
    pValue: number | null,
    r2: number | null,
    qCutoff: number | null,
  }): PmpoValidationResult {
    const errors = new Map<PmpoInputId, TooltipContent>();
    const {descriptors, desirability, threshold, sign, desirableCategories, pValue, r2, qCutoff} = params;

    // Settings null or out of range
    if (pValue == null || r2 == null || qCutoff == null)
      return {valid: false, errors};

    if ((pValue <= 0) || (pValue > 1) || (r2 < 0) || (r2 > 1) || (qCutoff <= 0) || (qCutoff > 1))
      return {valid: false, errors};

    // Column inputs null
    if (descriptors == null || desirability == null)
      return {valid: false, errors};

    // At least one descriptor
    if (descriptors.length < 1) {
      errors.set('descriptors', 'Select at least one descriptor column.');
      return {valid: false, errors};
    }

    // Desirability column must not be among descriptors
    if (descriptors.includes(desirability)) {
      const msg = 'Desirability column cannot be used as a descriptor.';
      errors.set('descriptors', msg);
      errors.set('desirability', msg);
      return {valid: false, errors};
    }

    // No zero-variance descriptor columns
    const zeroStdevCols = descriptors.filter((col) => col.stats.stdev === 0).map((col) => col.name);
    if (zeroStdevCols.length > 0)
      errors.set('descriptors', () => ui.markdown(`Descriptor columns with zero variance cannot be used: **${zeroStdevCols.join(', ')}**`));

    // No all-null descriptor columns
    const nullCols = descriptors.filter((col) => col.stats.missingValueCount === col.length).map((col) => col.name);
    if (nullCols.length > 0)
      errors.set('descriptors', () => ui.markdown(`Descriptor columns with only missing values cannot be used: **${nullCols.join(', ')}**`));

    // Validate desirability column based on its type
    if (desirability.type === DG.COLUMN_TYPE.BOOL) {
      if (desirability.stats.stdev === 0)
        errors.set('desirability', 'All desirability values are the same - scoring is not feasible.');
    } else if (desirability.type === DG.COLUMN_TYPE.STRING) {
      const catsCount = desirability.categories.length;
      const selectedCatsCount = desirableCategories?.length ?? 0;

      if (catsCount < 2)
        errors.set('desirability', 'String desirability column must have at least 2 categories.');
      else if (selectedCatsCount === 0)
        errors.set('desirability', 'Select at least one preferable category.');
      else if (selectedCatsCount === catsCount)
        errors.set('desirability', 'At least one category must be non-preferable.');
    } else {
      // Numeric desirability
      if (desirability.stats.stdev === 0) {
        errors.set('desirability',
          desirability.stats.missingValueCount < desirability.length ?
            'All desirability values are the same - scoring is not feasible.' :
            'Empty column cannot be used as desirability column.',
        );
      } else if (threshold == null)
        errors.set('desirability', 'Specify non-null desirability threshold.');
      else if (!isDesirabilityValid(desirability, threshold, sign)) {
        errors.set('desirability', () => ui.markdown(`All compounds are either desired or non-desired for
          <div align="center">
          **${desirability.name} ${sign} ${threshold}.**
          </div>
          Adjust the threshold or condition to get both groups.`));
        errors.set('threshold', 'Adjust the threshold to get both desired and non-desired groups.');
      }
    }

    return {valid: !errors.size, errors};
  } // validateInputs

  /** Retrieves acceptable desirability columns (boolean or numerical with non-zero standard deviation) from the data frame */
  private getDesirabilityColumns(): DG.Column[] {
    const res: DG.Column[] = [];

    for (const col of this.table.columns) {
      if (((col.type === DG.COLUMN_TYPE.BOOL) || (col.isNumerical) || (col.type === DG.COLUMN_TYPE.STRING)))
        res.push(col);
    }

    return res;
  } // getDesirabilityColumns

  /** Retrieves valid (numerical, no missing values, non-zero standard deviation) numeric columns from the data frame */
  private getValidNumericCols(): DG.Column[] {
    const res: DG.Column[] = [];

    for (const col of this.table.columns) {
      if (col.isNumerical)
        res.push(col);
    }

    return res;
  } // getValidNumericCols

  /** Fits the pMPO model to the given data and updates the viewers accordingly */
  private async getOptimalSettings(descriptors: DG.ColumnList, desirability: DG.Column, useSigmoid: boolean): Promise<OptimalPoint> {
    const pi = DG.TaskBarProgressIndicator.create('Optimizing... ', {cancelable: true});

    try {
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
      if (selectedByPvalue.length < 1) {
        pi.close();

        return {
          pValTresh: 0,
          r2Tresh: 0,
          qCutoff: 0,
          state: 'failed',
          msg: 'No descriptors passed the p-value threshold filter.',
        };
      }

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
          state: 'success',
          msg: 'Optimization completed successfully.',
        };
      } else {
        return {
          pValTresh: 0,
          r2Tresh: 0,
          qCutoff: 0,
          state: 'canceled',
          msg: 'Auto-tuning was canceled by the user.',
        };
      }
    } catch (err) {
      pi.close();

      return {
        pValTresh: 0,
        r2Tresh: 0,
        qCutoff: 0,
        state: 'failed',
        msg: err instanceof Error ? err.message : 'Optimization failed due to an unexpected error.',
      };
    }
  } // getOptimalSettings
}; // Pmpo
