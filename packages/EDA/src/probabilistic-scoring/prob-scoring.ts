/* eslint-disable max-len */
// Probabilistic scoring (pMPO) features
// Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC4716604/

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

import '../../css/pmpo.css';

import {getDesiredTables, getDescriptorStatistics, gaussDesirabilityFunc, sigmoidS} from './stat-tools';
import {MIN_SAMPLES_COUNT, PMPO_NON_APPLICABLE, DescriptorStatistics, P_VAL_TRES_MIN, DESCR_TITLE,
  R2_MIN, Q_CUTOFF_MIN, PmpoParams, SCORES_TITLE, DESCR_TABLE_TITLE, PMPO_COMPUTE_FAILED, SELECTED_TITLE,
  P_VAL, DESIRABILITY_COL_NAME, STAT_GRID_HEIGHT, DESIRABILITY_COLUMN_WIDTH, WEIGHT_TITLE,
  P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT} from './pmpo-defs';
import {addSelectedDescriptorsCol, getDescriptorStatisticsTable, getFilteredByPvalue, getFilteredByCorrelations,
  getModelParams, getDescrTooltip, saveModel, getScoreTooltip, getDesirabilityProfileJson, getCorrelationTriples,
  addCorrelationColumns, setPvalColumnColorCoding, setCorrColumnColorCoding} from './pmpo-utils';
import {getOutputPalette} from '../pareto-optimization/utils';
import {OPT_TYPE} from '../pareto-optimization/defs';

export type PmpoTrainingResult = {
  params: Map<string, PmpoParams>,
  descrStatsTable: DG.DataFrame,
  selectedByPvalue: string[],
  selectedByCorr: string[],
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
    pValTresh: number, r2Tresh: number, qCutoff: number): PmpoTrainingResult {
    if (!Pmpo.isApplicable(descriptors, desirability, pValTresh, r2Tresh, qCutoff))
      throw new Error('Failed to train pMPO model: the method is not applicable to the inputs');

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
  static predict(df: DG.DataFrame, params: Map<string, PmpoParams>, predictionName: string): DG.Column {
    const count = df.rowCount;
    const scores = new Float64Array(count).fill(0);
    let x = 0;

    // Compute pMPO scores (see https://pmc.ncbi.nlm.nih.gov/articles/PMC4716604/
    params.forEach((param, name) => {
      const col = df.col(name);

      if (col == null)
        throw new Error(`Failed to apply pMPO: inconsistent data, no column "${name}" in the table "${df.name}"`);

      const vals = col.getRawData();
      for (let i = 0; i < count; ++i) {
        x = vals[i];
        scores[i] += param.weight * gaussDesirabilityFunc(x, param.desAvg, param.desStd) * sigmoidS(x, param.cutoff, param.b, param.c);
      }
    });

    return DG.Column.fromFloat64Array(predictionName, scores);
  } // predict

  private params: Map<string, PmpoParams> | null = null;

  private table: DG.DataFrame;
  private view: DG.TableView;
  private boolCols: DG.Column[];
  private numericCols: DG.Column[];

  private initTable = grok.data.demo.demog(10);

  private statGrid = DG.Viewer.grid(this.initTable, {showTitle: true, title: DESCR_TABLE_TITLE});

  private predictionName = SCORES_TITLE;

  private desirabilityProfileRoots = new Map<string, HTMLElement>();

  constructor(df: DG.DataFrame) {
    this.table = df;
    this.view = grok.shell.tableView(df.name) ?? grok.shell.addTableView(df);
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

            if ((colName === DESIRABILITY_COL_NAME) || (colName === WEIGHT_TITLE)) {
              const startText = (colName === WEIGHT_TITLE) ? 'No weight' : 'No chart shown';

              if (!this.desirabilityProfileRoots.has(descriptor)) {
                if (selectedByPvalue.includes(descriptor))
                  ui.tooltip.show(`${startText}: <b>${descriptor}</b> is excluded due to a high correlation with other descriptors.`, x, y);
                else
                  ui.tooltip.show(`${startText}: <b>${descriptor}</b> is excluded due to a high p-value.`, x, y);

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
      cell.element = this.desirabilityProfileRoots.get(descriptor) ?? ui.div();
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
  private updateDesirabilityProfileData(descrStatsTable: DG.DataFrame): void {
    if (this.params == null)
      return;

    // Clear existing roots
    this.desirabilityProfileRoots.forEach((root) => root.remove());
    this.desirabilityProfileRoots.clear();

    const desirabilityProfile = getDesirabilityProfileJson(this.params, '', '');

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

  /** Fits the pMPO model to the given data and updates the viewers accordingly */
  private fitAndUpdateViewers(df: DG.DataFrame, descriptors: DG.ColumnList, desirability: DG.Column,
    pValTresh: number, r2Tresh: number, qCutoff: number): void {
    const trainResult = Pmpo.fit(df, descriptors, desirability, pValTresh, r2Tresh, qCutoff);
    this.params = trainResult.params;
    const descrStatsTable = trainResult.descrStatsTable;
    const selectedByPvalue = trainResult.selectedByPvalue;
    const selectedByCorr = trainResult.selectedByCorr;

    const descriptorNames = descriptors.names();

    //const weightsTable = getWeightsTable(this.params);
    const prediction = Pmpo.predict(df, this.params, this.predictionName);

    // Mark predictions with a color
    prediction.colors.setLinear(getOutputPalette(OPT_TYPE.MAX), {min: prediction.stats.min, max: prediction.stats.max});

    df.columns.remove(this.predictionName);
    df.columns.add(prediction);

    // Update viewers
    this.updateGrid();

    // Update desirability profile roots map
    this.updateDesirabilityProfileData(descrStatsTable);

    // Update statistics grid
    this.updateStatisticsGrid(descrStatsTable, descriptorNames, selectedByPvalue, selectedByCorr);
  } // fitAndUpdateViewers

  /** Runs the pMPO model training application */
  public runTrainingApp(): void {
    const dockMng = this.view.dockManager;

    // Inputs form
    dockMng.dock(this.getInputForm(), DG.DOCK_TYPE.LEFT, null, undefined, 0.1);

    // Dock viewers
    const gridNode = dockMng.findNode(this.view.grid.root);

    if (gridNode == null)
      throw new Error('Failed to train pMPO: missing a grid in the table view.');

    dockMng.dock(this.statGrid, DG.DOCK_TYPE.DOWN, gridNode, undefined, 0.5);

    this.setRibbons();
  } // runTrainingApp

  /** Creates and returns the input form for pMPO model training */
  private getInputForm(): HTMLElement {
    const form = ui.form([]);
    form.append(ui.h2('Training data'));
    const numericColNames = this.numericCols.map((col) => col.name);

    // Function to run computations on input changes
    const runComputations = () => {
      try {
        this.fitAndUpdateViewers(
          this.table,
          DG.DataFrame.fromColumns(descrInput.value).columns,
        this.table.col(desInput.value!)!,
        pInput.value!,
        rInput.value!,
        qInput.value!,
        );
      } catch (err) {
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
        if (value != null)
          runComputations();
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
        if (value != null)
          runComputations();
      },
    });
    form.append(desInput.root);

    const header = ui.h2('Thresholds');
    ui.tooltip.bind(header, 'Settings of the pMPO model training.');
    form.append(header);

    // p-value threshold input
    const pInput = ui.input.float('p-value', {
      nullable: false,
      min: P_VAL_TRES_MIN,
      max: 1,
      step: 0.01,
      value: P_VAL_TRES_DEFAULT,
      tooltipText: 'Descriptors with p-values above this threshold are excluded.',
      onValueChanged: (value) => {
        if ((value != null) && (value >= P_VAL_TRES_MIN) && (value <= 1))
          runComputations();
      },
    });
    form.append(pInput.root);

    // R² threshold input
    const rInput = ui.input.float('R²', {
      nullable: false,
      min: R2_MIN,
      value: R2_DEFAULT,
      max: 1,
      step: 0.01,
      // eslint-disable-next-line max-len
      tooltipText: 'Descriptors with squared correlation above this threshold are considered highly correlated. Among them, the descriptor with the lower p-value is retained.',
      onValueChanged: (value) => {
        if ((value != null) && (value >= R2_MIN) && (value <= 1))
          runComputations();
      },
    });
    form.append(rInput.root);

    // q-cutoff input
    const qInput = ui.input.float('q-cutoff', {
      nullable: false,
      min: Q_CUTOFF_MIN,
      value: Q_CUTOFF_DEFAULT,
      max: 1,
      step: 0.01,
      tooltipText: 'Q-cutoff for the pMPO model computation.',
      onValueChanged: (value) => {
        if ((value != null) && (value >= Q_CUTOFF_MIN) && (value <= 1))
          runComputations();
      },
    });
    form.append(qInput.root);

    setTimeout(() => runComputations(), 10);

    // Save model button
    const saveBtn = ui.button('Save model', async () => {
      if (this.params == null) {
        grok.shell.warning('Failed to save pMPO model: null parameters.');
        return;
      }

      saveModel(this.params, this.table.name);
    }, 'Save model as platform file.');
    form.append(saveBtn);

    const div = ui.div([form]);
    div.classList.add('eda-pmpo-input-form');

    return div;
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
}; // Pmpo
