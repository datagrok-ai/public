import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/pmpo.css';

import {getDesiredTables, getDescriptorStatistics, normalPdf, sigmoidS} from './stat-tools';
import {MIN_SAMPLES_COUNT, PMPO_NON_APPLICABLE, DescriptorStatistics, P_VAL_TRES_MIN, DESCR_TITLE,
  R2_MIN, Q_CUTOFF_MIN, WEIGHT_TABLE_TITLE, PmpoParams, SCORES_TITLE, DESCR_TABLE_TITLE,
  PMPO_COMPUTE_FAILED, SELECTED_TITLE, P_VAL} from './pmpo-defs';
import {addSelectedDescriptorsCol, getDescriptorStatisticsTable, getFilteredByPvalue, getFilteredByCorrelations,
  getModelParams, getWeightsTable, getDescrTooltip, saveModel} from './pmpo-utils';
import {getOutputPalette} from '../pareto-optimization/utils';
import {OPT_TYPE} from '../pareto-optimization/defs';


export class Pmpo {
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
  }

  static isTableValid(df: DG.DataFrame, toShowMsg: boolean = true): boolean {
    if (df.rowCount < 2) {
      if (toShowMsg)
        grok.shell.warning(PMPO_NON_APPLICABLE + `. Not enough of samples: ${df.rowCount}, minimum: 2.`);
      return false;
    }

    let boolColsCount = 0;
    let validNumericColsCount = 0;

    for (const col of df.columns) {
      if (col.isNumerical) {
        if ((col.stats.missingValueCount < 1) && (col.stats.stdev > 0))
          ++validNumericColsCount;
      } else if (col.type == DG.COLUMN_TYPE.BOOL)
        ++boolColsCount;
    }

    if (boolColsCount < 1) {
      if (toShowMsg)
        grok.shell.warning(PMPO_NON_APPLICABLE + ': no boolean columns.');
      return false;
    }

    if (validNumericColsCount < 1) {
      if (toShowMsg)
        grok.shell.warning(PMPO_NON_APPLICABLE + ': no numeric columns without missing values and non-zero variance.');
      return false;
    }

    return true;
  } // isTableValid

  static predict(df: DG.DataFrame, params: Map<string, PmpoParams>, predictionName: string): DG.Column {
    const count = df.rowCount;
    const scores = new Float64Array(count).fill(0);
    let x = 0;

    params.forEach((param, name) => {
      const col = df.col(name);

      if (col == null)
        throw new Error(`Filed to apply pMPO: inconsistent data, no column "${name}" in the table "${df.name}"`);

      const vals = col.getRawData();
      for (let i = 0; i < count; ++i) {
        x = vals[i];
        scores[i] += param.weight * normalPdf(x, param.desAvg, param.desStd) * sigmoidS(x, param.x0, param.b, param.c);
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
  private corrPlot = DG.Viewer.correlationPlot(this.initTable, {showTitle: true, title: `${DESCR_TITLE} Correlations`});
  private barChart = DG.Viewer.barChart(this.initTable, {
    valueAggrType: DG.AGG.AVG,
    showCategorySelector: false,
    showValueSelector: false,
    showTitle: true,
    title: WEIGHT_TABLE_TITLE,
    showStackSelector: false,
  });
  private boxPlot = DG.Viewer.boxPlot(this.initTable, {
    showPValue: false,
    showSizeSelector: false,
    showTitle: true,
    showColorSelector: false,
    legendVisibility: 'Always',
  });
  private hist = DG.Viewer.histogram(this.initTable, {
    normalizeValues: false,
    splitStack: true,
    showSplitSelector: false,
    title: 'Distribution of pMPO scores',
    showTitle: true,
  });

  private predictionName = SCORES_TITLE;

  constructor(df: DG.DataFrame) {
    this.table = df;
    this.view = grok.shell.tableView(df.name) ?? grok.shell.addTableView(df);
    this.boolCols = this.getBoolCols();
    this.numericCols = this.getValidNumericCols();
    this.predictionName = df.columns.getUnusedName(SCORES_TITLE);
  };

  private updateGrid(table: DG.DataFrame, selectedByPvalue: string[], selectedByCorr: string[]): void {
    this.statGrid.dataFrame = table;
    this.statGrid.setOptions({
      showTitle: true,
      title: table.name,
    });

    this.statGrid.sort([SELECTED_TITLE, P_VAL], [false, true]);
        this.statGrid.col(P_VAL)!.format = 'scientific';

        // set tooltips
        this.statGrid.onCellTooltip(function(cell, x, y) {
          if (cell.isColHeader) {
            const cellCol = cell.tableColumn;
            if (cellCol) {
              if (cell.tableColumn.name === DESCR_TITLE) {
                ui.tooltip.show(getDescrTooltip(), x, y);

                return true;
              }

              return false;
            }
          } else {
            if (cell.isTableCell) {
              const cellCol = cell.tableColumn;
              if (cellCol) {
                if (cell.tableColumn.name === DESCR_TITLE) {
                  const value = cell.value;
                  if (selectedByCorr.includes(value))
                    ui.tooltip.show('Selected for model construction.', x, y);
                  else if (selectedByPvalue.includes(value))
                    ui.tooltip.show('Excluded due to a high correlation with other descriptors.', x, y);
                  else
                    ui.tooltip.show('Excluded due to a high p-value.', x, y);

                  return true;
                }

                return false;
              }
            }
          }
        });
  }

  private fitAndUpdateViewers(df: DG.DataFrame, descriptors: DG.ColumnList, desirability: DG.Column,
    pValTresh: number, r2Tresh: number, qCutoff: number): void {
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

    // Filter by p-value
    const selectedByPvalue = getFilteredByPvalue(descrStatsTable, pValTresh);

    // Filter by correlations
    const selectedByCorr = getFilteredByCorrelations(descriptors, selectedByPvalue, descrStats, r2Tresh);

    addSelectedDescriptorsCol(descrStatsTable, selectedByCorr);

    // Compute pMPO parameters - training
    this.params = getModelParams(desired, nonDesired, selectedByCorr, qCutoff);

    const weightsTable = getWeightsTable(this.params);
    const prediction = Pmpo.predict(df, this.params, this.predictionName);

    // Mark predictions with a color
    prediction.colors.setLinear(getOutputPalette(OPT_TYPE.MAX), {min: prediction.stats.min, max: prediction.stats.max});

    df.columns.remove(this.predictionName);
    df.columns.add(prediction);

    // Sort with respect to the scores
    this.view.grid.sort([prediction.name], [false]);

    this.updateGrid(descrStatsTable, selectedByPvalue, selectedByCorr);

    this.corrPlot.dataFrame = df;
    this.corrPlot.setOptions({
      xColumnNames: descriptorNames,
      yColumnNames: descriptorNames,
      showTitle: true,
      title: `${DESCR_TITLE} Correlations`,
    });

    this.barChart.dataFrame = weightsTable;

    this.hist.dataFrame = df;
    this.hist.setOptions({
      valueColumnName: prediction.name,
      splitColumnName: desirability.name,
      normalizeValues: false,
      splitStack: true,
      showSplitSelector: false,
      title: 'Distribution of pMPO scores',
      showTitle: true,
    });

    this.boxPlot.dataFrame = df;
    this.boxPlot.setOptions({
      valueColumnName: prediction.name,
      category1ColumnName: desirability.name,
      showPValue: false,
      showSizeSelector: false,
      showTitle: true,
      title: `pMPO by ${desirability.name} label`,
    });
  } // fit

  public runTrainingApp(): void {
    const dockMng = this.view.dockManager;

    // Inputs form
    dockMng.dock(this.getInputForm(), DG.DOCK_TYPE.LEFT, null, undefined, 0.25);

    // Dock viewers
    const gridNode = dockMng.findNode(this.view.grid.root);

    if (gridNode == null)
      throw new Error('Failed to train pMPO: missing a grid in the table view.');

    const statGridNode = dockMng.dock(this.statGrid, DG.DOCK_TYPE.RIGHT, gridNode, undefined, 0.5);
    const corrNode = dockMng.dock(this.corrPlot, DG.DOCK_TYPE.RIGHT, statGridNode, undefined, 0.25);
    dockMng.dock(this.barChart, DG.DOCK_TYPE.DOWN, corrNode, undefined, 0.5);
    dockMng.dock(this.hist, DG.DOCK_TYPE.DOWN, statGridNode, undefined, 0.5);
    dockMng.dock(this.boxPlot, DG.DOCK_TYPE.DOWN, gridNode, undefined, 0.5);
  } // runTrainingApp

  private getInputForm(): HTMLElement {
    const form = ui.form([]);
    form.append(ui.h2('Training data'));
    const numericColNames = this.numericCols.map((col) => col.name);

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

    const pInput = ui.input.float('p-value', {
      nullable: false,
      min: P_VAL_TRES_MIN,
      max: 1,
      step: 0.01,
      value: 0.05,
      tooltipText: 'Descriptors with p-values above this threshold are excluded.',
      onValueChanged: (value) => {
        if ((value != null) && (value >= P_VAL_TRES_MIN) && (value <= 1))
          runComputations();
      },
    });
    form.append(pInput.root);

    const rInput = ui.input.float('R²', {
      nullable: false,
      min: R2_MIN,
      value: 0.5,
      max: 1,
      step: 0.01,
      // eslint-disable-next-line max-len
      tooltipText: 'Descriptors with correlation coefficients above this threshold are considered highly correlated; the one with the lower p-value is retained.',
      onValueChanged: (value) => {
        if ((value != null) && (value >= R2_MIN) && (value <= 1))
          runComputations();
      },
    });
    form.append(rInput.root);

    const qInput = ui.input.float('q-cutoff', {
      nullable: false,
      min: Q_CUTOFF_MIN,
      value: 0.05,
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

  private getBoolCols(): DG.Column[] {
    const res: DG.Column[] = [];

    for (const col of this.table.columns) {
      if ((col.type === DG.COLUMN_TYPE.BOOL) && (col.stats.stdev > 0))
        res.push(col);
    }

    return res;
  } // getBoolCols

  private getValidNumericCols(): DG.Column[] {
    const res: DG.Column[] = [];

    for (const col of this.table.columns) {
      if ((col.isNumerical) && (col.stats.missingValueCount < 1) && (col.stats.stdev > 0))
        res.push(col);
    }

    return res;
  } // getValidNumericCols
}; // Pmpo
