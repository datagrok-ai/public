import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getDesiredTables, getDescriptorStatistics, normalPdf, sigmoidS} from './stat-tools';
import {MIN_SAMPLES_COUNT, PMPO_NON_APPLICABLE, DescriptorStatistics, P_VAL_TRES_MIN, DESCR_TITLE,
  R2_MIN, Q_CUTOFF_MIN, WEIGHT_TABLE_TITLE, PmpoParams, SCORES_TITLE} from './pmpo-defs';
import {addSelectedDescriptorsCol, getDescriptorStatisticsGrid, getDescriptorStatisticsTable,
  getFilteredByPvalue, getFilteredByCorrelations, getModelParams, getWeightsTable} from './pmpo-utils';


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

  private params: Map<string, PmpoParams> | null = null;

  constructor() {};

  public fit(df: DG.DataFrame, descriptors: DG.ColumnList, desirability: DG.Column,
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
    const prediction = this.predict(df);
    df.columns.add(prediction);
    //grok.shell.addTableView(weightsTable);

    // Add viewers
    const grid = getDescriptorStatisticsGrid(descrStatsTable, selectedByPvalue, selectedByCorr);
    const view = grok.shell.tableView(df.name) ?? grok.shell.addTableView(df);
    const gridNode = view.dockManager.dock(grid, DG.DOCK_TYPE.RIGHT, null, undefined, 0.7);

    const corrPlot = DG.Viewer.correlationPlot(df, {
      xColumnNames: descriptorNames,
      yColumnNames: descriptorNames,
      showTitle: true,
      title: `${DESCR_TITLE} Correlations`,
    });
    const corrNode = view.dockManager.dock(corrPlot, DG.DOCK_TYPE.RIGHT, gridNode, undefined, 0.3);

    const bars = DG.Viewer.barChart(weightsTable, {
      valueAggrType: DG.AGG.AVG,
      showCategorySelector: false,
      showValueSelector: false,
      showTitle: true,
      title: WEIGHT_TABLE_TITLE,
      showStackSelector: false,
    });
    view.dockManager.dock(bars, DG.DOCK_TYPE.DOWN, gridNode, undefined, 0.5);
  } // fit

  public predict(df: DG.DataFrame): DG.Column {
    if (this.params == null)
      throw new Error('Filed to apply pMPO: non-trained model');

    const count = df.rowCount;
    const scores = new Float64Array(count).fill(0);
    let x = 0;

    this.params.forEach((param, name) => {
      const col = df.col(name);

      if (col == null)
        throw new Error(`'Filed to apply pMPO: inconsistent data, no column "${name}" in the table "${df.name}"`);

      const vals = col.getRawData();
      for (let i = 0; i < count; ++i) {
        x = vals[i];
        scores[i] += param.weight * normalPdf(x, param.desAvg, param.desStd) * sigmoidS(x, param.x0, param.b, param.c);
      }
    });

    return DG.Column.fromFloat64Array(SCORES_TITLE, scores);
  }
}; // Pmpo
