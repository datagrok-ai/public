import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getDesiredTables, getDescriptorStatistics} from './stat-tools';
import {MIN_SAMPLES_COUNT, PMPO_NON_APPLICABLE, DescriptorStatistics, P_VAL_TRES_MIN, DESCR_TITLE} from './pmpo-defs';
import {addSelectedDescriptorsCol, getDescriptorStatisticsGrid, getDescriptorStatisticsTable,
  getFilteredByPvalue} from './pmpo-utils';


export class Pmpo {
  static isApplicable(descriptors: DG.ColumnList, desirability: DG.Column, pValThresh: number,
    toShowWarning: boolean = false): boolean {
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

  constructor() {};

  public fit(df: DG.DataFrame, descriptors: DG.ColumnList, desirability: DG.Column, pValTresh: number): void {
    if (!Pmpo.isApplicable(descriptors, desirability, pValTresh))
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
    const selected = getFilteredByPvalue(descrStatsTable, pValTresh);
    addSelectedDescriptorsCol(descrStatsTable, selected);

    // Add viewers
    const grid = getDescriptorStatisticsGrid(descrStatsTable);
    const view = grok.shell.tableView(df.name) ?? grok.shell.addTableView(df);
    const node = view.dockManager.dock(grid, DG.DOCK_TYPE.RIGHT, null, undefined, 0.7);
    const corrPlot = DG.Viewer.correlationPlot(df, {
      xColumnNames: descriptorNames,
      yColumnNames: descriptorNames,
      showTitle: true,
      title: `${DESCR_TITLE} Correlations`,
    });
    view.dockManager.dock(corrPlot, DG.DOCK_TYPE.RIGHT, node, undefined, 0.3);
  } // fit
}; // Pmpo
