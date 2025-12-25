import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getDesiredTables, DescriptorStatistics, getDescriptorStatistics} from './utils';

export class Pmpo {
  static isApplicable(descriptors: DG.ColumnList, desirability: DG.Column): boolean {
    if (descriptors.length < 1)
      return false;

    for (const col of descriptors) {
      if (!col.isNumerical)
        return false;

      if (col.stats.missingValueCount > 0)
        return false;

      if (col.stats.stdev === 0)
        return false;
    }

    return (desirability.type === DG.COLUMN_TYPE.BOOL) && (desirability.stats.missingValueCount < 1);
  }

  constructor() {};

  public fit(df: DG.DataFrame, descriptors: DG.ColumnList, desirability: DG.Column): void {
    if (!Pmpo.isApplicable(descriptors, desirability))
      throw new Error('Failed to train pMPO model: the method is not applicable to the inputs');

    const descriptorNames = descriptors.names();
    const {desired, nonDesired} = getDesiredTables(df, desirability);

    grok.shell.addTableView(desired);
    grok.shell.addTableView(nonDesired);

    const descrStats = new Map<string, DescriptorStatistics>();
    descriptorNames.forEach((name) => {
      descrStats.set(name, getDescriptorStatistics(desired.col(name)!, nonDesired.col(name)!));
    });

    console.log(descrStats);
  } // fit
}; // Pmpo
