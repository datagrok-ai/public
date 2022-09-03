import * as DG from 'datagrok-api/dg';

export function names(columns: Iterable<DG.Column>): string[] {
  return Array.from(columns).map((c: any) => c.name);
}

export interface SummarySettingsBase {
  columnNames: string[];
}

export function getSettingsBase(gc: DG.GridColumn): SummarySettingsBase {
  return gc.settings ??= {
    columnNames: names(gc.grid.dataFrame.columns.numerical),
  };
}

export enum SparklineType {
  BarChart = 'barchart',
  PieChart = 'piechart',
  Radar = 'radar',
  Sparkline = 'sparkline'
}

export const sparklineTypes: string[] = [
  SparklineType.BarChart,
  SparklineType.PieChart,
  SparklineType.Radar,
  SparklineType.Sparkline,
];

// function getDataColumns<Type extends SummarySettingsBase>(
//   gc: DG.GridColumn, getSettings: getSettingsFunc<Type>,
// ): DG.Column[] {
//   return gc.grid.dataFrame.columns.byNames(getSettings(gc).columnNames);
// }
