import * as DG from 'datagrok-api/dg';
import wu from "wu";

type getSettingsFunc<Type extends SummarySettingsBase> = (gs: DG.GridColumn) => Type;

export function names(columns: Iterable<DG.Column>): string[] {
  return Array.from(columns).map((c: any) => c.name);
}

export interface SummarySettingsBase {
  columnNames: string[];
}

export function getSettingsBase<Type extends SummarySettingsBase>(gc: DG.GridColumn): Type {
  return gc.settings ??= {
    columnNames: names(wu(gc.grid.dataFrame.columns.numerical)
      .filter((c: DG.Column) => c.type != DG.TYPE.DATE_TIME)),
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
