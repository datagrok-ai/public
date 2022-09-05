import * as DG from 'datagrok-api/dg';
import wu from 'wu';
import * as ui from 'datagrok-api/ui';

type getSettingsFunc<Type extends SummarySettingsBase> = (gs: DG.GridColumn) => Type;

export interface CustomMouseEvent extends MouseEvent {
  layerX: number;
  layerY: number;
}

export function names(columns: Iterable<DG.Column>): string[] {
  return Array.from(columns).map((c: DG.Column) => c.name);
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

export function distance(p1: DG.Point, p2: DG.Point): number {
  return Math.sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

export class Hit {
  activeColumn: number = -1;
  cols: DG.Column[] = [];
  row: number = -1;
  isHit: boolean = false;
}

export function createTooltip(cols: DG.Column[], activeColumn: number, row: number): HTMLDivElement[] {
  let arr: HTMLDivElement[] = [];
  for (let i = 0; i < cols.length; i++) {
    arr.push(ui.divH([ui.divText(`${cols[i].name}:`, {
          style: {
            margin: '0 10px 0 0',
            fontWeight: (activeColumn == i) ? 'bold' : 'normal',
          }
        }), ui.divText(`${Math.floor(cols[i].get(row) * 100) / 100}`, {
          style: {
            fontWeight: (activeColumn == i) ? 'bold' : 'normal',
          }
        })]
      )
    );
  }
  return arr;
}

interface staticSettingsBarChart {
  minH: number;
}

export const renderSettingsBarChart: staticSettingsBarChart = {
  minH: 0.05,
};

interface staticSettingsPieChart {
  minRadius: number;
}

export const renderSettingsPieChart: staticSettingsPieChart = {
  minRadius: 10,
};

interface staticSettingsSparkline {
  minDistance: number;
}

export const renderSettingsSparkline: staticSettingsSparkline = {
  minDistance: 5,
};

// function getDataColumns<Type extends SummarySettingsBase>(
//   gc: DG.GridColumn, getSettings: getSettingsFunc<Type>,
// ): DG.Column[] {
//   return gc.grid.dataFrame.columns.byNames(getSettings(gc).columnNames);
// }
