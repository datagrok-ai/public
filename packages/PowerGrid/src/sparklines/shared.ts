import * as DG from 'datagrok-api/dg';
import wu from 'wu';
import * as ui from 'datagrok-api/ui';

type getSettingsFunc<Type extends SummarySettingsBase> = (gs: DG.GridColumn) => Type;


export function names(columns: Iterable<DG.Column>): string[] {
  return Array.from(columns).map((c: DG.Column) => c.name);
}

export enum SummaryColumnColoringType {
  Off = 'Off',
  Bins = 'Bins',
  Values = 'Values'
}

export interface SummarySettingsBase {
  columnNames: string[];
  colorCode: SummaryColumnColoringType;
}


/// Utility method for old summary columns format support
export function isSummarySettingsBase(obj: any): obj is SummarySettingsBase {
  return (obj as SummarySettingsBase).columnNames !== undefined;
}

export function getSettingsBase<Type extends SummarySettingsBase>(gc: DG.GridColumn,
  sparklineType: SparklineType): Type {
  return isSummarySettingsBase(gc.settings) ? gc.settings :
    gc.settings[sparklineType] ??= {
      columnNames: names(wu(gc.grid.dataFrame.columns.numerical)
        .filter((c: DG.Column) => c.type != DG.TYPE.DATE_TIME)),
    };
}

export enum SparklineType {
  BarChart = 'barchart',
  PieChart = 'piechart',
  Radar = 'radar',
  Sparkline = 'sparkline',
  Form = 'smartform'
}

export const sparklineTypes: string[] = [
  SparklineType.BarChart,
  SparklineType.PieChart,
  SparklineType.Radar,
  SparklineType.Sparkline,
  SparklineType.Form
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

export function createTooltip(cols: DG.Column[], activeColumn: number, row: number): HTMLDivElement {
  const keysDiv = ui.divV([], { style: { marginRight: '10px', fontWeight: 'bold', textAlign: 'right' } });
  const valuesDiv = ui.divV([], { style: { fontWeight: 'normal' } });

  for (let i = 0; i < cols.length; i++) {
    if (cols[i] === null) continue;

    const isActive = activeColumn === i;
    keysDiv.appendChild(ui.divText(`${cols[i].name}`, { 
      style: { fontWeight: isActive ? 'bold' : 'normal' } 
    }));
    valuesDiv.appendChild(ui.divText(`${Math.floor(cols[i].getNumber(row) * 100) / 100}`, { 
      style: { fontWeight: isActive ? 'bold' : 'normal' } 
    }));
  }

  return ui.divH([keysDiv, valuesDiv], { style: { display: 'flex' } });
}

export function createBaseInputs(gridColumn: DG.GridColumn, settings: SummarySettingsBase): DG.InputBase[] {
  const columnNames = settings?.columnNames ?? names(gridColumn.grid.dataFrame.columns.numerical);
  return [
    ui.input.columns('Columns', {
      value: gridColumn.grid.dataFrame.columns.byNames(columnNames),
      table: gridColumn.grid.dataFrame,
      onValueChanged: (value) => {
        settings.columnNames = names(value);
        gridColumn.grid.invalidate();
      },
      available: names(gridColumn.grid.dataFrame.columns.numerical)
    }),
    ui.input.choice<SummaryColumnColoringType>('Color Code', {
      value: settings.colorCode,
      items: [SummaryColumnColoringType.Off, SummaryColumnColoringType.Bins, SummaryColumnColoringType.Values],
      onValueChanged: (value) => {
        settings.colorCode = value;
        gridColumn.grid.invalidate();
      },
      tooltipText: 'Activates color rendering',
      nullable: false
    }),
  ];
}

export function getRenderColor(settings: SummarySettingsBase, baseColor: number,
   options: {column: DG.Column, colIdx: number, rowIdx: number}): number {
  return settings.colorCode === SummaryColumnColoringType.Off ? baseColor : settings.colorCode === SummaryColumnColoringType.Bins ?
    DG.Color.getCategoricalColor(options.colIdx) : options.column.meta.colors.getType() === DG.COLOR_CODING_TYPE.OFF ?
    DG.Color.getRowColor(options.column, options.rowIdx) : options.column.meta.colors.getColor(options.rowIdx);
}


// function getDataColumns<Type extends SummarySettingsBase>(
//   gc: DG.GridColumn, getSettings: getSettingsFunc<Type>,
// ): DG.Column[] {
//   return gc.grid.dataFrame.columns.byNames(getSettings(gc).columnNames);
// }
