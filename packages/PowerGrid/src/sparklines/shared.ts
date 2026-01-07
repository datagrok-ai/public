/* eslint-disable max-len */
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
  Values = 'Values',
  Auto = 'Auto',
}

export enum NormalizationType {
  Row = 'Row',
  Column = 'Column',
  Global = 'Global'
}

export type ScaleSettings = {
  normalization: NormalizationType;
  zeroScale?: boolean;
  invertScale?: boolean;
  logScale?: boolean;
};

export interface SummarySettingsBase {
  columnNames: string[];
  logColumnNames: string[];
  invertColumnNames: string[];
  minValues: Map<string, number>;
  maxValues: Map<string, number>;
  colorCode: SummaryColumnColoringType;
  normalization: NormalizationType;
  mode: MinMaxMode;
}

export enum MinMaxMode {
  Auto = 'Auto',
  Manual = 'Manual',
}


/// Utility method for old summary columns format support
export function isSummarySettingsBase(obj: any): obj is SummarySettingsBase {
  return (obj as SummarySettingsBase).columnNames !== undefined;
}

export function getSettingsBase<T extends SummarySettingsBase>(
  gc: DG.GridColumn,
  sparklineType: SparklineType
): T {
  const settings = isSummarySettingsBase(gc.settings) ?
    (gc.settings as unknown as T) :
    (gc.settings[sparklineType] ??= {
      columnNames: names(
        wu(gc.grid.dataFrame.columns.numerical)
          .filter((c: DG.Column) => c.type !== DG.TYPE.DATE_TIME)
      ).slice(0, 10),
    } as unknown as T);

  if (!settings.minValues || !settings.maxValues) {
    settings.minValues = new Map<string, number>();
    settings.maxValues = new Map<string, number>();

    for (const col of gc.grid.dataFrame.columns) {
      settings.minValues.set(col.name, col.min);
      settings.maxValues.set(col.name, col.max);
    }
  }

  return settings;
}

export enum SparklineType {
  BarChart = 'barchart',
  PieChart = 'piechart',
  Radar = 'radar',
  Sparkline = 'sparkline',
  Form = 'smartform',
  Tags = 'tags',
}

export const sparklineTypes: string[] = [
  SparklineType.BarChart,
  SparklineType.PieChart,
  SparklineType.Radar,
  SparklineType.Sparkline,
  SparklineType.Form
];

type AxisScaleSettings = ScaleSettings & {
  mode?: MinMaxMode;
  minValues?: Map<string, number>;
  maxValues?: Map<string, number>;
};

export function distance(p1: DG.Point, p2: DG.Point): number {
  return Math.sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

const FLOAT_NONE = 2.6789344063684636e-34;

function toLogSafe(x: number): number {
  return x > 0 ? Math.log(x) : FLOAT_NONE;
}

export function getScaledNumber(
  cols: DG.Column[],
  row: number,
  activeColumn: DG.Column,
  settings: AxisScaleSettings
): number {
  const {
    normalization,
    zeroScale = false,
    invertScale = false,
    logScale = false,
    mode = MinMaxMode.Auto,
    minValues,
    maxValues,
  } = settings;

  const colNumbers: number[] = [];
  const colMins: number[] = [];
  const colMaxs: number[] = [];

  for (const c of cols) {
    const num = c?.getNumber(row);
    colNumbers.push(logScale ? toLogSafe(num) : num);

    const min = mode === MinMaxMode.Manual ? minValues?.get(c.name) : c.min;
    const max = mode === MinMaxMode.Manual ? maxValues?.get(c.name) : c.max;
    if (min != null)
      colMins.push(logScale ? toLogSafe(min) : min);
    if (max != null)
      colMaxs.push(logScale ? toLogSafe(max) : max);
  }

  const value = logScale ? toLogSafe(activeColumn.getNumber(row)) : activeColumn.getNumber(row);
  let normalized = 0;

  if (normalization === NormalizationType.Global || normalization === NormalizationType.Row) {
    const values = normalization === NormalizationType.Global ? colMins : colNumbers;
    const ranges = normalization === NormalizationType.Global ? colMaxs : colNumbers;

    const gmin = zeroScale ? 0 : Math.min(...values);
    const gmax = Math.max(...ranges);

    normalized = gmax === gmin ? 0 : (value - gmin) / (gmax - gmin);
  } else {
    const min = mode === MinMaxMode.Manual ? minValues?.get(activeColumn.name) ?? activeColumn.min : activeColumn.min;
    const max = mode === MinMaxMode.Manual ? maxValues?.get(activeColumn.name) ?? activeColumn.max : activeColumn.max;

    const minVal = logScale ? toLogSafe(min) : min;
    const maxVal = logScale ? toLogSafe(max) : max;

    normalized = maxVal === minVal ? 0 : (value - minVal) / (maxVal - minVal);
  }

  if (invertScale)
    normalized = 1 - normalized;

  return normalized;
}

export function getSparklinesContextPanel(gridCell: DG.GridCell, colNames: string[]): HTMLDivElement {
  const df = gridCell.grid.dataFrame;
  const row = gridCell.cell.row.idx;
  const cols = df.columns.byNames(colNames).filter((c) => c !== null);
  const columnName = ui.div(gridCell.gridColumn.name, {style: {textAlign: 'center', marginBottom: '10px'}});
  const values = ui.divV(cols.map((col) => ui.divText(col.name + ': ' + col.getString(row))), {style: {marginTop: '20px'}});
  return ui.div([columnName, DG.GridCellWidget.fromGridCell(gridCell).root, values]);
}

export class Hit {
  activeColumn: number = -1;
  cols: DG.Column[] = [];
  row: number = -1;
  isHit: boolean = false;
}

export function createTooltip(cols: DG.Column[], activeColumn: number, row: number): HTMLDivElement {
  const keysDiv = ui.divV([], {style: {marginRight: '10px', fontWeight: 'bold', textAlign: 'right'}});
  const valuesDiv = ui.divV([], {style: {fontWeight: 'normal'}});

  for (let i = 0; i < cols.length; i++) {
    if (cols[i] === null) continue;

    const isActive = activeColumn === i;
    keysDiv.appendChild(ui.divText(`${cols[i].name}`, {
      style: {fontWeight: isActive ? 'bold' : 'normal'}
    }));
    valuesDiv.appendChild(ui.divText(`${Math.floor(cols[i].getNumber(row) * 100) / 100}`, {
      style: {fontWeight: isActive ? 'bold' : 'normal'}
    }));
  }

  return ui.divH([keysDiv, valuesDiv], {style: {display: 'flex'}});
}

export function createBaseInputs(gridColumn: DG.GridColumn, settings: SummarySettingsBase, isSmartForm: boolean = false): DG.InputBase[] {
  const df = gridColumn.grid.dataFrame;
  const invalidate = () => gridColumn.grid.invalidate();

  function createNormalizationInput(): DG.InputBase | null {
    if (isSmartForm) return null;

    return ui.input.choice<NormalizationType>('Normalization', {
      value: settings.normalization,
      items: [NormalizationType.Row, NormalizationType.Column, NormalizationType.Global],
      onValueChanged: (value) => {
        settings.normalization = value;
        invalidate();
      },
      tooltipText: 'Defines how values are scaled:<br>' +
        '- ROW: Scales each row individually (row minimum to row maximum). Use for comparing values within a row.<br>' +
        '- COLUMN: Scales each column individually (column minimum to column maximum).' +
        'Use when columns have different units.<br>' +
        '- GLOBAL: Applies a single scale across all values.' +
        'Use for comparing values across columns with the same units (e.g., tracking changes over time).',
      nullable: false
    });
  }

  function getMinMaxProperties(): DG.Property[] | null {
    if (isSmartForm || settings.mode === MinMaxMode.Auto || !settings.minValues?.size || !settings.maxValues?.size)
      return null;

    return [
      DG.Property.create('min', DG.TYPE.FLOAT,
        (col: string) => settings.minValues.get(col),
        (col: string, value: number) => settings.minValues.set(col, value)
      ),
      DG.Property.create('max', DG.TYPE.FLOAT,
        (col: string) => settings.maxValues.get(col),
        (col: string, value: number) => settings.maxValues.set(col, value)
      ),
    ];
  }

  function createModeInput(): DG.InputBase | null {
    if (isSmartForm)
      return null;

    return ui.input.choice<MinMaxMode>('Mode', {
      value: settings.mode ?? MinMaxMode.Auto,
      items: [MinMaxMode.Auto, MinMaxMode.Manual],
      onValueChanged: (value) => {
        settings.mode = value;
        invalidate();
      },
      tooltipText: 'Auto: system calculates min/max.<br>' +
      'Manual: user sets min/max explicitly.',
    });
  }

  function getAdditionalColumns() {
    if (isSmartForm)
      return null;

    return {
      additionalColumns: {
        'log': df.columns.byNames(settings.logColumnNames ?? []),
        'invert': df.columns.byNames(settings.invertColumnNames ?? []),
      },
      onAdditionalColumnsChanged: (values: { [key: string]: DG.Column[] }) => {
        settings.logColumnNames = names(values['log'] ?? []);
        settings.invertColumnNames = names(values['invert'] ?? []);
        invalidate();
      },
    };
  }

  function createColumnsInput(): DG.InputBase {
    const columnNames = settings?.columnNames ?? names(df.columns.numerical);
    const options: any = {
      value: df.columns.byNames(columnNames),
      table: df,
      available: isSmartForm ? names(df.columns) : names(df.columns.numerical),
      onValueChanged: (value: DG.Column[]) => {
        settings.columnNames = names(value);
        invalidate();
      },
    };

    const minMax = getMinMaxProperties();
    if (minMax)
      options.additionalColumnProperties = minMax;

    const additionalCols = getAdditionalColumns();
    if (additionalCols) {
      options.additionalColumns = additionalCols.additionalColumns;
      options.onAdditionalColumnsChanged = additionalCols.onAdditionalColumnsChanged;
    }
    return ui.input.columns('Columns', options);
  }

  function createColorCodeInput(): DG.InputBase {
    return ui.input.choice<SummaryColumnColoringType>('Color Code', {
      value: settings.colorCode,
      items: [
        SummaryColumnColoringType.Auto,
        SummaryColumnColoringType.Bins,
        SummaryColumnColoringType.Values,
        SummaryColumnColoringType.Off,
      ],
      onValueChanged: (value) => {
        settings.colorCode = value;
        invalidate();
      },
      tooltipText: 'Activates color rendering',
      nullable: false,
    });
  }

  return [createColumnsInput(), createModeInput(), createNormalizationInput(), createColorCodeInput()].filter(Boolean) as DG.InputBase[];
}

export function getRenderColor(settings: SummarySettingsBase, baseColor: number,
  options: {column: DG.Column, colIdx: number, rowIdx: number}): number {
  return settings.colorCode === SummaryColumnColoringType.Off ? baseColor : // off, base color
    settings.colorCode === SummaryColumnColoringType.Bins ?
      DG.Color.getCategoricalColor(options.colIdx) : settings.colorCode === SummaryColumnColoringType.Values ? (options.column.meta.colors.getType() === DG.COLOR_CODING_TYPE.OFF ?
        DG.Color.getRowColor(options.column, options.rowIdx) : options.column.meta.colors.getColor(options.rowIdx)) : (options.column.meta.colors.getType() === DG.COLOR_CODING_TYPE.OFF ? baseColor : options.column.meta.colors.getColor(options.rowIdx));
}


// function getDataColumns<Type extends SummarySettingsBase>(
//   gc: DG.GridColumn, getSettings: getSettingsFunc<Type>,
// ): DG.Column[] {
//   return gc.grid.dataFrame.columns.byNames(getSettings(gc).columnNames);
// }
