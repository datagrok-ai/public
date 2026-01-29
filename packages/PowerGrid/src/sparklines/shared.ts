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
  minValues?: Map<string, number>;
  maxValues?: Map<string, number>;
};

export function distance(p1: DG.Point, p2: DG.Point): number {
  return Math.sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

const FLOAT_NONE = 2.6789344063684636e-34;

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
    minValues,
    maxValues,
  } = settings;

  const toLogSafe = (v: number) => v > 0 ? Math.log(v) : FLOAT_NONE;

  const scaleValue = (v: number): number => logScale ? toLogSafe(v) : v;

  const resolveMinMax = (column: DG.Column): { min: number; max: number } => {
    const rawMin = minValues?.get(column.name);
    const rawMax = maxValues?.get(column.name);

    const min = rawMin != null && rawMin !== FLOAT_NONE ? rawMin : column.min;
    const max = rawMax != null && rawMax !== FLOAT_NONE ? rawMax : column.max;

    return {
      min: scaleValue(min),
      max: scaleValue(max),
    };
  };

  const normalize = (value: number, min: number, max: number): number =>
    max === min ? 0 : (value - min) / (max - min);

  const rowValues: number[] = [];
  const colMins: number[] = [];
  const colMaxs: number[] = [];

  for (const col of cols) {
    rowValues.push(scaleValue(col.getNumber(row)));

    const {min, max} = resolveMinMax(col);
    colMins.push(min);
    colMaxs.push(max);
  }

  const value = scaleValue(activeColumn.getNumber(row));
  let normalized: number;

  if (normalization === NormalizationType.Global || normalization === NormalizationType.Row) {
    const mins = normalization === NormalizationType.Global ? colMins : rowValues;
    const maxs = normalization === NormalizationType.Global ? colMaxs : rowValues;

    const globalMin = zeroScale ? 0 : Math.min(...mins);
    const globalMax = Math.max(...maxs);
    normalized = normalize(value, globalMin, globalMax);
  } else {
    const {min, max} = resolveMinMax(activeColumn);
    normalized = normalize(value, min, max);
  }

  return invertScale ? 1 - normalized : normalized;
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
    if (isSmartForm || !settings.minValues?.size || !settings.maxValues?.size)
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

  return [createColumnsInput(), createNormalizationInput(), createColorCodeInput()].filter(Boolean) as DG.InputBase[];
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
