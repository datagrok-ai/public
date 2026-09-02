/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import '../../css/powergrid.css';

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
  minValues: Record<string, number>;
  maxValues: Record<string, number>;
  colorCode: SummaryColumnColoringType;
  normalization: NormalizationType;
  useFilteredData: boolean;
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
      columnNames: names(gc.grid.dataFrame.columns.numericalNoDateTime).slice(0, 10),
    } as unknown as T);

  if (!settings.minValues || !settings.maxValues) {
    settings.minValues = {};
    settings.maxValues = {};

    for (const col of gc.grid.dataFrame.columns) {
      settings.minValues[col.name] = col.min;
      settings.maxValues[col.name] = col.max;
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
  minValues?: Record<string, number>;
  maxValues?: Record<string, number>;
  useFilteredData?: boolean;
};

export function scaleSettings(
  settings: SummarySettingsBase,
  column: DG.Column,
  overrides?: Partial<ScaleSettings>
): AxisScaleSettings {
  return {
    normalization: settings.normalization,
    invertScale: settings.invertColumnNames?.includes(column.name),
    logScale: settings.logColumnNames?.includes(column.name),
    minValues: settings.minValues,
    maxValues: settings.maxValues,
    useFilteredData: settings.useFilteredData,
    ...overrides,
  };
}

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
    useFilteredData = false,
  } = settings;

  const toLogSafe = (v: number) => v > 0 ? Math.log(v) : FLOAT_NONE;

  const scaleValue = (v: number): number => logScale ? toLogSafe(v) : v;

  const resolveMinMax = (column: DG.Column): { min: number; max: number } => {
    const rawMin = minValues?.[column.name];
    const rawMax = maxValues?.[column.name];

    let colMin: number;
    let colMax: number;

    if (useFilteredData) {
      const stats = DG.Stats.fromColumn(column, column.dataFrame.filter);
      colMin = stats.min;
      colMax = stats.max;
    } else {
      colMin = column.min;
      colMax = column.max;
    }

    const min = rawMin != null && rawMin !== FLOAT_NONE ? rawMin : colMin;
    const max = rawMax != null && rawMax !== FLOAT_NONE ? rawMax : colMax;

    return {min: scaleValue(min), max: scaleValue(max)};
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

export interface ColumnGroup {
  name: string;
  cols: DG.Column[];
  color?: string;
}

function appendGroups(root: HTMLElement, groups: ColumnGroup[], value: (col: DG.Column) => string, active?: DG.Column): void {
  for (const group of groups) {
    const name = ui.divText(group.name, 'power-grid-sparkline-group-name');
    if (group.color)
      name.prepend(ui.span([], {classes: 'power-grid-sparkline-swatch', style: {background: group.color}}));
    const host = group.name ? ui.div([name], 'power-grid-sparkline-group') : root;
    for (const col of group.cols) {
      const cls = col === active ? ' power-grid-sparkline-active' : '';
      host.append(ui.divText(col.name, 'power-grid-sparkline-name' + cls), ui.divText(value(col), 'power-grid-sparkline-value' + cls));
    }
    if (host !== root)
      root.append(host);
  }
}

export function getSparklinesContextPanel(gridCell: DG.GridCell, colNames: string[], groups?: ColumnGroup[]): HTMLDivElement {
  const row = gridCell.cell.row.idx;
  const columnName = ui.div(gridCell.gridColumn.name, 'power-grid-sparkline-panel-title');
  const values = ui.div([], 'power-grid-sparkline-panel-values');
  if (groups) {
    values.classList.add('power-grid-sparkline-grid');
    appendGroups(values, groups, (col) => col.getString(row));
  } else {
    const cols = gridCell.grid.dataFrame.columns.byNames(colNames).filter((c) => c != null);
    values.append(...cols.map((col) => ui.divText(col.name + ': ' + col.getString(row))));
  }
  return ui.div([columnName, DG.GridCellWidget.fromGridCell(gridCell).root, values]);
}

export class Hit {
  activeColumn: number = -1;
  cols: DG.Column[] = [];
  row: number = -1;
  isHit: boolean = false;
}

export function createTooltip(cols: DG.Column[], activeColumn: number, row: number, groups?: ColumnGroup[]): HTMLDivElement {
  const root = ui.div([], 'power-grid-sparkline-grid power-grid-sparkline-tooltip');
  appendGroups(root, groups ?? [{name: '', cols}],
    (col) => col.isNumerical ? `${Math.floor(col.getNumber(row) * 100) / 100}` : col.getString(row), cols[activeColumn]);
  return root;
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
    if (isSmartForm || !Object.keys(settings.minValues ?? {}).length || !Object.keys(settings.maxValues ?? {}).length)
      return null;

    return [
      DG.Property.create('min', DG.TYPE.FLOAT,
        (col: string) => settings.minValues[col],
        (col: string, value: number) => settings.minValues[col] = value
      ),
      DG.Property.create('max', DG.TYPE.FLOAT,
        (col: string) => settings.maxValues[col],
        (col: string, value: number) => settings.maxValues[col] = value
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
    const columnNames = settings?.columnNames ?? names(df.columns.numericalNoDateTime);
    const options: any = {
      value: df.columns.byNames(columnNames),
      table: df,
      available: isSmartForm ? names(df.columns) : names(df.columns.numericalNoDateTime),
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

  function createFilteredDataInput(): DG.InputBase | null {
    if (isSmartForm) return null;
    return ui.input.bool('Use Filtered Data', {
      value: settings.useFilteredData ?? false,
      onValueChanged: (value) => {
        settings.useFilteredData = value;
        invalidate();
      },
      tooltipText: 'When enabled, normalization uses min/max from filtered rows only',
    });
  }

  return [createColumnsInput(), createNormalizationInput(), createFilteredDataInput(), createColorCodeInput()].filter(Boolean) as DG.InputBase[];
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
