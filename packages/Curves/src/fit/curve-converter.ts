import * as DG from 'datagrok-api/dg';

import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';


/** Converter function type: takes a raw cell string value and returns native IFitChartData JSON string. */
export type CurveConverterFn = (value: string) => string;

/** Registry of converter functions, keyed by format name (e.g. '3dx', 'compact-dr', 'pzfx').
 * Local converters are registered at init time.
 * External converters are discovered via initExternalConverters(). */
const converterRegistry = new Map<string, CurveConverterFn>();

/** Value-based LRU cache for converted curve data.
 * Keys include format name + raw cell value; values are the converted JSON strings. */
const convertedValuesCache = new DG.LruCache<string, string>(2000);


/** Registers a curve converter function under a format name.
 * Local converters call this at init time. */
export function registerCurveConverter(formatName: string, fn: CurveConverterFn): void {
  converterRegistry.set(formatName, fn);
}

/** Discovers and registers curve converter functions from external packages.
 * Finds DG functions that have a meta.curveFormat annotation,
 * skipping any format names that are already registered (i.e. local converters). */
export async function initExternalConverters(): Promise<void> {
  const funcs = DG.Func.find({meta: {'role': 'curveConverter'}});
  for (const func of funcs) {
    const formatName = func.options['curveFormat'];
    if (formatName && !converterRegistry.has(formatName)) {
      const converterFn = await func.apply() as CurveConverterFn;
      registerCurveConverter(formatName, converterFn);
    }
  }
}

/** Returns the curve format name for a column, or null if native format. */
export function getCurveFormat(column: DG.Column): string | null {
  const tag = column.getTag(FitConstants.TAG_CURVE_FORMAT);
  return tag && tag !== '' ? tag : null;
}

/** Converts a raw cell value to native IFitChartData JSON string using the column's converter.
 * Returns the original value if no converter is set (native format).
 * Results are cached by value. */
export function convertCellValue(cellValue: string, column: DG.Column): string {
  const formatName = getCurveFormat(column);
  if (!formatName)
    return cellValue;

  const cacheKey = `${formatName}||${cellValue}`;
  return convertedValuesCache.getOrCreate(cacheKey, () => {
    const fn = converterRegistry.get(formatName);
    if (!fn) {
      console.warn(`Curve converter for format '${formatName}' not found in registry.`);
      return cellValue;
    }
    return fn(cellValue);
  });
}

/** Parses a cell value into IFitChartData, applying converter if needed.
 * This is the single entry point that replaces all hardcoded format checks. */
export function parseCellValue(cellValue: string, column?: DG.Column | null): IFitChartData {
  if (!cellValue || cellValue === '')
    return {series: []};

  const jsonStr = column ? convertCellValue(cellValue, column) : cellValue;
  return JSON.parse(jsonStr ?? '{}') ?? {};
}

/** Returns true if the column uses a non-native format (has a converter). */
export function hasConverter(column: DG.Column): boolean {
  return getCurveFormat(column) !== null;
}

/** Returns true if the column uses native JSON format (no converter needed). */
export function isNativeFormat(column: DG.Column): boolean {
  return !hasConverter(column);
}
