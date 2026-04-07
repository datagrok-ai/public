import {Plate, IPlateWellFilter} from '../../plate';
import {FitSeries} from '@datagrok-libraries/statistics/src/fit/new-fit-API';

interface ISeriesData {
  x: number[];
  y: number[];
  meta: number[];
  outlier: boolean[];
}

/**
 * Auto-detects DRC column mappings from plate column names using heuristics.
 * Matches are case-insensitive substring matches, checked in priority order.
 */
export function autoDetectDrcMappings(plate: Plate): Map<string, string> {
  const mappings = new Map<string, string>();
  const columnNames = plate.data.columns.names();

  const detect = (target: string, candidates: string[]) => {
    const lower = columnNames.map((c) => c.toLowerCase());
    for (const cand of candidates) {
      const idx = lower.findIndex((name) => name.includes(cand));
      if (idx !== -1) {
        mappings.set(target, columnNames[idx]);
        break;
      }
    }
  };

  detect('Activity', ['activity', 'response', 'readout', 'value', 'signal', 'raw data', 'raw_data']);
  detect('Concentration', ['concentration', 'conc', 'dose']);
  detect('SampleID', ['sample', 'compound', 'layout', 'plate layout', 'plate_layout', 'sampleid', 'sample_id']);

  return mappings;
}

export function getDoseResponseSeries(plate: Plate, options?: IPlateWellFilter & {
  concentration?: string;
  value?: string;
  groupBy?: string;
}): Record<string, FitSeries> {
  const valueOptions = {
    includeEmpty: options?.includeEmpty ?? false,
    exclude: options?.exclude ?? {'role': ['High Control', 'Low Control']}
  };
  const concKey = options?.concentration ?? 'concentration';
  const valueKey = options?.value ?? 'value';
  const values = plate.values([concKey, valueKey, ...(options?.groupBy ? [options.groupBy] : [])], valueOptions);

  const series: Record<string, ISeriesData & Record<string, any>> = {};
  for (const v of values) {
    const group = options?.groupBy ? v[options.groupBy] : '0';
    if (!series[group])
      series[group] = {x: [], y: [], meta: [], outlier: []};
    series[group].x.push(v[concKey]);
    series[group].y.push(v[valueKey]);
    series[group].meta.push(v.innerDfRow);

    const isOutlier = plate._isOutlier(v.innerDfRow);
    series[group].outlier.push(isOutlier);
  }

  return Object.fromEntries(Object.entries(series).map(([k, v]) => {
    const fitSeries = new FitSeries(v.x.map((_, i) => ({
      x: v.x[i],
      y: v.y[i],
      outlier: v.outlier[i],
      meta: v.meta[i]
    })).sort((a, b) => a.x - b.x));
    fitSeries.name = k;
    return [k, fitSeries];
  }));
}
