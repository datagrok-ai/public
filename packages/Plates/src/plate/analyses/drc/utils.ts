import {Plate, IPlateWellFilter} from '../../plate';
import {FitSeries} from '@datagrok-libraries/statistics/src/fit/new-fit-API';

interface ISeriesData {
  x: number[];
  y: number[];
  meta: number[];
  outlier: boolean[];
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
