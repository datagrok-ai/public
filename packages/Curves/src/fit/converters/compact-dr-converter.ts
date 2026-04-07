import {
  IFitChartData,
  IFitPoint,
  IFitSeries,
  FIT_FUNCTION_4PL_DOSE_RESPONSE,
} from '@datagrok-libraries/statistics/src/fit/fit-curve';


/** Compact dose-response JSON format:
 * {
 *   "rdt": "2025-05-31 00:00:00.0",   // run date
 *   "mnr": -4.28,                       // min response
 *   "mxr": 109.134,                     // max response
 *   "yu": "%",                          // y unit
 *   "hs": -1.5498,                      // hill slope
 *   "poi": 0.00647580,                  // point of inflection (e.g. IC50)
 *   "xu": "uM",                         // x unit
 *   "ctp": "dose-response",             // curve type
 *   "p": [[0, 3], [0.3, 4], [0.5, 7, 1]]  // points [x, y, outlier?]
 * }
 */
interface ICompactDoseResponse {
  rdt?: string;
  mnr?: number;
  mxr?: number;
  yu?: string;
  hs?: number;
  poi?: number;
  xu?: string;
  ctp?: string;
  p?: number[][];
}

/** Converts compact dose-response JSON to native IFitChartData JSON string. */
export function convertCompactDrToJson(value: string): string {
  const data: ICompactDoseResponse = JSON.parse(value);
  const chartData = convertCompactDrToChartData(data);
  return JSON.stringify(chartData);
}

function convertCompactDrToChartData(data: ICompactDoseResponse): IFitChartData {
  const points: IFitPoint[] = (data.p ?? []).map((arr) => ({
    x: arr[0],
    y: arr[1],
    outlier: arr.length > 2 ? arr[2] === 1 : false,
  }));

  const hasParams = data.mxr !== undefined && data.hs !== undefined &&
    data.poi !== undefined && data.mnr !== undefined;

  const series: IFitSeries = {
    points: points,
    fitFunction: FIT_FUNCTION_4PL_DOSE_RESPONSE,
    showFitLine: true,
    showPoints: 'points',
    clickToToggle: true,
    name: data.rdt ? `${data.rdt}` : 'Dose-Response',
    droplines: ['IC50']
  };

  // Map to sigmoid parameters: [max, slope, IC50, min]
  if (hasParams)
    series.parameters = [data.mxr!, data.hs!, data.poi!, data.mnr!];

  return {
    chartOptions: {
      xAxisName: data.xu,
      yAxisName: data.yu,
      logX: true,
    },
    series: [series],
  };
}
