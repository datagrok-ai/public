import {convertXMLToIFitChartData} from '../fit-parser';


/** Converts XML 3DX curve format to native IFitChartData JSON string.
 * Registered as a Datagrok function so it can be referenced by nqName. */
export function convertXmlCurveToJson(value: string): string {
  const chartData = convertXMLToIFitChartData(value);
  return JSON.stringify(chartData);
}
