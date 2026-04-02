import {parsePzfxXml, pzfxTableToFitChartData} from '../../formats/pzfx/pzfx-parser';


/** Converts a PZFX XML string (single table or full document) to native IFitChartData JSON string.
 * Finds the first XY table in the parsed PZFX data and converts it. */
export function convertPzfxToJson(value: string): string {
  const tables = parsePzfxXml(value);
  const xyTable = tables.find((t) => t.tableType === 'XY');
  if (!xyTable)
    return JSON.stringify({series: []});
  const chartData = pzfxTableToFitChartData(xyTable);
  return JSON.stringify(chartData ?? {series: []});
}
