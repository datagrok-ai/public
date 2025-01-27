import { BaselineEndpointRenderLines, HysLawRenderLines, KaplanMeierRenderLines } from "./event-handlers";
import * as DG from "datagrok-api/dg";
import * as grok from 'datagrok-api/grok';

export function createHysLawScatterPlot(dataframe, xColumn, yColumn, colorColumn, APColumnName = '') {
  const sp = DG.Viewer.scatterPlot(dataframe, {
    x: xColumn,
    y: yColumn,
    color: dataframe.col(colorColumn) ? colorColumn : '',
    legendPosition: DG.FlexPosition.Right,
  });

  sp.onAfterDrawScene.subscribe(_ => HysLawRenderLines(sp, 3, 1, 2, APColumnName));
  return sp;
}


export function createBaselineEndpointScatterPlot(dataframe, xColumn, yColumn, colorColumn, lowLimit, highLimit) {
  const sp = DG.Viewer.scatterPlot(dataframe, {
    x: xColumn,
    y: yColumn,
    color: dataframe.col(colorColumn) ? colorColumn : '',
    legendPosition: DG.FlexPosition.Right,
  });

  if (!lowLimit || !highLimit) {
    const message = `information on limits hasn't been provided for: ${!lowLimit ? 'lower limit': ''}
      ${!highLimit ? !lowLimit ? 'and high limit' : 'high limit' : ''}`;
    grok.shell.warning(message);
    return;
  }

  sp.onAfterDrawScene.subscribe(_ => BaselineEndpointRenderLines(sp, lowLimit, highLimit));
  return sp;
}


export function createKaplanMeierScatterPlot(dataframe, subjIDColumn, xColumn, yColumn, groupColumn) {
  const sp = DG.Viewer.scatterPlot(dataframe, {
    x: xColumn,
    y: yColumn,
    markerType: 'dot'
  });

  sp.onAfterDrawScene.subscribe(_ => KaplanMeierRenderLines(subjIDColumn, xColumn, yColumn, groupColumn, sp));
  return sp;
}