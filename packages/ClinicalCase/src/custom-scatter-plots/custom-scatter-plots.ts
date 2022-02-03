import { BaselineEndpointRenderLines, HysLawRenderLines, KaplanMeierRenderLines } from "./event-handlers";
import * as DG from "datagrok-api/dg";

export function createHysLawScatterPlot(dataframe, xColumn, yColumn, colorColumn, APColumnName = '') {
  const sp = DG.Viewer.scatterPlot(dataframe, {
    x: xColumn,
    y: yColumn,
    color: dataframe.col(colorColumn) ? colorColumn : '',
    legendPosition: 'Right'
  });

  sp.onAfterDrawScene.subscribe(_ => HysLawRenderLines(sp, 3, 1, 2, APColumnName));
  return sp;
}


export function createBaselineEndpointScatterPlot(dataframe, xColumn, yColumn, colorColumn, lowLimit, highLimit) {
  const sp = DG.Viewer.scatterPlot(dataframe, {
    x: xColumn,
    y: yColumn,
    color: dataframe.col(colorColumn) ? colorColumn : '',
    legendPosition: 'Right'
  });

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