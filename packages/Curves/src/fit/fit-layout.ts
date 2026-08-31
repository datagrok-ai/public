/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';

import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';

/** Chart geometry and what a cell of a given size has room for. */

export function inflateScreenBounds(rect: DG.Rect): DG.Rect {
  return rect.inflate(rect.width < FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH ? FitConstants.MIN_INFLATE_SIZE : FitConstants.INFLATE_SIZE,
    rect.height < FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT ? FitConstants.MIN_INFLATE_SIZE : FitConstants.INFLATE_SIZE);
}

/** Performs a chart layout, returning [viewport, xAxis, yAxis] */
export function layoutChart(rect: DG.Rect, showXAxisName: boolean, showYAxisName: boolean, showTitle: boolean):
  [DG.Rect, DG.Rect?, DG.Rect?] {
  if (rect.width < FitConstants.MIN_AXES_CELL_PX_WIDTH || rect.height < FitConstants.MIN_AXES_CELL_PX_HEIGHT)
    return [rect, undefined, undefined];
  // room is made for the name that is drawn, not for both because one of them is
  const axesLeftPxMargin = showYAxisName ? FitConstants.AXES_LEFT_PX_MARGIN_WITH_AXES_LABELS : FitConstants.AXES_LEFT_PX_MARGIN;
  const axesBottomPxMargin = showXAxisName ? FitConstants.AXES_BOTTOM_PX_MARGIN_WITH_AXES_LABELS : FitConstants.AXES_BOTTOM_PX_MARGIN;
  const axesTopPxMargin = showTitle ? FitConstants.AXES_TOP_PX_MARGIN_WITH_TITLE : FitConstants.AXES_TOP_PX_MARGIN;
  const rightPxMargin = FitConstants.AXES_RIGHT_PX_MARGIN;
  return [
    rect.cutLeft(axesLeftPxMargin).cutBottom(axesBottomPxMargin).cutTop(axesTopPxMargin).cutRight(rightPxMargin),
    rect.getBottom(axesBottomPxMargin).getTop(axesTopPxMargin).cutLeft(axesLeftPxMargin).cutRight(rightPxMargin),
    rect.getLeft(axesLeftPxMargin).getRight(FitConstants.AXES_RIGHT_PX_MARGIN).cutBottom(axesBottomPxMargin).cutTop(axesTopPxMargin),
  ];
}

export function areAxesShown(screenBounds: DG.Rect): boolean {
  return screenBounds.width >= FitConstants.MIN_AXES_CELL_PX_WIDTH && screenBounds.height >= FitConstants.MIN_AXES_CELL_PX_HEIGHT;
}

/** Naming one axis is not naming the other, so each is asked about on its own. */
export function shownAxesNames(screenBounds: DG.Rect, data: IFitChartData): {x: boolean, y: boolean} {
  const roomy = screenBounds.width >= FitConstants.MIN_X_AXIS_NAME_VISIBILITY_PX_WIDTH &&
    screenBounds.height >= FitConstants.MIN_Y_AXIS_NAME_VISIBILITY_PX_HEIGHT;
  return {x: roomy && !!data.chartOptions?.xAxisName, y: roomy && !!data.chartOptions?.yAxisName};
}

export function areAxesLabelsShown(screenBounds: DG.Rect, data: IFitChartData): boolean {
  const shown = shownAxesNames(screenBounds, data);
  return shown.x || shown.y;
}

export function isTitleShown(screenBounds: DG.Rect, data: IFitChartData): boolean {
  return screenBounds.width >= FitConstants.MIN_TITLE_PX_WIDTH && screenBounds.height >= FitConstants.MIN_TITLE_PX_HEIGHT &&
    !!data.chartOptions?.title;
}

export function isLegendShown(screenBounds: DG.Rect): boolean {
  return screenBounds.width >= FitConstants.MIN_LEGEND_PX_WIDTH && screenBounds.height >= FitConstants.MIN_LEGEND_PX_HEIGHT;
}

export function areDroplineLabelsShown(screenBounds: DG.Rect): boolean {
  return screenBounds.width >= FitConstants.MIN_DROPLINE_LABELS_PX_WIDTH &&
    screenBounds.height >= FitConstants.MIN_DROPLINE_LABELS_PX_HEIGHT;
}

export function areDroplinesShown(screenBounds: DG.Rect): boolean {
  return screenBounds.width >= FitConstants.MIN_DROPLINES_VISIBILITY_PX_WIDTH &&
    screenBounds.height >= FitConstants.MIN_DROPLINES_VISIBILITY_PX_HEIGHT;
}
