import * as DG from 'datagrok-api/dg';
import {
    FitConfidenceIntervals,
    FitFunction,
    IFitChartData,
    IFitSeries,
    statisticsProperties
} from "@datagrok-libraries/statistics/src/fit/fit-curve";
import {
    getSeriesConfidenceInterval,
    getSeriesStatistics,
    LogOptions
} from "@datagrok-libraries/statistics/src/fit/fit-data";
import {Viewport} from "@datagrok-libraries/utils/src/transform";
import {FitConstants} from "./const";
import {BoxPlotStatistics, calculateBoxPlotStatistics} from "@datagrok-libraries/statistics/src/box-plot-statistics";
import {StringUtils} from "@datagrok-libraries/utils/src/string-utils";


interface FitRenderOptions {
    viewport: Viewport;
    ratio?: number;
}

interface FitPointRenderOptions extends FitRenderOptions {
    x: number;
    color: string;
    boxPlotStats: BoxPlotStatistics;
}

interface FitLineRenderOptions extends FitRenderOptions {
    logOptions: LogOptions;
    showAxes: boolean;
    showAxesLabels: boolean;
    screenBounds: DG.Rect;
    curveFunc?: (x: number) => number;
}

interface FitConfidenceIntervalRenderOptions extends FitLineRenderOptions {
    fitFunc?: FitFunction;
    userParamsFlag?: boolean;
    confidenceIntervals?: FitConfidenceIntervals;
    confidenceType?: string;
}

interface FitDroplineRenderOptions extends FitRenderOptions {
    showDroplines?: boolean;
    xValue: number;
    dataBounds: DG.Rect;
    curveFunc: (x: number) => number;
    logOptions: LogOptions;
}

interface FitStatisticsRenderOptions {
    statistics?: string[];
    fitFunc: FitFunction;
    logOptions: LogOptions;
    dataBox: DG.Rect;
}

interface FitTitleRenderOptions {
    showTitle: boolean;
    title?: string;
    dataBox: DG.Rect;
    screenBounds: DG.Rect;
}

interface FitAxesLabelsRenderOptions extends FitTitleRenderOptions {
    showAxesLabels: boolean;
    xAxisName?: string;
    yAxisName?: string;
}

interface FitLegendRenderOptions {
    showLegend?: boolean;
    dataBox: DG.Rect;
    ratio?: number;
}

interface FitLegendColumnlabelSeriesRenderOptions extends FitLegendRenderOptions {
    columnIdx: number;
    drawnCurvesInLegend: number;
    showColumnLabel?: boolean;
}


export function assignSeriesColors(series: IFitSeries, seriesIdx: number): void {
    const pointColor = series.pointColor ? DG.Color.fromHtml(series.pointColor) ? series.pointColor :
      DG.Color.toHtml(DG.Color.getCategoricalColor(seriesIdx)) : DG.Color.toHtml(DG.Color.getCategoricalColor(seriesIdx));
    const outlierColor = series.outlierColor ? DG.Color.fromHtml(series.outlierColor) ?
      series.outlierColor : DG.Color.toHtml(DG.Color.red) : DG.Color.toHtml(DG.Color.red);
    const fitLineColor = series.fitLineColor ? DG.Color.fromHtml(series.fitLineColor) ? series.fitLineColor :
      DG.Color.toHtml(DG.Color.getCategoricalColor(seriesIdx)) : DG.Color.toHtml(DG.Color.getCategoricalColor(seriesIdx));
    series.pointColor = pointColor;
    series.outlierColor = outlierColor;
    series.fitLineColor = fitLineColor;
}

export function renderConnectDots(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitRenderOptions): void {
    if (series.connectDots ?? false) {
        g.strokeStyle = series.pointColor!;
        g.lineWidth = 2 * renderOptions.ratio!;
        g.beginPath();
        for (let i = 0; i < series.points.length; i++) {
            const x = series.points[i].x;
            const y = series.points[i].y;
            const screenX = renderOptions.viewport.xToScreen(x);
            const screenY = renderOptions.viewport.yToScreen(y);
            if (i === 0)
                g.moveTo(screenX, screenY);
            else
                g.lineTo(screenX, screenY);
        }
        g.stroke();
    }
}

export function renderPoints(g: CanvasRenderingContext2D, series: IFitSeries, options: FitRenderOptions) {
    if (series.showPoints ?? 'points') {
        g.strokeStyle = series.pointColor!;
        if ((series.connectDots && series.showPoints !== '') || series.showPoints === 'points')
            drawPoints(g, series, options);
        else if (['candlesticks', 'both'].includes(series.showPoints!))
            drawCandles(g, series, options);
    }
}

/** Performs points drawing */
function drawPoints(g: CanvasRenderingContext2D, series: IFitSeries, options: FitRenderOptions): void {
    for (let i = 0; i < series.points.length!; i++) {
        const p = series.points[i];
        const color = series.connectDots ? series.pointColor! :
          p.outlier ? (p.outlierColor ? DG.Color.fromHtml(p.outlierColor) ? p.outlierColor : series.outlierColor! : series.outlierColor!) :
          (p.color ? DG.Color.fromHtml(p.color) ? p.color : series.pointColor! : series.pointColor!);
        const marker = p.marker ? p.marker as DG.MARKER_TYPE : series.markerType as DG.MARKER_TYPE;
        const outlierMarker = p.outlierMarker ? p.outlierMarker as DG.MARKER_TYPE : series.outlierMarkerType as DG.MARKER_TYPE;
        const size = !series.connectDots ? p.outlier ? FitConstants.OUTLIER_PX_SIZE * options.ratio! :
          p.size ? p.size : FitConstants.POINT_PX_SIZE * options.ratio! : FitConstants.POINT_PX_SIZE * options.ratio!;
        DG.Paint.marker(g, !series.connectDots ? p.outlier ? outlierMarker : marker : marker,
          options.viewport.xToScreen(p.x), options.viewport.yToScreen(p.y), color, size);
        if (p.stdev && !p.outlier) {
            g.strokeStyle = color;
            g.beginPath();
            g.moveTo(options.viewport.xToScreen(p.x), options.viewport.yToScreen(p.y + p.stdev));
            g.lineTo(options.viewport.xToScreen(p.x), options.viewport.yToScreen(p.y - p.stdev));
            g.stroke();
        }
    }
}

/** Performs candles drawing */
function drawCandles(g: CanvasRenderingContext2D, series: IFitSeries, options: FitRenderOptions) : void {
    for (let i = 0, candleStart = null; i < series.points.length!; i++) {
        const p = series.points[i];
        if (p.outlier)
            continue;
        const nextSame = i + 1 < series.points.length && series.points[i + 1].x === p.x;
        if (!candleStart && nextSame)
            candleStart = i;
        else if (candleStart !== null && !nextSame) {
            const values: number[] = [];
            for (let j = candleStart, ind = 0; j <= i; j++, ind++)
                values[ind] = series.points[j].y;
            const boxPlotStats = calculateBoxPlotStatistics(values);

            g.beginPath();
            drawCandlestick(g, {viewport: options.viewport, ratio: options.ratio!, x: p.x,
              boxPlotStats: boxPlotStats, color: series.pointColor!});
            g.stroke();

            if (series.showPoints === 'both') {
                for (let ind = 0; ind < values.length; ind++) {
                    if (values[ind] < boxPlotStats.lowerAdjacentValue || values[ind] > boxPlotStats.upperAdjacentValue) {
                        DG.Paint.marker(g, DG.MARKER_TYPE.OUTLIER,
                          options.viewport.xToScreen(p.x), options.viewport.yToScreen(values[ind]),
                          series.pointColor!, FitConstants.CANDLESTICK_OUTLIER_PX_SIZE * options.ratio!);
                    }
                }
            }

            candleStart = null;
        }
    }
}

/** Performs candlestick drawing */
function drawCandlestick(g: CanvasRenderingContext2D, renderOptions: FitPointRenderOptions): void {
    drawCandlestickBorder(g, renderOptions.x, renderOptions.boxPlotStats.lowerAdjacentValue, renderOptions.viewport);
    g.moveTo(renderOptions.viewport.xToScreen(renderOptions.x), renderOptions.viewport.yToScreen(renderOptions.boxPlotStats.lowerAdjacentValue));
    g.lineTo(renderOptions.viewport.xToScreen(renderOptions.x), renderOptions.viewport.yToScreen(renderOptions.boxPlotStats.upperAdjacentValue));
    drawCandlestickBorder(g, renderOptions.x, renderOptions.boxPlotStats.upperAdjacentValue, renderOptions.viewport);
    DG.Paint.marker(g, DG.MARKER_TYPE.CIRCLE, renderOptions.viewport.xToScreen(renderOptions.x),
      renderOptions.viewport.yToScreen(renderOptions.boxPlotStats.q2), renderOptions.color,
      FitConstants.CANDLESTICK_MEDIAN_PX_SIZE * renderOptions.ratio!);
}

/** Performs candlestick border drawing */
function drawCandlestickBorder(g: CanvasRenderingContext2D, x: number, adjacentValue: number, transform: Viewport): void {
    g.moveTo(transform.xToScreen(x) - (FitConstants.CANDLESTICK_BORDER_PX_SIZE / 2), transform.yToScreen(adjacentValue));
    g.lineTo(transform.xToScreen(x) + (FitConstants.CANDLESTICK_BORDER_PX_SIZE / 2), transform.yToScreen(adjacentValue));
}

export function renderFitLine(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitLineRenderOptions): void {
    if (series.showFitLine ?? true) {
        g.save();
        g.strokeStyle = series.fitLineColor!;
        g.lineWidth = 2 * renderOptions.ratio!;
        g.beginPath();
        if (series.lineStyle)
            g.setLineDash(FitConstants.LINE_STYLES[series.lineStyle]);

        const axesLeftPxMargin = renderOptions.showAxes ? renderOptions.showAxesLabels ?
          FitConstants.AXES_LEFT_PX_MARGIN_WITH_AXES_LABELS : FitConstants.AXES_LEFT_PX_MARGIN : 0;
        const axesRightPxMargin = renderOptions.showAxes ? FitConstants.AXES_RIGHT_PX_MARGIN : 0;
        for (let j = axesLeftPxMargin; j <= renderOptions.screenBounds.width - axesRightPxMargin; j++) {
            const x = renderOptions.screenBounds.x + j;
            const xForY = renderOptions.logOptions.logX ? Math.log10(renderOptions.viewport.xToWorld(x)) :
              renderOptions.viewport.xToWorld(x);
            const y = renderOptions.logOptions.logY ? renderOptions.viewport.yToScreen(Math.pow(10, renderOptions.curveFunc!(xForY))) :
              renderOptions.viewport.yToScreen(renderOptions.curveFunc!(xForY));
            if (j === axesLeftPxMargin)
                g.moveTo(x, y);
            else
                g.lineTo(x, y);
        }
        g.stroke();
        g.restore();
    }
}

export function renderConfidenceIntervals(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitConfidenceIntervalRenderOptions): void {
    if ((series.showFitLine ?? true) && (series.showCurveConfidenceInterval ?? false)) {
        // TODO: improve confidence intervals colors
        g.strokeStyle = series.confidenceIntervalColor ?? FitConstants.CONFIDENCE_INTERVAL_STROKE_COLOR;
        g.fillStyle = series.confidenceIntervalColor ?? FitConstants.CONFIDENCE_INTERVAL_FILL_COLOR;

        const confidenceIntervals = getSeriesConfidenceInterval(series, renderOptions.fitFunc!,
          renderOptions.userParamsFlag!, renderOptions.logOptions);
        drawConfidenceInterval(g, {viewport: renderOptions.viewport, logOptions: renderOptions.logOptions,
            showAxes: renderOptions.showAxes, showAxesLabels: renderOptions.showAxesLabels, screenBounds: renderOptions.screenBounds,
            confidenceIntervals: confidenceIntervals, confidenceType: FitConstants.CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP});
        drawConfidenceInterval(g, {viewport: renderOptions.viewport, logOptions: renderOptions.logOptions,
            showAxes: renderOptions.showAxes, showAxesLabels: renderOptions.showAxesLabels, screenBounds: renderOptions.screenBounds,
            confidenceIntervals: confidenceIntervals, confidenceType: FitConstants.CURVE_CONFIDENCE_INTERVAL_BOUNDS.BOTTOM});
        fillConfidenceInterval(g, {viewport: renderOptions.viewport, logOptions: renderOptions.logOptions,
            showAxes: renderOptions.showAxes, showAxesLabels: renderOptions.showAxesLabels, screenBounds: renderOptions.screenBounds,
            confidenceIntervals: confidenceIntervals});
    }
}

/** Performs a curve confidence interval drawing */
function drawConfidenceInterval(g: CanvasRenderingContext2D, renderOptions: FitConfidenceIntervalRenderOptions): void {
    g.beginPath();
    const axesLeftPxMargin = renderOptions.showAxes ? renderOptions.showAxesLabels ?
      FitConstants.AXES_LEFT_PX_MARGIN_WITH_AXES_LABELS : FitConstants.AXES_LEFT_PX_MARGIN : 0;
    const axesRightPxMargin = renderOptions.showAxes ? FitConstants.AXES_RIGHT_PX_MARGIN : 0;
    const confIntervalFunc = renderOptions.confidenceType === FitConstants.CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP ?
      renderOptions.confidenceIntervals!.confidenceTop : renderOptions.confidenceIntervals!.confidenceBottom;
    for (let i = axesLeftPxMargin; i <= renderOptions.screenBounds.width - axesRightPxMargin; i++) {
        const x = renderOptions.screenBounds.x + i;
        const xForY = renderOptions.logOptions.logX ? Math.log10(renderOptions.viewport.xToWorld(x)) : renderOptions.viewport.xToWorld(x);
        const y = renderOptions.logOptions.logY ? renderOptions.viewport.yToScreen(Math.pow(10, confIntervalFunc(xForY))) :
          renderOptions.viewport.yToScreen(confIntervalFunc(xForY));
        if (i === axesLeftPxMargin)
            g.moveTo(x, y);
        else
            g.lineTo(x, y);
    }
    g.stroke();
}

/** Performs a curve confidence interval color filling */
function fillConfidenceInterval(g: CanvasRenderingContext2D, renderOptions: FitConfidenceIntervalRenderOptions): void {
    g.beginPath();
    const axesLeftPxMargin = renderOptions.showAxes ? renderOptions.showAxesLabels ?
      FitConstants.AXES_LEFT_PX_MARGIN_WITH_AXES_LABELS : FitConstants.AXES_LEFT_PX_MARGIN : 0;
    const axesRightPxMargin = renderOptions.showAxes ? FitConstants.AXES_RIGHT_PX_MARGIN : 0;
    for (let i = axesLeftPxMargin; i <= renderOptions.screenBounds.width - axesRightPxMargin; i++) {
        const x = renderOptions.screenBounds.x + i;
        const xForY = renderOptions.logOptions.logX ? Math.log10(renderOptions.viewport.xToWorld(x)) :
          renderOptions.viewport.xToWorld(x);
        const y = renderOptions.logOptions.logY ? renderOptions.viewport.yToScreen(
          Math.pow(10, renderOptions.confidenceIntervals!.confidenceTop(xForY))) :
          renderOptions.viewport.yToScreen(renderOptions.confidenceIntervals!.confidenceTop(xForY));
        if (i === axesLeftPxMargin)
            g.moveTo(x, y);
        else
            g.lineTo(x, y);
    }

    // reverse traverse to make a shape of confidence interval to fill it
    for (let i = renderOptions.screenBounds.width - axesRightPxMargin; i >= axesLeftPxMargin; i--) {
        const x = renderOptions.screenBounds.x + i;
        const xForY = renderOptions.logOptions.logX ? Math.log10(renderOptions.viewport.xToWorld(x)) :
          renderOptions.viewport.xToWorld(x);
        const y = renderOptions.logOptions.logY ? renderOptions.viewport.yToScreen(
          Math.pow(10, renderOptions.confidenceIntervals!.confidenceBottom(xForY))) :
          renderOptions.viewport.yToScreen(renderOptions.confidenceIntervals!.confidenceBottom(xForY));
        g.lineTo(x, y);
    }
    g.closePath();
    g.fill();
}

export function renderDroplines(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitDroplineRenderOptions): void {
    if ((series.showFitLine ?? true) && series.droplines && renderOptions.showDroplines!) {
        g.save();
        g.strokeStyle = 'blue';
        g.lineWidth = renderOptions.ratio!;
        g.beginPath();
        g.setLineDash([5, 5]);
        for (let j = 0; j < series.droplines.length; j++) {
            const droplineName = series.droplines[j];
            if (droplineName === 'IC50')
                drawDropline(g, {viewport: renderOptions.viewport, xValue: series.parameters![2],
                  dataBounds: renderOptions.dataBounds, curveFunc: renderOptions.curveFunc, logOptions: renderOptions.logOptions});
        }
        g.stroke();
        g.restore();
    }
}

/** Performs a dropline drawing */
function drawDropline(g: CanvasRenderingContext2D, renderOptions: FitDroplineRenderOptions): void {
    if (renderOptions.logOptions.logX)
        renderOptions.xValue = Math.pow(10, renderOptions.xValue);
    const xForY = renderOptions.logOptions.logX ? Math.log10(renderOptions.xValue) : renderOptions.xValue;
    const y = renderOptions.logOptions.logY ? Math.pow(10, renderOptions.curveFunc(xForY)) : renderOptions.curveFunc(xForY);
    g.moveTo(renderOptions.viewport.xToScreen(renderOptions.dataBounds.minX), renderOptions.viewport.yToScreen(y));
    g.lineTo(renderOptions.viewport.xToScreen(renderOptions.xValue), renderOptions.viewport.yToScreen(y));
    g.lineTo(renderOptions.viewport.xToScreen(renderOptions.xValue), renderOptions.viewport.yToScreen(renderOptions.dataBounds.minY));
}

export function renderStatistics(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitStatisticsRenderOptions): void {
    if ((series.showFitLine ?? true) && renderOptions.statistics) {
        const statistics = getSeriesStatistics(series, renderOptions.fitFunc, renderOptions.logOptions);
        for (let i = 0; i < renderOptions.statistics.length; i++) {
            const statName = renderOptions.statistics[i];
            const prop = statisticsProperties.find((p) => p.name === statName);
            if (prop) {
                const s = StringUtils.formatNumber(prop.get(statistics));
                g.fillStyle = series.fitLineColor!;
                g.textAlign = 'left';
                g.fillText(prop.name + ': ' + s, renderOptions.dataBox.x + 5, renderOptions.dataBox.y + 20 + 20 * i);
            }
        }
    }
}

export function renderTitle(g: CanvasRenderingContext2D, renderOptions: FitTitleRenderOptions): void {
    if (renderOptions.showTitle) {
        g.font = '12px Roboto, "Roboto Local"';
        g.textAlign = 'center';
        g.fillStyle = 'black';
        g.fillText(renderOptions.title!, renderOptions.dataBox.midX - 5, renderOptions.screenBounds.y + 15);
    }
}

export function renderAxesLabels(g: CanvasRenderingContext2D, renderOptions: FitAxesLabelsRenderOptions): void {
    if (renderOptions.showAxesLabels) {
        g.font = '11px Roboto, "Roboto Local"';
        g.textAlign = 'center';
        g.fillStyle = 'black';
        g.fillText(renderOptions.xAxisName!, renderOptions.dataBox.midX - 5,
          renderOptions.screenBounds.maxY - FitConstants.X_AXIS_LABEL_BOTTOM_PX_MARGIN);
        g.translate(renderOptions.screenBounds.x, renderOptions.screenBounds.y);
        g.rotate(-Math.PI / 2);
        const axesTopPxMargin = renderOptions.showTitle ? FitConstants.AXES_TOP_PX_MARGIN_WITH_TITLE : FitConstants.AXES_TOP_PX_MARGIN;
        g.fillText(renderOptions.yAxisName!, -(renderOptions.dataBox.height / 2 + axesTopPxMargin + 15), 15);
        g.restore();
    }
}

export function renderLegend(g: CanvasRenderingContext2D, data: IFitChartData, renderOptions: FitLegendRenderOptions): void {
    if (renderOptions.showLegend) {
        g.font = '11px Roboto, "Roboto Local"';
        const columnNames = [...new Set(data.series?.map((series) => series.columnName))]
          .filter((colName) => colName !== null && colName !== undefined);
        let drawnCurvesInLegend = 0;
        for (let i = 0; i < columnNames.length; i++) {
            const colName = columnNames[i];
            if (data.chartOptions?.showColumnLabel)
              renderColumnLabel(g, colName!, {dataBox: renderOptions.dataBox, columnIdx: i, drawnCurvesInLegend});

            const series = data.series?.filter((series) => series.columnName === colName);
            for (let j = 0; j < series?.length!; j++) {
                const currentSeries = series![j];
                if (currentSeries.name === '' || currentSeries.name === null || currentSeries.name === undefined)
                    continue;
                renderLegendSeries(g, currentSeries, {dataBox: renderOptions.dataBox, ratio: renderOptions.ratio,
                  columnIdx: i, drawnCurvesInLegend, showColumnLabel: data.chartOptions?.showColumnLabel});
                drawnCurvesInLegend++;
            }
        }
        g.restore();
    }
}

function renderColumnLabel(g: CanvasRenderingContext2D, colName: string, renderOptions: FitLegendColumnlabelSeriesRenderOptions): void {
    g.beginPath();
    g.fillStyle = 'black';
    const colNameWidth = g.measureText(colName!).width;
    g.fillText(colName!, renderOptions.dataBox.maxX - colNameWidth, renderOptions.dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN +
      (renderOptions.columnIdx * FitConstants.LEGEND_RECORD_PX_HEIGHT + FitConstants.LEGEND_RECORD_PX_HEIGHT * renderOptions.drawnCurvesInLegend));
}

function renderLegendSeries(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitLegendColumnlabelSeriesRenderOptions): void {
    g.beginPath();
    g.strokeStyle = series.fitLineColor!;
    g.lineWidth = 2 * renderOptions.ratio!;
    const textWidth = g.measureText(series.name!).width;
    g.moveTo(renderOptions.dataBox.maxX - textWidth - FitConstants.LEGEND_RECORD_LINE_PX_WIDTH - FitConstants.LEGEND_RECORD_LINE_RIGHT_PX_MARGIN,
      renderOptions.dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN - FitConstants.LEGEND_RECORD_LINE_BOTTOM_PX_MARGIN +
      (renderOptions.showColumnLabel ? FitConstants.LEGEND_RECORD_PX_HEIGHT * (renderOptions.columnIdx + 1) : 0) +
      FitConstants.LEGEND_RECORD_PX_HEIGHT * renderOptions.drawnCurvesInLegend);
    g.lineTo(renderOptions.dataBox.maxX - textWidth - FitConstants.LEGEND_RECORD_LINE_RIGHT_PX_MARGIN,
      renderOptions.dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN - FitConstants.LEGEND_RECORD_LINE_BOTTOM_PX_MARGIN +
      (renderOptions.showColumnLabel ? FitConstants.LEGEND_RECORD_PX_HEIGHT * (renderOptions.columnIdx + 1) : 0) +
      FitConstants.LEGEND_RECORD_PX_HEIGHT * renderOptions.drawnCurvesInLegend);
    const marker = series.markerType ? series.markerType as DG.MARKER_TYPE : DG.MARKER_TYPE.CIRCLE;
    DG.Paint.marker(g, marker, renderOptions.dataBox.maxX - textWidth - FitConstants.LEGEND_RECORD_LINE_RIGHT_PX_MARGIN
      - FitConstants.LEGEND_RECORD_LINE_PX_WIDTH / 2, renderOptions.dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN -
      FitConstants.LEGEND_RECORD_LINE_BOTTOM_PX_MARGIN + (renderOptions.showColumnLabel ? FitConstants.LEGEND_RECORD_PX_HEIGHT
      * (renderOptions.columnIdx + 1) : 0) + FitConstants.LEGEND_RECORD_PX_HEIGHT * renderOptions.drawnCurvesInLegend,
      series.pointColor!, FitConstants.POINT_PX_SIZE * renderOptions.ratio!);
    g.fillStyle = series.fitLineColor!;
    g.fillText(series.name!, renderOptions.dataBox.maxX - textWidth,
      renderOptions.dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN + (renderOptions.showColumnLabel ?
      FitConstants.LEGEND_RECORD_PX_HEIGHT * (renderOptions.columnIdx + 1) : 0) + FitConstants.LEGEND_RECORD_PX_HEIGHT * renderOptions.drawnCurvesInLegend);
    g.stroke();
}
