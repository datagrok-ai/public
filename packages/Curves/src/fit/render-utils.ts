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
    dataPoints?: {x: number[], y: number[]};
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
    dataPoints?: {x: number[], y: number[]};
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
    series.pointColor = series.pointColor ? DG.Color.fromHtml(series.pointColor) ? series.pointColor :
      DG.Color.toHtml(DG.Color.getCategoricalColor(seriesIdx)) : DG.Color.toHtml(DG.Color.getCategoricalColor(seriesIdx));;
    series.outlierColor = series.outlierColor ? DG.Color.fromHtml(series.outlierColor) ?
      series.outlierColor : DG.Color.toHtml(DG.Color.red) : DG.Color.toHtml(DG.Color.red);
    series.fitLineColor = series.fitLineColor ? DG.Color.fromHtml(series.fitLineColor) ? series.fitLineColor :
      DG.Color.toHtml(DG.Color.getCategoricalColor(seriesIdx)) : DG.Color.toHtml(DG.Color.getCategoricalColor(seriesIdx));
}

export function renderConnectDots(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitRenderOptions): void {
    if (series.connectDots ?? false) {
        const viewport = renderOptions.viewport;
        g.strokeStyle = series.pointColor!;
        g.lineWidth = 2 * renderOptions.ratio!;
        g.beginPath();
        for (let i = 0; i < series.points.length; i++) {
            const x = series.points[i].x;
            const y = series.points[i].y;
            const screenX = viewport.xToScreen(x);
            const screenY = viewport.yToScreen(y);
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
    const ratio = options.ratio!;
    const defaultSize = FitConstants.POINT_PX_SIZE * ratio;
    const viewport = options.viewport;
    const connectDots = series.connectDots;
    const pointColor = series.pointColor!;
    const outlierColor = series.outlierColor!;

    for (let i = 0; i < series.points.length!; i++) {
        const p = series.points[i];
        const xScreen = viewport.xToScreen(p.x);
        const yScreen = viewport.yToScreen(p.y);
        const color = connectDots ? pointColor :
          p.outlier ? (p.outlierColor ? DG.Color.fromHtml(p.outlierColor) ? p.outlierColor : outlierColor : outlierColor) :
          (p.color ? DG.Color.fromHtml(p.color) ? p.color : pointColor : pointColor);
        const marker = p.marker ? p.marker as DG.MARKER_TYPE : series.markerType as DG.MARKER_TYPE;
        const outlierMarker = p.outlierMarker ? p.outlierMarker as DG.MARKER_TYPE : series.outlierMarkerType as DG.MARKER_TYPE;
        const size = !connectDots ? p.outlier ? FitConstants.OUTLIER_PX_SIZE * ratio : p.size ? p.size : defaultSize : defaultSize;
        const markerToDraw = !connectDots ? p.outlier ? outlierMarker : marker : marker;

        DG.Paint.marker(g, markerToDraw, xScreen, yScreen, color, size);
        if (p.stdev && !p.outlier) {
            g.strokeStyle = color;
            g.beginPath();
            g.moveTo(xScreen, viewport.yToScreen(p.y + p.stdev));
            g.lineTo(xScreen, viewport.yToScreen(p.y - p.stdev));
            g.stroke();
        }
    }
}

/** Performs candles drawing */
function drawCandles(g: CanvasRenderingContext2D, series: IFitSeries, options: FitRenderOptions) : void {
    const ratio = options.ratio!;
    const viewport = options.viewport;
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
            drawCandlestick(g, {viewport: viewport, ratio: ratio, x: p.x,
              boxPlotStats: boxPlotStats, color: series.pointColor!});
            g.stroke();

            if (series.showPoints === 'both') {
                for (let ind = 0; ind < values.length; ind++) {
                    if (values[ind] < boxPlotStats.lowerAdjacentValue || values[ind] > boxPlotStats.upperAdjacentValue) {
                        DG.Paint.marker(g, DG.MARKER_TYPE.OUTLIER,
                          viewport.xToScreen(p.x), viewport.yToScreen(values[ind]),
                          series.pointColor!, FitConstants.CANDLESTICK_OUTLIER_PX_SIZE * ratio);
                    }
                }
            }

            candleStart = null;
        }
    }
}

/** Performs candlestick drawing */
function drawCandlestick(g: CanvasRenderingContext2D, renderOptions: FitPointRenderOptions): void {
    const viewport = renderOptions.viewport;
    const x = renderOptions.x;
    const screenX = viewport.xToScreen(x);
    const boxPlotStats = renderOptions.boxPlotStats;

    drawCandlestickBorder(g, x, boxPlotStats.lowerAdjacentValue, viewport);
    g.moveTo(screenX, viewport.yToScreen(boxPlotStats.lowerAdjacentValue));
    g.lineTo(screenX, viewport.yToScreen(boxPlotStats.upperAdjacentValue));
    drawCandlestickBorder(g, x, boxPlotStats.upperAdjacentValue, viewport);
    DG.Paint.marker(g, DG.MARKER_TYPE.CIRCLE, screenX, viewport.yToScreen(boxPlotStats.q2), renderOptions.color,
      FitConstants.CANDLESTICK_MEDIAN_PX_SIZE * renderOptions.ratio!);
}

/** Performs candlestick border drawing */
function drawCandlestickBorder(g: CanvasRenderingContext2D, x: number, adjacentValue: number, transform: Viewport): void {
    const xScreen = transform.xToScreen(x);
    const yScreen = transform.yToScreen(adjacentValue);
    g.moveTo(xScreen - (FitConstants.CANDLESTICK_BORDER_PX_SIZE / 2), yScreen);
    g.lineTo(xScreen + (FitConstants.CANDLESTICK_BORDER_PX_SIZE / 2), yScreen);
}

/** Gets rendering variables */
function getRenderingVariables(renderOptions: FitConfidenceIntervalRenderOptions | FitLineRenderOptions) {
    const logX = renderOptions.logOptions.logX;
    const logY = renderOptions.logOptions.logY;
    const screenBounds = renderOptions.screenBounds;
    const axesLeftPxMargin = renderOptions.showAxes ? renderOptions.showAxesLabels ?
      FitConstants.AXES_LEFT_PX_MARGIN_WITH_AXES_LABELS : FitConstants.AXES_LEFT_PX_MARGIN : 0;
    const axesRightPxMargin = renderOptions.showAxes ? FitConstants.AXES_RIGHT_PX_MARGIN : 0;
    const xMin = screenBounds.x + axesLeftPxMargin;
    const xMax = screenBounds.x + screenBounds.width - axesRightPxMargin;
    const viewport = renderOptions.viewport;
    return {logX, logY, xMin, xMax, viewport};
}

export function renderFitLine(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitLineRenderOptions): void {
    if (series.showFitLine ?? true) {
        g.save();
        g.strokeStyle = series.fitLineColor!;
        g.lineWidth = 2 * renderOptions.ratio!;
        g.beginPath();
        if (series.lineStyle)
            g.setLineDash(FitConstants.LINE_STYLES[series.lineStyle]);
        const {logX, logY, xMin, xMax, viewport} = getRenderingVariables(renderOptions);
        const curveFunc = renderOptions.curveFunc!;

        for (let i = xMin; i <= xMax; i++) {
            const xForY = logX ? Math.log10(viewport.xToWorld(i)) : viewport.xToWorld(i);
            const y = logY ? viewport.yToScreen(Math.pow(10, curveFunc(xForY))) : viewport.yToScreen(curveFunc(xForY));
            if (i === xMin)
                g.moveTo(i, y);
            else
                g.lineTo(i, y);
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
        const showAxes = renderOptions.showAxes;
        const showAxesLabels = renderOptions.showAxesLabels;
        const logOptions = renderOptions.logOptions;
        const viewport = renderOptions.viewport;
        const screenBounds = renderOptions.screenBounds;
        const dataPoints = series.dataPoints;

        const confidenceIntervals = getSeriesConfidenceInterval(series, renderOptions.fitFunc!, renderOptions.userParamsFlag!, dataPoints, logOptions);
        drawConfidenceInterval(g, {viewport: viewport, logOptions: logOptions, showAxes: showAxes, showAxesLabels: showAxesLabels,
            screenBounds: screenBounds, confidenceIntervals: confidenceIntervals, confidenceType: FitConstants.CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP});
        drawConfidenceInterval(g, {viewport: viewport, logOptions: logOptions, showAxes: showAxes, showAxesLabels: showAxesLabels,
            screenBounds: screenBounds, confidenceIntervals: confidenceIntervals, confidenceType: FitConstants.CURVE_CONFIDENCE_INTERVAL_BOUNDS.BOTTOM});
        fillConfidenceInterval(g, {viewport: viewport, logOptions: logOptions, showAxes: showAxes, showAxesLabels: showAxesLabels,
            screenBounds: screenBounds, confidenceIntervals: confidenceIntervals});
    }
}

/** Performs a curve confidence interval drawing */
function drawConfidenceInterval(g: CanvasRenderingContext2D, renderOptions: FitConfidenceIntervalRenderOptions): void {
    g.beginPath();

    const {logX, logY, xMin, xMax, viewport} = getRenderingVariables(renderOptions);
    const confIntervalFunc = renderOptions.confidenceType === FitConstants.CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP ?
      renderOptions.confidenceIntervals!.confidenceTop : renderOptions.confidenceIntervals!.confidenceBottom;

    for (let i = xMin; i <= xMax; i++) {
        const xForY = logX ? Math.log10(viewport.xToWorld(i)) : viewport.xToWorld(i);
        const y = logY ? viewport.yToScreen(Math.pow(10, confIntervalFunc(xForY))) : viewport.yToScreen(confIntervalFunc(xForY));
        if (i === xMin)
            g.moveTo(i, y);
        else
            g.lineTo(i, y);
    }
    g.stroke();
}

/** Performs a curve confidence interval color filling */
function fillConfidenceInterval(g: CanvasRenderingContext2D, renderOptions: FitConfidenceIntervalRenderOptions): void {
    g.beginPath();

    const {logX, logY, xMin, xMax, viewport} = getRenderingVariables(renderOptions);
    const confTop = renderOptions.confidenceIntervals!.confidenceTop;
    const confBottom = renderOptions.confidenceIntervals!.confidenceBottom;

    for (let i = xMin; i <= xMax; i++) {
        const xForY = logX ? Math.log10(viewport.xToWorld(i)) : viewport.xToWorld(i);
        const y = logY ? viewport.yToScreen(Math.pow(10, confTop(xForY))) : viewport.yToScreen(confTop(xForY));
        if (i === xMin)
            g.moveTo(i, y);
        else
            g.lineTo(i, y);
    }

    // reverse traverse to make a shape of confidence interval to fill it
    for (let i = xMax; i >= xMin; i--) {
        const xForY = logX ? Math.log10(viewport.xToWorld(i)) : viewport.xToWorld(i);
        const y = logY ? viewport.yToScreen(Math.pow(10, confBottom(xForY))) : viewport.yToScreen(confBottom(xForY));
        g.lineTo(i, y);
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
        const viewport = renderOptions.viewport;
        const dataBounds = renderOptions.dataBounds;
        const curveFunc = renderOptions.curveFunc;
        const logOptions = renderOptions.logOptions;
        for (let j = 0; j < series.droplines.length; j++) {
            const droplineName = series.droplines[j];
            if (droplineName === 'IC50')
                drawDropline(g, {viewport: viewport, xValue: series.parameters![2], dataBounds: dataBounds,
                  curveFunc: curveFunc, logOptions: logOptions});
        }
        g.stroke();
        g.restore();
    }
}

/** Performs a dropline drawing */
function drawDropline(g: CanvasRenderingContext2D, renderOptions: FitDroplineRenderOptions): void {
    const logX = renderOptions.logOptions.logX;
    const logY = renderOptions.logOptions.logY;
    const xValue = logX ? Math.pow(10, renderOptions.xValue) : renderOptions.xValue;
    const viewport = renderOptions.viewport;
    const dataBounds = renderOptions.dataBounds;

    const xForY = logX ? Math.log10(xValue) : xValue;
    const y = logY ? Math.pow(10, renderOptions.curveFunc(xForY)) : renderOptions.curveFunc(xForY);
    const screenX = viewport.xToScreen(xValue);
    const screenY = viewport.yToScreen(y);
    g.moveTo(viewport.xToScreen(dataBounds.minX), screenY);
    g.lineTo(screenX, screenY);
    g.lineTo(screenX, viewport.yToScreen(dataBounds.minY));
}

export function renderStatistics(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitStatisticsRenderOptions): void {
    if ((series.showFitLine ?? true) && renderOptions.statistics && renderOptions.statistics.length > 0) {
        const dataBox = renderOptions.dataBox;
        const stats = renderOptions.statistics;
        const dataPoints = renderOptions.dataPoints;
        const statistics = getSeriesStatistics(series, renderOptions.fitFunc, dataPoints, renderOptions.logOptions);
        for (let i = 0; i < stats.length; i++) {
            const statName = stats[i];
            const prop = statisticsProperties.find((p) => p.name === statName);
            if (prop) {
                const s = StringUtils.formatNumber(prop.get(statistics));
                g.fillStyle = series.fitLineColor!;
                g.textAlign = 'left';
                g.fillText(prop.name + ': ' + s, dataBox.x + 5, dataBox.y + 20 + 20 * i);
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
        const screenBounds = renderOptions.screenBounds;
        const dataBox = renderOptions.dataBox;
        g.font = '11px Roboto, "Roboto Local"';
        g.textAlign = 'center';
        g.fillStyle = 'black';
        g.fillText(renderOptions.xAxisName!, dataBox.midX - 5, screenBounds.maxY - FitConstants.X_AXIS_LABEL_BOTTOM_PX_MARGIN);
        g.translate(screenBounds.x, screenBounds.y);
        g.rotate(-Math.PI / 2);
        const axesTopPxMargin = renderOptions.showTitle ? FitConstants.AXES_TOP_PX_MARGIN_WITH_TITLE : FitConstants.AXES_TOP_PX_MARGIN;
        g.fillText(renderOptions.yAxisName!, -(dataBox.height / 2 + axesTopPxMargin + 15), 15);
        g.restore();
    }
}

export function renderLegend(g: CanvasRenderingContext2D, data: IFitChartData, renderOptions: FitLegendRenderOptions): void {
    if (renderOptions.showLegend) {
        const dataBox = renderOptions.dataBox;
        const ratio = renderOptions.ratio;
        g.font = '11px Roboto, "Roboto Local"';
        const columnNames = [...new Set(data.series?.map((series) => series.columnName))]
          .filter((colName) => colName !== null && colName !== undefined);
        let drawnCurvesInLegend = 0;
        for (let i = 0; i < columnNames.length; i++) {
            const colName = columnNames[i];
            if (data.chartOptions?.showColumnLabel)
              renderColumnLabel(g, colName!, {dataBox: dataBox, columnIdx: i, drawnCurvesInLegend});

            const series = data.series?.filter((series) => series.columnName === colName);
            for (let j = 0; j < series?.length!; j++) {
                const currentSeries = series![j];
                if (currentSeries.name === '' || currentSeries.name === null || currentSeries.name === undefined)
                    continue;
                renderLegendSeries(g, currentSeries, {dataBox: dataBox, ratio: ratio,
                  columnIdx: i, drawnCurvesInLegend, showColumnLabel: data.chartOptions?.showColumnLabel});
                drawnCurvesInLegend++;
            }
        }
        g.restore();
    }
}

function renderColumnLabel(g: CanvasRenderingContext2D, colName: string, renderOptions: FitLegendColumnlabelSeriesRenderOptions): void {
    const dataBox = renderOptions.dataBox;
    g.beginPath();
    g.fillStyle = 'black';
    const colNameWidth = g.measureText(colName!).width;
    g.fillText(colName!, dataBox.maxX - colNameWidth, dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN +
      (renderOptions.columnIdx * FitConstants.LEGEND_RECORD_PX_HEIGHT + FitConstants.LEGEND_RECORD_PX_HEIGHT * renderOptions.drawnCurvesInLegend));
}

function renderLegendSeries(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitLegendColumnlabelSeriesRenderOptions): void {
    const dataBox = renderOptions.dataBox;
    const columnIdx = renderOptions.columnIdx;
    const ratio = renderOptions.ratio!;
    const showColumnLabel = renderOptions.showColumnLabel;
    const drawnCurvesInLegend = renderOptions.drawnCurvesInLegend;
    g.beginPath();
    g.strokeStyle = series.fitLineColor!;
    g.lineWidth = 2 * ratio;
    const textWidth = g.measureText(series.name!).width;
    g.moveTo(dataBox.maxX - textWidth - FitConstants.LEGEND_RECORD_LINE_PX_WIDTH - FitConstants.LEGEND_RECORD_LINE_RIGHT_PX_MARGIN,
      dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN - FitConstants.LEGEND_RECORD_LINE_BOTTOM_PX_MARGIN +
      (showColumnLabel ? FitConstants.LEGEND_RECORD_PX_HEIGHT * (columnIdx + 1) : 0) +
      FitConstants.LEGEND_RECORD_PX_HEIGHT * drawnCurvesInLegend);
    g.lineTo(dataBox.maxX - textWidth - FitConstants.LEGEND_RECORD_LINE_RIGHT_PX_MARGIN,
      dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN - FitConstants.LEGEND_RECORD_LINE_BOTTOM_PX_MARGIN +
      (showColumnLabel ? FitConstants.LEGEND_RECORD_PX_HEIGHT * (columnIdx + 1) : 0) +
      FitConstants.LEGEND_RECORD_PX_HEIGHT * drawnCurvesInLegend);
    const marker = series.markerType ? series.markerType as DG.MARKER_TYPE : DG.MARKER_TYPE.CIRCLE;
    DG.Paint.marker(g, marker, dataBox.maxX - textWidth - FitConstants.LEGEND_RECORD_LINE_RIGHT_PX_MARGIN
      - FitConstants.LEGEND_RECORD_LINE_PX_WIDTH / 2, dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN -
      FitConstants.LEGEND_RECORD_LINE_BOTTOM_PX_MARGIN + (showColumnLabel ? FitConstants.LEGEND_RECORD_PX_HEIGHT
      * (columnIdx + 1) : 0) + FitConstants.LEGEND_RECORD_PX_HEIGHT * drawnCurvesInLegend,
      series.pointColor!, FitConstants.POINT_PX_SIZE * ratio);
    g.fillStyle = series.fitLineColor!;
    g.fillText(series.name!, dataBox.maxX - textWidth,
      dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN + (showColumnLabel ?
      FitConstants.LEGEND_RECORD_PX_HEIGHT * (columnIdx + 1) : 0) + FitConstants.LEGEND_RECORD_PX_HEIGHT * drawnCurvesInLegend);
    g.stroke();
}
