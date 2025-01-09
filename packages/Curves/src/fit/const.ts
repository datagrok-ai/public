import {IFitSeries} from '@datagrok-libraries/statistics/src/fit/fit-curve';

export class FitConstants {
    static TAG_FIT_CHART_FORMAT = '.fitChartFormat';
    static TAG_FIT_CHART_FORMAT_3DX = '3dx';
    static FIT_SEM_TYPE = 'fit';
    static FIT_CELL_TYPE = 'fit';
    static TAG_FIT = '.fit';

    static CELL_DEFAULT_WIDTH = 230;
    static CELL_DEFAULT_HEIGHT = 160;
    static MIN_CELL_RENDERER_PX_WIDTH = 10;
    static MIN_CELL_RENDERER_PX_HEIGHT = 5;
    static MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH = 70;
    static MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT = 45;
    static MIN_AXES_CELL_PX_WIDTH = 100;
    static MIN_AXES_CELL_PX_HEIGHT = 55;
    static MIN_TITLE_PX_WIDTH = 275;
    static MIN_TITLE_PX_HEIGHT = 225;
    static MIN_X_AXIS_NAME_VISIBILITY_PX_WIDTH = 180;
    static MIN_Y_AXIS_NAME_VISIBILITY_PX_HEIGHT = 140;
    static MIN_LEGEND_PX_WIDTH = 325;
    static MIN_LEGEND_PX_HEIGHT = 275;
    static MIN_DROPLINES_VISIBILITY_PX_WIDTH = 120;
    static MIN_DROPLINES_VISIBILITY_PX_HEIGHT = 110;
    static AXES_LEFT_PX_MARGIN = 38;
    static AXES_LEFT_PX_MARGIN_WITH_AXES_LABELS = 45;
    static AXES_TOP_PX_MARGIN = 5;
    static AXES_TOP_PX_MARGIN_WITH_TITLE = 25;
    static AXES_RIGHT_PX_MARGIN = 18;
    static AXES_BOTTOM_PX_MARGIN = 15;
    static AXES_BOTTOM_PX_MARGIN_WITH_AXES_LABELS = 30;
    static X_AXIS_LABEL_BOTTOM_PX_MARGIN = 4;
    static LEGEND_TOP_PX_MARGIN = 10;

    static OUTLIER_PX_SIZE = 12;
    static POINT_PX_SIZE = 4;
    static OUTLIER_HITBOX_RADIUS = 2;
    static CANDLESTICK_BORDER_PX_SIZE = 4;
    static CANDLESTICK_MEDIAN_PX_SIZE = 3.5;
    static CANDLESTICK_OUTLIER_PX_SIZE = 6;
    static INFLATE_SIZE = -12 / 2;
    static MIN_INFLATE_SIZE = -4 / 2;
    static LEGEND_RECORD_PX_HEIGHT = 18;
    static LEGEND_RECORD_LINE_PX_WIDTH = 20;
    static LEGEND_RECORD_LINE_RIGHT_PX_MARGIN = 5;
    static LEGEND_RECORD_LINE_BOTTOM_PX_MARGIN = 4;

    static LINE_STYLES: {[key: string]: number[]} = {
        'solid': [],
        'dotted': [1, 1],
        'dashed': [5, 5],
        'dashdotted': [5, 5, 2, 5],
    };

    static CURVE_CONFIDENCE_INTERVAL_BOUNDS = {
        TOP: 'top',
        BOTTOM: 'bottom',
    };

    static CONFIDENCE_INTERVAL_STROKE_COLOR = 'rgba(255,191,63,0.4)';
    static CONFIDENCE_INTERVAL_FILL_COLOR = 'rgba(255,238,204,0.3)';

    static CONDITION_MAP: {[key: string]: ((series?: IFitSeries[]) => boolean)} = {
        'No series to show': (series) => series === undefined || series?.length === 0,
        'There were series with no points': (series) =>
          series!.some((series) => series.points?.length === 0),
        'There were series with all x zeroes': (series) =>
          series!.some((series) => series.points?.every((point) => point.x === 0)),
        'There were series with all y zeroes': (series) =>
          series!.some((series) => series.points?.every((point) => point.y === 0)),
    };
}
