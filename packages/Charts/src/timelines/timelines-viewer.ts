import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import * as echarts from 'echarts';
import { format } from 'echarts/lib/util/time';
import $ from 'cash-dom';

import { EChartViewer } from '../echart-viewer';
import { options, deepCopy } from './echarts-options';
import { ColumnData, ColumnsData, Indexable, markerPosition, markerType, timePoint, visibilityMode, VISIBILITY_MODE } from './constants';


export class TimelinesViewer extends EChartViewer {
  splitByColumnName: string;
  startColumnName: string;
  endColumnName: string;
  colorByColumnName: string;
  showOpenIntervals: boolean;
  eventColumnName: string;
  eventsColumnNames: string[];
  showEventInTooltip: boolean;
  marker: markerType;
  markerSize: number;
  markerPosition: markerPosition;
  lineWidth: number;
  dateFormat: string;
  axisPointer: string;
  showZoomSliders: boolean;
  legendVisibility: visibilityMode;
  splitByRegexps: RegExp[];
  colorByRegexps: RegExp[];
  eventRegexps: any[];
  startRegexps: RegExp[];
  endRegexps: RegExp[];
  defaultDateFormat: string;
  data: any[]; //(string | timePoint | string[] | boolean)[];
  columnData: ColumnsData;
  count: number;
  selectionColor: string;
  defaultColor: string;
  zoomState: number[][];
  tooltipOffset: number;
  initialized: boolean;
  titleDiv: HTMLDivElement = ui.div();
  legendDiv: HTMLDivElement = ui.div();
  colorMap: Indexable | null = null;
  dataMax: any;

  constructor() {
    super();

    this.splitByColumnName = this.string('splitByColumnName');
    this.startColumnName = this.string('startColumnName');
    this.endColumnName = this.string('endColumnName');
    this.colorByColumnName = this.string('colorByColumnName');
    this.showOpenIntervals = this.bool('showOpenIntervals', false);
    this.eventColumnName = this.string('eventColumnName');
    this.eventsColumnNames = this.addProperty('eventsColumnNames', DG.TYPE.COLUMN_LIST);
    this.showEventInTooltip = this.bool('showEventInTooltip', true);

    this.marker = <markerType>this.string('marker', 'circle', { choices: ['circle', 'rect', 'ring', 'diamond'] });
    this.markerSize = this.int('markerSize', 6);
    this.markerPosition = <markerPosition>this.string('markerPosition', 'main line',
      { choices: ['main line', 'above main line', 'scatter'] });
    this.lineWidth = this.int('lineWidth', 3);
    this.dateFormat = this.string('dateFormat'); // TODO: add an extendable dropdown
    this.axisPointer = this.string('axisPointer', 'shadow',
      { choices: ['cross', 'line', 'shadow', 'none'] });
    this.showZoomSliders = this.bool('showZoomSliders', true);
    this.legendVisibility = <visibilityMode>this.string('legendVisibility', VISIBILITY_MODE.AUTO,
      { choices: Object.values(VISIBILITY_MODE) });

    this.splitByRegexps = [/^USUBJID$/, /id/i];
    this.colorByRegexps = [/^([A-Z]{2}(TERM|TEST|TRT|VAL)|(ACT)?ARM|MIDS(TYPE)?|VISIT)$/, /event/i];
    this.eventRegexps = [/^event$/, ...this.colorByRegexps];
    this.startRegexps = [/^((VISIT|[A-Z]{2}(ST)?)DY)$/, /start|begin/i];
    this.endRegexps = [/^((VISIT|[A-Z]{2}(EN)?)DY)$/, /stop|end/i];

    this.defaultDateFormat = '{MMM} {d}';
    this.data = [];
    this.columnData = {};
    this.count = 0;
    this.selectionColor = DG.Color.toRgb(DG.Color.selectedRows);
    this.defaultColor = DG.Color.toRgb(DG.Color.histogramBar);
    this.zoomState = [[0, 100], [0, 100], [0, 100], [0, 100]];
    this.tooltipOffset = 10;
    this.initialized = false;
    this.option = deepCopy(options);
  }

  init() {
    if (!this.initialized) {
      this.helpUrl = 'https://raw.githubusercontent.com/datagrok-ai/public/master/packages/Charts/README.md#timelines';

      this.updateZoom();
      this.chart.on('dataZoom', () => {
        this.chart.getOption().dataZoom!.forEach((z, i) => {
          this.zoomState[i][0] = z.start!;
          this.zoomState[i][1] = z.end!;
          if (z.type === 'slider' && Object.keys(z).includes('yAxisIndex')) {
            this.lineWidth = z.end! - z.start! < 60 ? z.end! - z.start! < 30 ? 3 : 2 : 1;
            this.markerSize = z.end! - z.start! < 60 ? z.end! - z.start! < 30 ? 6 : 4 : 3;
          }
        });
      });

      this.subs.push(this.onEvent('d4-context-menu').subscribe((data) => {
          data.args.menu.item('Reset View', () => {
            this.zoomState = [[0, 100], [0, 100], [0, 100], [0, 100]];
            this.render();
          });
      }));

      this.chart.on('rendered', () => {
        this.count = 0;
      });

      this.chart.on('click', (params: Indexable) => this.dataFrame.selection.handleClick((i) => {
        if (params.componentType === 'yAxis')
          return this.getStrValue(this.columnData.splitByColumnName!, i) === params.value;
        if (params.componentType === 'series') {
          return params.value[0] === this.getStrValue(this.columnData.splitByColumnName!, i) &&
                 this.columnData.startColumnName && this.columnData.endColumnName ?
                 (this.isSameDate(params.value[1], this.getSafeValue(this.columnData.startColumnName, i)) &&
                 this.isSameDate(params.value[2], this.getSafeValue(this.columnData.endColumnName, i)))
                 : Object.values(this.columnData.eventsColumnNames!).some((c) =>
                 this.isSameDate(params.value[1], this.getSafeValue(c, i)));
        }
        return false;
      }, params.event!.event));

      this.chart.on('mouseover', (params: Indexable) => this.getTooltip(params));

      this.chart.on('mouseout', () => ui.tooltip.hide());

      this.option.tooltip.axisPointer.type = this.axisPointer;

      this.root.appendChild(this.titleDiv);
      this.root.appendChild(this.legendDiv);
      this.initialized = true;
    }
  }

  onPropertyChanged(property: DG.Property) {
    if (!this.initialized) return;
    if (property.name === 'axisPointer')
      this.option.tooltip.axisPointer.type = property.get(this);
    else if (property.name === 'showZoomSliders') {
      (this.option.dataZoom as echarts.EChartOption.DataZoom[]).forEach((z) => {
        if (z.type === 'slider') (<echarts.EChartOption.DataZoom.Slider>z).show = this.showZoomSliders;
      });
    } else if (property.name.endsWith('ColumnName') || property.name.endsWith('ColumnNames')) {
      if (property.get(this)) {
        const columnData = this.updateColumnData(property);
        if (property.name === 'colorByColumnName') {
          this.colorMap = this.getColorMap(columnData!.categories);
          this.updateLegend(columnData!.column);
          this.switchLegendVisibility(this.legendVisibility);
        }

        if (property.name === 'eventsColumnNames' || this.eventsColumnNames?.length > 0) {
          this.columnData.startColumnName = null;
          this.columnData.endColumnName = null;
        }
      } else {
        this.columnData[property.name] = null;
        if (property.name === 'colorByColumnName') {
          this.hideLegend();
          this.colorMap = null;
        }
      }
    } else if (property.name === 'legendVisibility' && this.colorByColumnName)
      this.switchLegendVisibility(property.get(this));

    this.render();
  }

  formatDate(value: timePoint): string {
    return value instanceof Date ? format(value, this.dateFormat, false) : value;
  }

  getTooltip(params: Indexable): void {
    const x = params.event.event.x + this.tooltipOffset;
    const y = params.event.event.y + this.tooltipOffset;
    if (this.showEventInTooltip) {
      const tooltipContent = params.componentType === 'yAxis' ? ui.div(`${params.value}`) :
        ui.divV([ui.div(`key: ${params.value[0]}`),
          ui.div(`event: ${params.value[4]}`),
          ui.div(`start: ${this.formatDate(params.value[1])}`),
          ui.div(`end: ${this.formatDate(params.value[2])}`),
        ]);
      ui.tooltip.show(tooltipContent, x, y);
    } else {
      ui.tooltip.showRowGroup(this.dataFrame, (i) => {
        if (params.componentType === 'yAxis')
          return this.getStrValue(this.columnData.splitByColumnName!, i) === params.value;
        if (params.componentType === 'series') {
          return params.value[0] === this.getStrValue(this.columnData.splitByColumnName!, i) &&
            this.columnData.startColumnName && this.columnData.endColumnName ?
            (this.isSameDate(params.value[1], this.getSafeValue(this.columnData.startColumnName, i)) &&
            this.isSameDate(params.value[2], this.getSafeValue(this.columnData.endColumnName, i)))
            : Object.values(this.columnData.eventsColumnNames!).some((c) =>
            this.isSameDate(params.value[1], this.getSafeValue(c, i)));
        }
        return false;
      }, x, y);
    }
  }

  getColorMap(categories: string[]) {
    return categories.reduce((colorMap, c, i) => {
      colorMap[c] = DG.Color.toRgb(DG.Color.getCategoricalColor(i));
      return colorMap;
    }, <Indexable>{});
  }

  onTableAttached(): void {
    this.init();

    const columns = this.dataFrame.columns.toList();

    const strColumns = columns.filter((col) => col.type === DG.COLUMN_TYPE.STRING)
      .sort((a, b) => a.categories.length - b.categories.length);

    const numColumns = [...this.dataFrame.columns.numerical].sort((a, b) => a.stats.avg - b.stats.avg);
    const numericalTypes = [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.DATE_TIME];

    if (strColumns.length < 1 || numColumns.length < 1) {
      this.showErrorMessage('Not enough data to produce the result.');
      return;
    }

    this.splitByColumnName = (this.findColumn(columns, this.splitByRegexps) || strColumns[strColumns.length - 1]).name;
    this.startColumnName = (this.findColumn(columns, this.startRegexps, numericalTypes) || numColumns[0]).name;
    this.endColumnName = (this.findColumn(columns, this.endRegexps, numericalTypes) ||
      numColumns[numColumns.length - 1]).name;
    this.colorByColumnName = (this.findColumn(columns, this.colorByRegexps, [DG.COLUMN_TYPE.STRING]) ||
      strColumns[0]).name;
    this.eventColumnName = (this.findColumn(columns, this.eventRegexps, [DG.COLUMN_TYPE.STRING]) || strColumns[0]).name;

    const columnPropNames = [
      'splitByColumnName', 'startColumnName', 'endColumnName', 'colorByColumnName', 'eventColumnName',
    ];

    const columnNames = [
      this.splitByColumnName, this.startColumnName, this.endColumnName, this.colorByColumnName, this.eventColumnName,
    ];

    this.columnData = columnPropNames.reduce((map, v, i) => {
      const column = this.dataFrame.getCol(columnNames[i]);
      map[v] = this.getColumnData(column);
      return map;
    }, <Indexable>{});
    if (this.eventsColumnNames?.length > 0) {
      this.columnData['eventsColumnNames'] = {};
      for (const columnName of this.eventsColumnNames) {
        const column = this.dataFrame.col(columnName);
        if (column == null)
          continue;
        this.columnData['eventsColumnNames'][columnName] = this.getColumnData(column);
      }
      this.columnData.startColumnName = null;
      this.columnData.endColumnName = null;
    }

    this.colorMap = this.getColorMap(this.columnData.colorByColumnName!.categories!);
    this.updateLegend(this.columnData.colorByColumnName!.column);
    this.switchLegendVisibility(this.legendVisibility);

    let prevSubj: string | null = null;

    (this.option.series[0] as echarts.EChartOption.SeriesCustom).renderItem = (params, api) => {
      let overlap = false;
      if (params.dataIndex! > 0) {
        const prev = this.data[params.dataIndex! - 1];
        const curSubj = this.data[params.dataIndex!][0];
        if (curSubj === prev[0] &&
            prev[1] && prev[2] && prev[1] !== prev[2] &&
            api.value!(1) < prev[2]) {
          overlap = true;
          if (prevSubj !== curSubj) {
            this.count = 0;
            prevSubj = curSubj;
          }
        }
      }

      const categoryIndex = api.value!(0);
      const start = api.coord!([api.value!(1), categoryIndex]);
      const end = api.coord!([api.value!(2), categoryIndex]);
      const width = end[0] - start[0];

      const group: echarts.EChartOption.SeriesCustom.RenderItemReturnGroup = {
        type: 'group',
        children: [],
      };

      if (isNaN(api.value!(1)) || isNaN(api.value!(2)) || this.markerSize > width) {
        const xPos = (shift: number): number => isNaN(start[0]) ? end[0] : start[0] - shift;
        const yPos = (shift: number): number => end[1] - (this.markerPosition === 'main line' ? shift :
          this.markerPosition === 'above main line' ? Math.max(this.markerSize, this.lineWidth) + shift :
            ((params.dataIndex! % 2) * 2 - 1)*(this.markerSize * 3));

        // echarts.EChartOption.SeriesCustom.RenderItemReturnCircle | echarts.EChartOption.SeriesCustom.RenderItemReturnRect | echarts.EChartOption.SeriesCustom.RenderItemReturnRing
        const marker: any = {
          type: this.marker,
          shape: this.marker === 'circle' ? {
            cx: xPos(0),
            cy: yPos(0),
            r: this.markerSize / 2,
          } : this.marker === 'ring' ? {
            cx: xPos(0),
            cy: yPos(0),
            r: this.markerSize / 2,
            r0: this.markerSize / 4,
          } : {
            x: xPos(this.markerSize / 2),
            y: yPos(this.markerSize / 2),
            width: this.markerSize,
            height: this.markerSize,
          },
          style: {
            fill: api.value!(5) ? this.selectionColor :
              this.colorMap ? this.colorMap[isNaN(api.value!(3)) ? this.data[params.dataIndex!][3][0] : api.value!(3)] :
                this.defaultColor,
          },
        };

        if (this.marker === 'diamond') {
          marker.type = 'rect';
          marker.x = xPos(0);
          marker.y = yPos(0);
          marker.shape.x = -this.markerSize / 2;
          marker.shape.y = -this.markerSize / 2;
          marker.shape.r = this.markerSize / 4;
          marker.rotation = 0.785398;
        } else if (this.marker === 'rect') {
          marker.x = 0;
          marker.y = 0;
          marker.shape.x = xPos(this.markerSize / 2);
          marker.shape.y = yPos(this.markerSize / 2);
          marker.shape.r = 0;
          marker.rotation = 0;
        }

        group.children!.push(marker);
      } else {
        const rectShape = echarts.graphic.clipRectByRect({
          x: start[0],
          y: start[1] - this.lineWidth / 2,
          width: width,
          height: this.lineWidth,
        }, {
          x: params.coordSys!.x!,
          y: params.coordSys!.y!,
          width: params.coordSys!.width!,
          height: params.coordSys!.height!,
        });

        if (overlap) {
          const height = api.size!([0, 1])[1];
          const offset = Math.max(this.markerSize * 2, this.lineWidth);
          // Shift along the Y axis
          rectShape.y += (this.count % 3) ? (this.count % 3 === 2) ?
            0 : offset-height/2 : height/2-offset;
          this.count += 1;
        }

        group.children!.push({
          type: 'rect',
          transition: ['shape'],
          shape: rectShape,
          style: {
            fill: api.value!(5) ? this.selectionColor :
              this.colorMap ? this.colorMap[isNaN(api.value!(3)) ? this.data[params.dataIndex!][3][0] : api.value!(3)] :
                this.defaultColor,
          },
        });
      }

      return group;
    };

    super.onTableAttached();
  }

  /** Find a column based on provided name patterns and types.
   * The list of patterns should be sorted by priority in descending order. */
  findColumn(columns: DG.Column[], regexps: RegExp[], types: DG.ColumnType[] = []): DG.Column | null {
    if (types.length)
      columns = columns.filter((c) => types.includes(c.type));

    for (const regex of regexps) {
      const column = columns.find((c) => c.name.match(regex));
      if (column) return column;
    }
    return null;
  }

  getColumnData(column: DG.Column): ColumnData {
    return {
      column,
      data: column.getRawData(),
      categories: column.type === DG.COLUMN_TYPE.STRING ? column.categories : null,
    };
  }

  updateColumnData(prop: DG.Property): ColumnData | ColumnsData | null {
    if (prop.name.endsWith('ColumnNames')) {
      this.columnData[prop.name] = {};
      for (const columnName of prop.get(this)) {
        const column = this.dataFrame.col(columnName);
        if (column == null)
          return null;
        this.columnData[prop.name][columnName] = this.getColumnData(column);
      }
    } else {
      const column = this.dataFrame.col(prop.get(this));
      if (column == null)
        return null;
      this.columnData[prop.name] = this.getColumnData(column);
    }
    return this.columnData[prop.name];
  }

  updateZoom(): void {
    (this.option.dataZoom as echarts.EChartOption.DataZoom[]).forEach((z, i) => {
      z.start = this.zoomState[i][0];
      z.end = this.zoomState[i][1];
    });
  }

  updateLegend(column: DG.Column): void {
    $(this.legendDiv).empty();
    const legend = DG.Legend.create(column);
    this.legendDiv.appendChild(legend.root);
    $(legend.root).addClass('charts-legend');
  }

  showLegend(): void {
    $(this.legendDiv).show();
    $(this.chart.getDom()).css('marginRight', '100px');
  }

  hideLegend(): void {
    $(this.legendDiv).hide();
    $(this.chart.getDom()).css('marginRight', '');
  }

  switchLegendVisibility(mode: visibilityMode): void {
    const { column, categories } = this.columnData.colorByColumnName as ColumnData;
    const autoShow = column.matches(DG.TYPE.CATEGORICAL) && categories!.length < 100;
    if (mode === VISIBILITY_MODE.ALWAYS || (mode === VISIBILITY_MODE.AUTO && autoShow))
      this.showLegend();
    else
      this.hideLegend();
  }

  getStrValue(columnData: ColumnData, idx: number): string {
    const { column, categories, data } = columnData;
    return column.type === DG.COLUMN_TYPE.STRING ? categories![data[idx]] : column.getString(idx);
  }

  getSafeValue(columnData: ColumnData, idx: number): timePoint {
    const { column, data } = columnData;
    return column.isNone(idx) ? null : column.type === DG.COLUMN_TYPE.DATE_TIME ?
      new Date(data[idx] * 1e-3) : data[idx];
  }

  isSameDate(x: timePoint, y: timePoint): boolean {
    if (x instanceof Date && y instanceof Date)
      return x.getTime() === y.getTime();
    else if ((typeof x === typeof y && typeof x === 'number') || (x == null || y == null))
      return x === y;
    grok.shell.warning('The columns of different types cannot be used for representing dates.');
    return false;
  }

  getColumnMin(column: DG.Column): number {
    return column.type === DG.COLUMN_TYPE.DATE_TIME ? new Date(column.min * 1e-3).getTime() : column.min;
  }

  getColumnMax(column: DG.Column): number {
    return column.type === DG.COLUMN_TYPE.DATE_TIME ? new Date(column.max * 1e-3).getTime() : column.max;
  }

  getSeriesData() {
    this.data.length = 0;
    const tempObj: Indexable = {};

    const colorCategories = this.columnData.colorByColumnName?.categories;
    const colorBuf = this.columnData.colorByColumnName?.data;
    const eventCategories = this.columnData.eventColumnName?.categories;
    const eventBuf = this.columnData.eventColumnName?.data;
    const startColumn = this.columnData.startColumnName?.column;
    const endColumn = this.columnData.endColumnName?.column;
    const eventColumns = this.columnData.eventsColumnNames ? Object.keys(this.columnData.eventsColumnNames)
      .map((columnName) => this.columnData.eventsColumnNames![columnName]!.column) : [];

    if ([startColumn, endColumn, ...eventColumns].some((column) =>
      column ? column.type !== DG.COLUMN_TYPE.DATE_TIME : false))
      this.removeTimeOptions();
    else
      this.addTimeOptions();

    for (const i of this.dataFrame.filter.getSelectedIndexes()) {
      const id = this.getStrValue(this.columnData.splitByColumnName!, i);
      const color = this.colorByColumnName ? colorCategories![colorBuf![i]] : this.defaultColor;
      let start = startColumn ? this.getSafeValue(this.columnData.startColumnName!, i) : null;
      let end = endColumn ? this.getSafeValue(this.columnData.endColumnName!, i) : null;
      if (start === end && end === null && (!this.eventsColumnNames || this.eventsColumnNames.length === 0)) continue;
      if (this.eventsColumnNames?.length > 0) {
        for (const columnName of Object.keys(this.columnData.eventsColumnNames!)) {
          const start = this.getSafeValue(this.columnData.eventsColumnNames![columnName], i);
          const end = this.showOpenIntervals ? this.dataMax : null;
          const key = `${id}-${columnName}-${start}-${end}`;
          if (key in tempObj) {
            tempObj[key][3].push(color);
            tempObj[key][4].push(columnName);
          } else
            tempObj[key] = [id, start, end, [color], [columnName], this.dataFrame.selection.get(i)];
        }
      }
      if (startColumn && endColumn) {
        if (this.showOpenIntervals) {
          // TODO: handle edge case of different column types
          if (start == null)
            start = Math.min(this.getColumnMin(startColumn), this.getColumnMin(endColumn));
          if (end == null)
            end = Math.max(this.getColumnMax(startColumn), this.getColumnMax(endColumn));
        }
        const event = this.eventColumnName ? eventCategories![eventBuf![i]] : null;
        const key = `${id}-${event}-${start}-${end}`;
        if (key in tempObj) {
          tempObj[key][3].push(color);
          tempObj[key][4].push(event);
        } else
          tempObj[key] = [id, start, end, [color], [event], this.dataFrame.selection.get(i)];
      }
    }

    this.data = Object.values(tempObj).sort((a, b) => {
      if (a[0] > b[0]) return 1;
      if (a[0] < b[0]) return -1;

      // Items with the same id are sorted based on `start` value
      if (a[1] > b[1]) return 1;
      if (a[1] < b[1]) return -1;
      return 0;
    });

    return this.data;
  }

  addTimeOptions(): void {
    if (this.dateFormat === null)
      (this.props as Indexable).dateFormat = this.defaultDateFormat;

    this.option.xAxis = {
      type: 'time',
      boundaryGap: ['5%', '5%'],
      axisLabel: { formatter: this.dateFormat },
    };
  }

  removeTimeOptions(): void {
    this.option.xAxis = {
      type: 'value',
      axisLabel: { formatter: null },
      min: 'dataMin',
      max: (value: { min: number; max: number; }) => this.dataMax = value.max,
    };
  }

  showErrorMessage(msg: string): void {
    this.titleDiv.innerText = msg;
    $(this.titleDiv).addClass('d4-viewer-error');
    $(this.legendDiv).hide();
    $(this.chart.getDom()).hide();
  }

  updateContainers(): void {
    $(this.titleDiv).removeClass().empty();
    if (this.colorByColumnName)
      this.switchLegendVisibility(this.legendVisibility);
    $(this.chart.getDom()).show();
  }

  render(): void {
    this.updateContainers();
    if (!this.splitByColumnName || ((!this.startColumnName && !this.endColumnName) &&
        (!this.eventsColumnNames || this.eventsColumnNames.length === 0))) {
      this.showErrorMessage('Not enough data to produce the result.');
      return;
    }
    this.option.series[0].data = this.getSeriesData();
    this.updateZoom();
    this.chart.setOption(this.option);
  }
}
