import * as echarts from 'echarts';
import { EChartViewer } from '../echart-viewer';
import { options, deepCopy } from './echarts-options';


export class TimelinesViewer extends EChartViewer {
  constructor() {
    super();

    this.splitByColumnName = this.string('splitByColumnName');
    this.startColumnName = this.string('startColumnName');
    this.endColumnName = this.string('endColumnName');
    this.colorByColumnName = this.string('colorByColumnName');
    this.showOpenIntervals = this.bool('showOpenIntervals', false);
    this.eventColumnName = this.string('eventColumnName');
    this.showEventInTooltip = this.bool('showEventInTooltip', true);

    this.marker = this.string('marker', 'circle', { choices: ['circle', 'rect', 'ring', 'diamond'] });
    this.markerSize = this.int('markerSize', 6);
    this.markerPosition = this.string('markerPosition', 'main line',
      { choices: ['main line', 'above main line', 'scatter'] });
    this.lineWidth = this.int('lineWidth', 3);
    this.dateFormat = this.string('dateFormat', null, { choices: [
      '{yyyy}-{MM}-{dd}', '{M}/{d}/{yyyy}', '{MMM} {d}', '{dd}', '{d}'
    ]});
    this.axisPointer = this.string('axisPointer', 'shadow',
      { choices: ['cross', 'line', 'shadow', 'none'] });
    this.showZoomSliders = this.bool('showZoomSliders', true);

    this.splitByRegexps = [/^USUBJID$/, /id/i];
    this.colorByRegexps = [/^([A-Z]{2}(TERM|TEST|TRT|VAL)|(ACT)?ARM|MIDS(TYPE)?|VISIT)$/, /event/i];
    this.eventRegexps = [/^event$/, ...this.colorByRegexps];
    this.startRegexps = [/^((VISIT|[A-Z]{2}(ST)?)DY)$/, /start|begin/i];
    this.endRegexps = [/^((VISIT|[A-Z]{2}(EN)?)DY)$/, /stop|end/i];

    this.defaultDateFormat = '{MMM} {d}';
    this.data = [];
    this.count = 0;
    this.selectionColor = DG.Color.toRgb(DG.Color.selectedRows);
    this.zoomState = [[0, 100], [0, 100], [0, 100], [0, 100]];
    this.tooltipOffset = 10;
    this.initialized = false;
    this.option = deepCopy(options);
  }

  init() {
    if (!this.initialized) {
      this.updateZoom();
      this.chart.on('dataZoom', () => {
        this.chart.getOption().dataZoom.forEach((z, i) => {
          this.zoomState[i][0] = z.start;
          this.zoomState[i][1] = z.end;
          if(z.type === 'slider' && Object.keys(z).includes('yAxisIndex')){
            this.lineWidth = z.end - z.start < 60 ? z.end - z.start < 30 ? 3 : 2 : 1;
            this.markerSize = z.end - z.start < 60 ? z.end - z.start < 30 ? 6 : 4 : 3;
          }
        });
      });

      this.chart.on('rendered', () => {
        this.count = 0;
      });

      this.chart.on('click', params => this.dataFrame.selection.handleClick((i) => {
        if (params.componentType === 'yAxis') return this.getStrValue(this.columnData.splitByColumnName, i) === params.value;
        if (params.componentType === 'series') {
          return params.value[0] === this.getStrValue(this.columnData.splitByColumnName, i) &&
                 this.isSameDate(params.value[1], this.getSafeValue(this.columnData.startColumnName, i)) &&
                 this.isSameDate(params.value[2], this.getSafeValue(this.columnData.endColumnName, i));
        }
        return false;
      }, params.event.event));

      this.chart.on('mouseover', params => this.getTooltip(params));

      this.chart.on('mouseout', () => ui.tooltip.hide());

      this.option.tooltip.axisPointer.type = this.axisPointer;
      this.initialized = true;
    }
  }

  onPropertyChanged(property) {
    if (!this.initialized) return;
    if (property.name === 'axisPointer') {
      this.option.tooltip.axisPointer.type = property.get(this);
    } else if (property.name === 'showZoomSliders') {
      this.option.dataZoom.forEach((z) => {
        if (z.type === 'slider') z.show = this.showZoomSliders;
      });
    } else if (property.name.endsWith('ColumnName')) {
      const columnData = this.updateColumnData(property);
      if (property.name === 'colorByColumnName') {
        this.colorMap = this.getColorMap(columnData.categories);
      }
    }
    this.render();
  }

  getTooltip(params) {
    const x = params.event.event.x + this.tooltipOffset;
    const y = params.event.event.y + this.tooltipOffset;
    if (this.showEventInTooltip) {
      let tooltipContent = params.componentType === 'yAxis' ? ui.div(`${params.value}`) :
        ui.divV([ui.div(`key: ${params.value[0]}`),
        ui.div(`event: ${params.value[4]}`),
        ui.div(`start: ${params.value[1]}`),
        ui.div(`end: ${params.value[2]}`),
        ])
      ui.tooltip.show(tooltipContent, x, y);
    } else {
      ui.tooltip.showRowGroup(this.dataFrame, (i) => {
        if (params.componentType === 'yAxis') return this.getStrValue(this.columnData.splitByColumnName, i) === params.value;
        if (params.componentType === 'series') {
          return params.value[0] === this.getStrValue(this.columnData.splitByColumnName, i) &&
            this.isSameDate(params.value[1], this.getSafeValue(this.columnData.startColumnName, i)) &&
            this.isSameDate(params.value[2], this.getSafeValue(this.columnData.endColumnName, i));
        }
        return false;
      }, x, y);
    }
  }

  getColorMap(categories) {
    return categories.reduce((colorMap, c, i) => {
      colorMap[c] = DG.Color.toRgb(DG.Color.getCategoricalColor(i));
      return colorMap;
    }, {});
  }

  onTableAttached() {
    this.init();

    const columns = this.dataFrame.columns.toList();

    const strColumns = columns.filter(col => col.type === DG.COLUMN_TYPE.STRING)
      .sort((a, b) => a.categories.length - b.categories.length);

    const intColumns = columns.filter(col => col.type === DG.COLUMN_TYPE.INT)
      .sort((a, b) => a.stats.avg - b.stats.avg);

    this.splitByColumnName = (this.findColumn(columns, this.splitByRegexps) || strColumns[strColumns.length - 1]).name;
    this.startColumnName = (this.findColumn(columns, this.startRegexps, [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.DATE_TIME]) || intColumns[0]).name;
    this.endColumnName = (this.findColumn(columns, this.endRegexps, [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.DATE_TIME]) || intColumns[intColumns.length - 1]).name;
    this.colorByColumnName = (this.findColumn(columns, this.colorByRegexps, [DG.COLUMN_TYPE.STRING]) || strColumns[0]).name;
    this.eventColumnName = (this.findColumn(columns, this.eventRegexps, [DG.COLUMN_TYPE.STRING]) || strColumns[0]).name;

    const columnPropNames = ['splitByColumnName', 'startColumnName', 'endColumnName', 'colorByColumnName', 'eventColumnName'];
    const columnNames = [this.splitByColumnName, this.startColumnName, this.endColumnName, this.colorByColumnName, this.eventColumnName];

    this.columnData = columnPropNames.reduce((map, v, i) => {
      const column = this.dataFrame.getCol(columnNames[i]);
      map[v] = {
        column,
        data: column.getRawData(),
        categories: column.type === DG.COLUMN_TYPE.STRING ? column.categories : null,
      };
      return map;
    }, {});

    this.colorMap = this.getColorMap(this.dataFrame.getCol(this.colorByColumnName).categories);

    let prevSubj = null;

    this.option.series[0].renderItem = (params, api) => {
      let overlap = false;
      if (params.dataIndex > 0) {
        const prev = this.data[params.dataIndex - 1];
        const curSubj = this.data[params.dataIndex][0];
        if (curSubj === prev[0] &&
            prev[1] && prev[2] && prev[1] !== prev[2] &&
            api.value(1) < prev[2]) {
              overlap = true;
              if (prevSubj !== curSubj) {
                this.count = 0;
                prevSubj = curSubj;
              }
            }
      }

      const categoryIndex = api.value(0);
      const start = api.coord([api.value(1), categoryIndex]);
      const end = api.coord([api.value(2), categoryIndex]);
      const width = end[0] - start[0];

      let group = {
        type: 'group',
        children: []
      };

      if (isNaN(api.value(1)) || isNaN(api.value(2)) || this.markerSize > width) {
        const xPos = shift => isNaN(start[0]) ? end[0] : start[0] - shift;
        const yPos = shift => end[1] - (this.markerPosition === 'main line' ? shift :
          this.markerPosition === 'above main line' ? Math.max(this.markerSize, this.lineWidth) + shift :
          ((params.dataIndex % 2) * 2 - 1)*(this.markerSize * 3));

        let marker = {
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
            height: this.markerSize
          },
          style: {
            fill: api.value(5) ? this.selectionColor : this.colorMap[isNaN(api.value(3)) ?
              this.data[params.dataIndex][3][0] : api.value(3)]
          }
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

        group.children.push(marker);
      } else {
        const rectShape = echarts.graphic.clipRectByRect({
          x: start[0],
          y: start[1] - this.lineWidth / 2,
          width: width,
          height: this.lineWidth
        }, {
          x: params.coordSys.x,
          y: params.coordSys.y,
          width: params.coordSys.width,
          height: params.coordSys.height
        });

        if (overlap) {
          const height = api.size([0, 1])[1];
          const offset = Math.max(this.markerSize * 2, this.lineWidth);
          // Shift along the Y axis
          rectShape.y += (this.count % 3) ? (this.count % 3 === 2) ?
            0 : offset-height/2 : height/2-offset;
          this.count += 1;
        }

        group.children.push({
          type: 'rect',
          transition: ['shape'],
          shape: rectShape,
          style: { fill: api.value(5) ? this.selectionColor : this.colorMap[isNaN(api.value(3)) ?
            this.data[params.dataIndex][3][0] : api.value(3)] }
        });
      }

      return group;
    };

    super.onTableAttached();
  }

  /** Find a column based on provided name patterns and types.
   * The list of patterns should be sorted by priority in descending order. */
  findColumn(columns, regexps, types = []) {
    if (types.length) {
      columns = columns.filter((c) => types.includes(c.type));
    }
    for (let regex of regexps) {
      const column = columns.find((c) => c.name.match(regex));
      if (column) return column;
    }
    return null;
  }

  updateColumnData(prop) {
    const column = this.dataFrame.col(prop.get(this));
    if (column === null)
      return null;
    if (prop.name === 'startColumnName' || prop.name === 'endColumnName') {
      column.type === DG.COLUMN_TYPE.DATE_TIME ? this.addTimeOptions() : this.removeTimeOptions();
    }
    this.columnData[prop.name] = {
      column,
      data: column.getRawData(),
      categories: column.type === DG.COLUMN_TYPE.STRING ? column.categories : null,
    };
    return this.columnData[prop.name];
  }

  updateZoom() {
    this.option.dataZoom.forEach((z, i) => {
      z.start = this.zoomState[i][0];
      z.end = this.zoomState[i][1];
    });
  }

  getStrValue(columnData, idx) {
    const { column, categories, data } = columnData;
    return column.type === DG.COLUMN_TYPE.STRING ? categories[data[idx]] : column.getString(idx);
  }

  getSafeValue(columnData, idx) {
    const { column, data } = columnData;
    return column.isNone(idx) ? null : column.type === DG.COLUMN_TYPE.DATE_TIME ?
      new Date(data[idx] * 1e-3) : data[idx];
  }

  isSameDate(x, y) {
    if (x instanceof Date && y instanceof Date) {
      return x.getTime() === y.getTime();
    } else if (typeof x === typeof y && typeof x === 'number') {
      return x === y;
    }
    grok.shell.warning('The columns of different types cannot be used for representing dates.');
  }

  getSeriesData() {
    this.data.length = 0;
    let tempObj = {};

    const getColumnMin = (column) => column.type === DG.COLUMN_TYPE.DATE_TIME ? new Date(column.min * 1e-3) : column.min;
    const getColumnMax = (column) => column.type === DG.COLUMN_TYPE.DATE_TIME ? new Date(column.max * 1e-3) : column.max;

    const { categories: colorCategories, data: colorBuf } = this.columnData.colorByColumnName;
    const { categories: eventCategories, data: eventBuf } = this.columnData.eventColumnName;
    const { column: startColumn } = this.columnData.startColumnName;
    const { column: endColumn } = this.columnData.endColumnName;

    for (const i of this.dataFrame.filter.getSelectedIndexes()) {
      const id = this.getStrValue(this.columnData.splitByColumnName, i);
      let start = this.getSafeValue(this.columnData.startColumnName, i);
      let end = this.getSafeValue(this.columnData.endColumnName, i);
      if (start === end && end === null) continue;
      if (this.showOpenIntervals) {
        // TODO: handle edge case of different column types
        if (start == null)
          start = Math.min(getColumnMin(startColumn), getColumnMin(endColumn));
        if (end == null)
          end = Math.max(getColumnMax(startColumn), getColumnMax(endColumn));
      }
      const color = colorCategories[colorBuf[i]];
      const event = eventCategories[eventBuf[i]];
      const key = `${id}-${event}-${start}-${end}`;
      if (tempObj.hasOwnProperty(key)) {
        tempObj[key][3].push(color);
        tempObj[key][4].push(event);
      } else {
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

  addTimeOptions() {
    if (this.dateFormat === null) {
      this.props.dateFormat = this.defaultDateFormat;
    }
    this.option.xAxis = {
      type: 'time',
      boundaryGap: ['5%', '5%'],
      axisLabel: { formatter: this.dateFormat },
    };
  }

  removeTimeOptions() {
    this.props.dateFormat = null;
    this.option.xAxis = {
      type: 'value',
      boundaryGap: ['0%', '0%'],
      axisLabel: { formatter: null },
    };
  }

  showErrorMessage(msg) {
    this.root.appendChild(ui.divText(msg, 'd4-viewer-error'));
  }

  render() {
    if (this.splitByColumnName == null || !(this.startColumnName || this.endColumnName)) {
      this.showErrorMessage('Not enough data to produce the result.');
      return;
    }
    this.option.series[0].data = this.getSeriesData();
    this.updateZoom();
    this.chart.setOption(this.option);
  }
}
