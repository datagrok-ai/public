import * as echarts from 'echarts';
import { EChartViewer } from './echart-viewer';

export class TimelinesViewer extends EChartViewer {
  constructor() {
    super();
    
    this.subjectColumnName = this.string('subjectColumnName', 'USUBJID');
    this.startColumnName = this.string('startColumnName', 'AESTDY');
    this.endColumnName = this.string('endColumnName', 'AEENDY');
    this.colorByColumnName = this.string('colorByColumnName', 'EVENT');

    this.marker = this.string('marker', 'circle', { choices: ['circle', 'rect'] });
    this.markerSize = this.int('markerSize', 6);
    this.markerPosition = this.string('markerPosition', 'main line',
      { choices: ['main line', 'above main line', 'scatter'] });
    this.lineWidth = this.int('lineWidth', 3);
    this.dateFormat = this.string('dateFormat', null, { choices: [
      '{yyyy}-{MM}-{dd}', '{M}/{d}/{yyyy}', '{MMM} {d}', '{dd}', '{d}'
    ]});
    this.axisPointer = this.string('axisPointer', 'shadow',
      { choices: ['cross', 'line', 'shadow', 'none'] });

    this.data = [];
    this.count = 0;
    this.selectionColor = '#FB8C28';
    this.zoomState = [[0, 100], [0, 100], [0, 100], [0, 100]];
    this.tooltipOffset = 10;
    this.initialized = false;

    this.option = {
      tooltip: {
        trigger: 'axis',
        showContent: false,
        axisPointer: { type: this.axisPointer }
      },
      grid: {
        left: '3%',
        right: '4%',
        top: '2%',
        bottom: '3%',
        containLabel: true
      },
      xAxis: {
        type: 'value',
        min: value => value.min - 1,
        max: value => value.max + 1
      },
      yAxis: {
        type: 'category',
        triggerEvent: true,
        axisTick: { show: false },
        axisLine: { show: false },
      },
      dataZoom: [
      {
        type: 'inside',
        xAxisIndex: [1, 2],
        start: this.zoomState[0][0],
        end: this.zoomState[0][1],
        filterMode: 'weakFilter',
      },
      {
        start: this.zoomState[1][0],
        end: this.zoomState[1][1],
        height: 10,
        bottom: '1%',
        filterMode: 'weakFilter',
      },
      { 
        type: 'inside',
        yAxisIndex: 0,
        start: this.zoomState[2][0],
        end: this.zoomState[2][1],
      },
      { 
        type: 'slider',
        yAxisIndex: 0,
        start: this.zoomState[3][0],
        end: this.zoomState[3][1],
        width: 10,
      }
      ],
      series: [
        {
          type: 'custom',
          progressive: 0,   // Disable progressive rendering
          animation: false,
          encode: {
            x: [1, 2], 
            y: 0,
            tooltip: 3
          }
        }
      ]
    };
  
    this.chart.on('dataZoom', () => {
      this.chart.getOption().dataZoom.forEach((z, i) => {
        this.zoomState[i][0] = z.start;
        this.zoomState[i][1] = z.end;
      });
    });
  }

  init() {
    if (!this.initialized) {
      // FIXME: shouldn't impact the context menu of another TimelinesViewer instance (unsubscribe)
      grok.events.onContextMenu.subscribe(args => {
        if (args.args.context.type === 'TimelinesViewer') {
          args.args.menu.item('Reset View', () => {
            this.zoomState = [[0, 100], [0, 100], [0, 100], [0, 100]];
            this.render();
          });
        }
      });

      this.initialized = true;
    }
  }

  onPropertyChanged(property) {
    if (property.name === 'axisPointer') {
      this.option.tooltip.axisPointer.type = property.get();
    }
    this.render();
  }

  onTableAttached() {
    this.init();

    this.columns = this.dataFrame.columns.byNames([
      this.subjectColumnName, this.startColumnName,
      this.endColumnName, this.colorByColumnName
    ]).map((col, ind) => {
      if (col === null) {
        const propNames = ['subjectColumnName', 'startColumnName',
          'endColumnName', 'colorByColumnName'];
        const newColumn = this.dataFrame.columns.byIndex(ind);
        this[propNames[ind]] = newColumn.name;
        return newColumn;
      }
      return col;
    });

    this.colorMap = this.columns[3].categories.reduce((colorMap, c, i) => {
      colorMap[c] = DG.Color.toRgb(DG.Color.getCategoricalColor(i));
      return colorMap;
    }, {});

    this.option.series[0].renderItem = (params, api) => {
      let overlap = false;
      if (params.dataIndex > 0) {
        const prev = this.data[params.dataIndex - 1];
        if (this.data[params.dataIndex][0] === prev[0] &&
            prev[1] && prev[2] && prev[1] !== prev[2] &&
            api.value(1) <= prev[2]) {
              overlap = true;
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
          (Math.round(Math.random()) * 2 - 1)*(this.markerSize * 3));

        group.children.push({
          type: this.marker,
          shape: this.marker === 'circle' ? {
            cx: xPos(0),
            cy: yPos(0),
            r: this.markerSize / 2
          } : {
            x: xPos(this.markerSize / 2),
            y: yPos(this.markerSize / 2),
            width: this.markerSize,
            height: this.markerSize
          },
          style: {
            fill: api.value(4) ? this.selectionColor : this.colorMap[isNaN(api.value(3)) ?
              this.data[params.dataIndex][3][0] : api.value(3)]
          }
        });
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
          let height = api.size([0, 1])[1];
          let offset = Math.max(this.markerSize * 2, this.lineWidth);
          // Shift along the Y axis
          rectShape.y += (this.count % 3) ? (this.count % 3 === 2) ?
            0 : offset-height/2 : height/2-offset;
          this.count += 1;
        }

        group.children.push({
          type: 'rect',
          transition: ['shape'],
          shape: rectShape,
          style: { fill: api.value(4) ? this.selectionColor : this.colorMap[isNaN(api.value(3)) ?
            this.data[params.dataIndex][3][0] : api.value(3)] }
        });
      }
  
      return group;
    }

    super.onTableAttached();
  }

  getSeriesData() {
    this.data.length = 0;
    let tempObj = {};

    let getTime = (i, j) => {
      if (this.columns[j].type === 'datetime') {
        if (this.dateFormat === null) {
          this.props.dateFormat = this.getProperty('dateFormat').choices[2];
        }
        this.option.xAxis = {
          type: 'time',
          boundaryGap: ['5%', '5%'],
          axisLabel: { formatter: this.dateFormat }
        };
      }
      return this.columns[j].type === 'datetime' ?
        new Date(`${this.columns[j].get(i)}`) : this.columns[j].isNone(i) ?
        null : this.columns[j].get(i);
    }

    let selectedIndexes = this.dataFrame.selection.getSelectedIndexes();

    for (let i of this.dataFrame.filter.getSelectedIndexes()) {
      let id = this.columns[0].get(i);
      let start = getTime(i, 1);
      let end = getTime(i, 2);
      if (start === end && end === null) continue;
      let event = this.columns[3].get(i);
      let key = `${id}-${start}-${end}`;
      if (tempObj.hasOwnProperty(key)) {
        tempObj[key][3].push(event);
      } else {
        tempObj[key] = [id, start, end, [event,], selectedIndexes.includes(i)];
      }
    }

    this.data = Object.values(tempObj).sort((a, b) => {
      if (a[0] > b[0]) return 1;
      if (a[0] < b[0]) return -1;

      // Items with the same id are sorted based on `start` value
      if (a[1] > b[1]) return 1;
      if (a[1] < b[1]) return -1;
    });

    this.chart.on('click', params => this.dataFrame.selection.handleClick( i => {
      if (params.componentType === 'yAxis') return this.columns[0].get(i) === params.value;
      if (params.componentType === 'series') {
        return params.value[0] === this.columns[0].get(i) &&
               params.value[1] === (this.columns[1].isNone(i) ? null : this.columns[1].get(i)) &&
               params.value[2] === (this.columns[2].isNone(i) ? null : this.columns[2].get(i))
      }
      return false;
    }, params.event.event));

    this.chart.on('mouseover', params => ui.tooltip.showRowGroup(this.dataFrame, i => {
      if (params.componentType === 'yAxis') return this.columns[0].get(i) === params.value;
      if (params.componentType === 'series') {
        return params.value[0] === this.columns[0].get(i) &&
               params.value[1] === (this.columns[1].isNone(i) ? null : this.columns[1].get(i)) &&
               params.value[2] === (this.columns[2].isNone(i) ? null : this.columns[2].get(i))
      }
      return false;
    }, params.event.event.x + this.tooltipOffset, params.event.event.y + this.tooltipOffset));

    this.chart.on('mouseout', () => ui.tooltip.hide());

    return this.data;
  }

  render() {
    this.option.series[0].data = this.getSeriesData();
    this.option.dataZoom.forEach((z, i) => {
      z.start = this.zoomState[i][0];
      z.end = this.zoomState[i][1];
    });
    this.chart.setOption(this.option);
  }
}
