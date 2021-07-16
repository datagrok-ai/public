
  export function initTimeLine(visibleIndex: number, viewerThis: any): any {
    let count = 0;
    const renderFailed = 1234;
    const timeLinesSeries = {
      type: 'custom',
      progressive: 0,
      animation: false,
      renderItem: ((t: any, echarts: any) => (params, echartAPI) => {
        const data = this.timeLinesData;
        if (!data || data.length == 0) return renderFailed;
        const customDebug = false;
        const av0 = echartAPI.value(0);
        const av1 = echartAPI.value(1);
        const av2 = echartAPI.value(2);
        const gridTopRaw = this.echartOptions.grid[visibleIndex].top;
        const gridTop = parseFloat(gridTopRaw);
        let overlap = false;
        if (params.dataIndex > 0 && data[params.dataIndex - 1][0] === data[params.dataIndex][0] &&
          echartAPI.value(1) <= data[params.dataIndex - 1][2]) {
          overlap = true;
        }
        const categoryIndex = echartAPI.value(0);
        const start = echartAPI.coord([echartAPI.value(1), categoryIndex + 0]);
        const end = echartAPI.coord([echartAPI.value(2), categoryIndex]);
        const height = echartAPI.size([0, 1])[1];
        const rect0 = {
          x: start[0],
          y: start[1] - this.timeLineWidth / 2,
          width: end[0] - start[0],
          height: this.timeLineWidth,
        };
        const rect1 = {
          x: params.coordSys.x,
          y: params.coordSys.y, // + gridTop,
          width: params.coordSys.width,
          height: params.coordSys.height,
        };
        let rectShape = echarts.graphic.clipRectByRect(rect0, rect1);
        if (!rectShape) {
          rectShape = {x: -220, y: 20, width: 100, height: 10};
        }
        let endColor = this.categoryColors[categoryIndex % this.categoryColors.length];
        if (this.selection.includes(params.dataIndex)) endColor = 'red';
        let group = {
          type: 'group',
          children: [{
            type: 'rect',
            transition: ['shape'],
            // shape: rectShape,
            shape: rectShape,
            style: {fill: echartAPI.value(3)},
          },
          {
            type: 'circle',
            shape: {
              cx: start[0], cy: end[1] + 0, r: this.timeLineRadius,
            },
            style: {fill: this.categoryColors[categoryIndex % this.categoryColors.length]},            
          },
          {
            type: 'circle',
            shape: {
              cx: end[0], cy: end[1] + 0, r: this.timeLineRadius,
            },
            style: {fill: endColor},
          },
          ],
        };
        if (overlap && (end[0] - start[0] > 20)) {
          // let shift = count ? 20 : -20;
          const shift = (count % 3) ? (count % 3 === 2) ?
            0 : this.timeLineOverlapShift - height / 2 : height / 2 - this.timeLineOverlapShift;
          //    shift = 0;
          count += 1;
          rectShape.y += shift;
          group = {
            type: 'group',
            children: [{
              type: 'rect',
              transition: ['shape'],
              shape: rectShape,
              style: {fill: echartAPI.value(3)},
            },
            {
              type: 'circle',
              shape: {
                cx: start[0], cy: end[1] + shift, r: this.timeLineRadius,
              },
              style: {fill: 'darkgreen'},
            },
            {
              type: 'circle',
              shape: {
                cx: end[0], cy: end[1] + shift, r: this.timeLineRadius,
              },
              style: {fill: 'red'},
            },
            ],
          };
        } // overlap
        return group;
      })(this, viewerThis.echarts).bind(viewerThis), // render item
      encode: {
        x: [1, 2],
        y: 0,
        tooltip: 3,
      },
      data: this.timeLinesData,
    }; // this.timeLineSeries
    return timeLinesSeries;
  } // initTimeLine