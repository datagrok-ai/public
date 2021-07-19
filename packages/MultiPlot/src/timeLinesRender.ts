
export function initTimeLineSandBox(visibleIndex: number, viewerThis: any): any {
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

function timeLinesD(params, api) {
  const overlap = false;
  /*
  if (params.dataIndex > 0) {
    const prev = this.data[params.dataIndex - 1];
    const curSubj = this.data[params.dataIndex][0];
    if (curSubj === prev[0] &&
        prev[1] && prev[2] && prev[1] !== prev[2] &&
        api.value(1) <= prev[2]) {
          overlap = true;
          if (prevSubj !== curSubj) {
            this.count = 0;
            prevSubj = curSubj;
          }
        }
  }
*/

  const categoryIndex = api.value(0);
  const start = api.coord([api.value(1), categoryIndex]);
  const end = api.coord([api.value(2), categoryIndex]);
  const width = end[0] - start[0];

  const group = {
    type: 'group',
    children: [],
  };

  // if no reason to draw line
  if (isNaN(api.value(1)) || isNaN(api.value(2)) || this.markerSize > width) {
    const xPos = (shift) => isNaN(start[0]) ? end[0] : start[0] - shift;
    const yPos = (shift) => end[1] - (this.markerPosition === 'main line' ? shift :
      this.markerPosition === 'above main line' ? Math.max(this.markerSize, this.lineWidth) + shift :
      ((params.dataIndex % 2) * 2 - 1)*(this.markerSize * 3));

    const marker = {
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
        fill: 'red',
        /*
        fill: api.value(4) ? this.selectionColor : this.colorMap[isNaN(api.value(3)) ?
          this.data[params.dataIndex][3][0] : api.value(3)]
          */
      },
      x: 0,
      y: 0,
      rotation: 0,
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
    // draw line
    const rectShape = this.echarts.graphic.clipRectByRect({
      x: start[0],
      y: start[1] - this.lineWidth / 2,
      width: width,
      height: this.lineWidth,
    //  type: 'rect'
    }, {
      x: params.coordSys.x,
      y: params.coordSys.y,
      width: params.coordSys.width,
      height: params.coordSys.height,
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
      // style: { fill: api.value(4) ? this.selectionColor : this.colorMap[isNaN(api.value(3)) ?
      // this.data[params.dataIndex][3][0] : api.value(3)] }
      style: {fill: 'green'},
    });
  }

  return group;
}

export function initTimeLine(visibleIndex: number, viewerThis: any) {
  const r = {
    type: 'custom',
    progressive: 0,
    animation: false,
    renderItem: timeLinesD.bind(viewerThis),
    encode: {
      x: [1, 2],
      y: 0,
      tooltip: 3,
    },
    data: this.timeLinesData,
  };
  return r;
}
