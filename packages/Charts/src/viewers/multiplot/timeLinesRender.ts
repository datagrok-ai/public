// // initTimeLine
// function timeLinesD(params: any, api: any) {
//   const overlap = false;

//   const categoryIndex = api.value(categoryCol);
//   const start = api.coord([api.value(startCol), categoryIndex]);
//   const end = api.coord([api.value(endCol), categoryIndex]);
//   const width = end[0] - start[0];

//   const group = {
//     type: 'group',
//     children: [],
//   };

//   // if no reason to draw line
//   if (isNaN(api.value(startCol)) || isNaN(api.value(endCol)) || this.markerSize > width) {
//     const xPos = (shift) => isNaN(start[0]) ? end[0] : start[0] - shift;
//     const yPos = (shift) => end[1] - (this.markerPosition === 'main line' ? shift :
//       this.markerPosition === 'above main line' ? Math.max(this.markerSize, this.lineWidth) + shift :
//       ((params.dataIndex % 2) * 2 - 1)*(this.markerSize * 3));

//     const marker = {
//       type: this.marker,
//       shape: this.marker === 'circle' ? {
//         cx: xPos(0),
//         cy: yPos(0),
//         r: this.markerSize / 2,
//       } : this.marker === 'ring' ? {
//         cx: xPos(0),
//         cy: yPos(0),
//         r: this.markerSize / 2,
//         r0: this.markerSize / 4,
//       } : {
//         x: xPos(this.markerSize / 2),
//         y: yPos(this.markerSize / 2),
//         width: this.markerSize,
//         height: this.markerSize,
//       },
//       style: {
//         fill: 'red',
//         /*
//         fill: api.value(4) ? this.selectionColor : this.colorMap[isNaN(api.value(3)) ?
//           this.data[params.dataIndex][3][0] : api.value(3)]
//           */
//       },
//       x: 0,
//       y: 0,
//       rotation: 0,
//     };

//     if (this.marker === 'diamond') {
//       marker.type = 'rect';
//       marker.x = xPos(0);
//       marker.y = yPos(0);
//       marker.shape.x = -this.markerSize / 2;
//       marker.shape.y = -this.markerSize / 2;
//       marker.shape.r = this.markerSize / 4;
//       marker.rotation = 0.785398;
//     } else if (this.marker === 'rect') {
//       marker.x = 0;
//       marker.y = 0;
//       marker.shape.x = xPos(this.markerSize / 2);
//       marker.shape.y = yPos(this.markerSize / 2);
//       marker.shape.r = 0;
//       marker.rotation = 0;
//     }

//     group.children.push(marker);
//   } else {
//     // draw line
//     if (!this.echarts) {
//       debugger;
//     }
//     const rectShape = this.echarts.graphic.clipRectByRect({
//       x: start[0],
//       y: start[1] - this.lineWidth / 2,
//       width: width,
//       height: this.lineWidth,
//     //  type: 'rect'
//     }, {
//       x: params.coordSys.x,
//       y: params.coordSys.y,
//       width: params.coordSys.width,
//       height: params.coordSys.height,
//     });

//     if (overlap) {
//       const height = api.size([0, 1])[1];
//       const offset = Math.max(this.markerSize * 2, this.lineWidth);
//       // Shift along the Y axis
//       rectShape.y += (this.count % 3) ? (this.count % 3 === 2) ?
//         0 : offset-height/2 : height/2-offset;
//       this.count += 1;
//     }

//     group.children.push({
//       type: 'rect',
//       transition: ['shape'],
//       shape: rectShape,
//       // style: { fill: api.value(4) ? this.selectionColor : this.colorMap[isNaN(api.value(3)) ?
//       // this.data[params.dataIndex][3][0] : api.value(3)] }
//       style: {fill: 'green'},
//     });
//   }

//   return group;
// }

// export function initTimeLine(visibleIndex: number, viewerThis: any) {
//   const r = {
//     type: 'custom',
//     progressive: 0,
//     animation: false,
//     renderItem: timeLinesD.bind(viewerThis),
//     encode: {
//       x: [startCol, endCol],
//       y: categoryCol,
//       tooltip: 3,
//     },
//     data: this.timeLinesData[visibleIndex],
//   };
//   return r;
// }
