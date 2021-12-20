/** Initial options for `echarts` library. */
export const options = {
  tooltip: {
    trigger: 'axis',
    showContent: false,
    axisPointer: {},
  },
  grid: {
    left: '3%',
    right: '4%',
    top: '2%',
    bottom: '3%',
    containLabel: true,
  },
  animation: false,
  xAxis: {
    type: 'value',
    // min: (value) => value.min - 1,
    // max: (value) => value.max + 1,
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
      filterMode: 'weakFilter',
    },
    {
      type: 'slider',
      // xAxisIndex: [1, 2],
      height: 10,
      bottom: '1%',
      filterMode: 'weakFilter',
    },
    {
      type: 'inside',
      yAxisIndex: 0,
    },
    {
      type: 'slider',
      yAxisIndex: 0,
      width: 10,
    }
  ],
  series: [
    {
      type: 'custom',
      progressive: 0,   // Disable progressive rendering
      encode: { x: [1, 2], y: 0 },
    },
  ],
};

export function deepCopy(object) {
  return JSON.parse(JSON.stringify(object));
}
