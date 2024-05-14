export const option: { [key: string]: any } = {
  animation: false,
  silent: false,
  legend: {
    show: true,
  },
  color: ['#b0d7ff', '#4fbcf7', '#4287cc'],
  radar: {
    name: {
      textStyle: {
        color: '#4d5261',
        backgroundColor: 'transparent',
        fontSize: '13px',
        fontFamily: 'Roboto',
        padding: [3, 5],
      },
    },
    radius: '60%',
    indicator: [],
  },
  tooltip: {
    show: false,
  },
  series: [{
    type: 'radar',
    data: [],
  }, {
    type: 'radar',
    data: [],
  }, {
    type: 'radar',
    data: [],
  }],
};

export const MAXIMUM_ROW_NUMBER = 1000;
export const MAXIMUM_SERIES_NUMBER = 25;
export const MAXIMUM_COLUMN_NUMBER = 10;
export interface RadarIndicator {
  name: string;
  max?: number;
  min?: number;
}

export const HIGHLIGHT_WIDTH = '2.5';
export const LINE_MIN_WIDTH = '0.5';
export const LINE_MAX_WIDTH = '2';
export const MOUSE_OVER_GROUP_COLOR = '#20CDCD';