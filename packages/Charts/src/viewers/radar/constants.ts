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
        //borderRadius: 3,
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

export const MAXIMUM_SERIES_NUMBER = 25;