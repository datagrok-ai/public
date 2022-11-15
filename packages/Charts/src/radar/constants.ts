export const option: { [key: string]: any } = {
    animation: false,
    silent: false,
    legend: {
      show: true
    },
    color: ['#b0d7ff', '#4fbcf7', '#4287cc'],
    radar: {
      name: {
        textStyle: {
          color: '#000000',
          backgroundColor: '#ffffff',
          borderRadius: 3,
          padding: [3, 5],
        },
      },
      indicator: [],
    },
    tooltip: {
      show: false
    },
    series: [{
      type: 'radar',
      data: [],
    },{
      type: 'radar',
      data: [],
    }, {
      type: 'radar',
      data: [],
    }],
  };
