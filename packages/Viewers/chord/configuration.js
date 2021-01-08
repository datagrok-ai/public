export const layoutConf = {
  innerRadius: 250,
  outerRadius: 300,
  cornerRadius: 0,
  gap: 0.04,
  labels: {
    display: true,
    position: 'center',
    size: '14px',
    color: '#000000',
    radialOffset: 35,
  },
  ticks: {
    display: false,
    color: 'grey',
    spacing: 10000000,
    labels: true,
    labelSpacing: 10,
    labelSuffix: '',
    labelDenominator: 1000000,
    labelDisplay0: false,
    labelSize: '10px',
    labelColor: '#000000',
    labelFont: 'default',
    majorSpacing: 5,
    size: {
      minor: 2,
      major: 5,
    }
  },
  events: {}
};

export const chordConf = {
  color: (datum, index) => (datum.source.id === datum.target.id) ? '#ff5500' : '#fd6a62',
  opacity: 0.7,
  radius: null,
  events: {}
};
