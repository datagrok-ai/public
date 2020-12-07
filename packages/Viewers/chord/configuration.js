export const layoutConf = {
    innerRadius: 250,
    outerRadius: 300,
    cornerRadius: 10,
    gap: 0.04,
    labels: {
        display: true,
        position: 'center',
        size: '14px',
        color: '#000000',
        radialOffset: 20,
    },
    ticks: {
        display: true,
        color: 'grey',
        spacing: 10000000,
        labels: true,
        labelSpacing: 10,
        labelSuffix: '',
        labelDenominator: 1000000,
        labelDisplay0: true,
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
    color: '#fd6a62',
    opacity: 0.7,
    tooltipContent: d => {
        return `<h3>${d.source.id} ➤ ${d.target.id}: ${d.value}</h3><i>(CTRL+C to copy to clipboard)</i>`;
    },
    radius: null,
    events: {}
};
