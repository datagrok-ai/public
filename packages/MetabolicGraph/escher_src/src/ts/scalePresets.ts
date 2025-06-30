const scalePresets: { [key: string]: { type: string, value?: number, color: string, size: number }[] } = {
  GaBuGeRd: [
    {type: 'min', color: '#c8c8c8', size: 12},
    {type: 'value', value: 0.01, color: '#9696ff', size: 16},
    {type: 'value', value: 20, color: '#209123', size: 20},
    {type: 'max', color: '#ff0000', size: 25}
  ],
  GaBuRd: [
    {type: 'min', color: '#c8c8c8', size: 12},
    {type: 'median', color: '#9696ff', size: 20},
    {type: 'max', color: '#ff0000', size: 25}
  ],
  RdYlBu: [
    {type: 'min', color: '#d7191c', size: 12},
    {type: 'median', color: '#ffffbf', size: 20},
    {type: 'max', color: '#2c7bb6', size: 25}
  ],
  GeGaRd: [
    {type: 'min', color: '#209123', size: 25},
    {type: 'value', value: 0, color: '#c8c8c8', size: 12},
    {type: 'max', color: '#ff0000', size: 25}
  ],
  WhYlRd: [
    {type: 'min', color: '#fffaf0', size: 20},
    {type: 'median', color: '#f1c470', size: 30},
    {type: 'max', color: '#800000', size: 40}
  ]
};

export default scalePresets;
