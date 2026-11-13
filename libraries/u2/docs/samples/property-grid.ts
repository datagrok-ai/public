import {Control, PropertyGrid, component, computed, span} from '@datagrok-libraries/u2';
import './property-grid.css';

export function chartSettings(): Control {
  return component(() => {
    const grid = new PropertyGrid();
    grid.setProperties([
      {name: 'title', type: 'string', description: 'Shown above the chart'},
      {name: 'markerSize', type: 'int', min: 1, max: 20, category: 'Markers'},
      {name: 'markerType', type: 'choice', choices: ['circle', 'square', 'cross'], category: 'Markers'},
      {name: 'opacity', type: 'double', min: 0, max: 1, category: 'Markers'},
      {name: 'showLegend', type: 'bool', category: 'Legend'},
      {name: 'legendPosition', type: 'choice', choices: ['Left', 'Right', 'Top', 'Bottom'], category: 'Legend'},
    ], {title: 'Potency vs. MW', markerSize: 6, markerType: 'circle', opacity: 0.8, showLegend: true, legendPosition: 'Right'});
    const v = grid.values;
    return [grid, span(computed(() => `"${v.value.title}": ${v.value.markerType} markers of size ${v.value.markerSize}, ` +
      `legend ${v.value.showLegend ? `on the ${String(v.value.legendPosition).toLowerCase()}` : 'hidden'}`), 'chart-settings-summary')];
  }, {name: 'chart-settings'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = chartSettings();
  host.append(demo.root);
  return demo;
}
