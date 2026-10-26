import * as DG from 'datagrok-api/dg';
import type {BandChartViewer} from './band-chart-viewer';
import {Box, Readings} from '../../utils/utils';

/** While the viewer's error covers the chart nothing of the previous frame is reported. */
export function bandChartStatus(v: BandChartViewer): DG.IWidgetStatus & {values: Readings} {
  const error = v.renderError;
  const canvas = error === null ? v.chart.getDom().querySelector('canvas') : null;
  const hitAreas: {[name: string]: Box} = {};
  const values: Readings = {};

  values['message'] = error ?? '';
  values['bands'] = canvas ? v.bands.map((b) => b.name).join(', ') : '';
  values['points'] = canvas ? v.bands.reduce((n, b) => n + b.x.length, 0) : 0;
  values['reference lines'] = canvas ? v.refLines.map((r) => `${r.vertical ? 'x' : 'y'} = ${r.value}`).join(', ') : '';

  if (canvas) {
    const box = canvas.getBoundingClientRect();
    hitAreas['view'] = {x: 0, y: 0, width: box.width, height: box.height};
  }

  return {
    parts: canvas ? {root: v.root, canvas} : {root: v.root},
    hitAreas, values, shortcuts: {}, events: [], description: null, error,
  };
}
