import * as DG from 'datagrok-api/dg';
import type {RadarViewer} from './radar-viewer';
import {WARNING_CLASS} from './constants';
import {Box, MessageHandler, Readings, drawnItems} from '../../utils/utils';

/** echarts builds the indicator axes on the frame it draws, so these are the axes on screen, not the
 * option that is about to replace them. */
function drawnAxes(chart: any): string[] {
  const radar = chart?.getModel?.()?.getComponent?.('radar', 0)?.coordinateSystem;
  return (radar?.getIndicatorAxes?.() ?? []).map((axis: any) => String(axis.name));
}

/** `createSeriesData` draws the row lines without a symbol; the current and mouse-over rows are drawn
 * again on top with one, whether or not the legend let them through, and are not counted. */
function drawnRows(chart: any): number {
  const data = chart?.getModel?.()?.getSeriesByIndex?.(2)?.getData?.();
  return drawnItems(data).filter((i) => data.getItemVisual(i, 'symbol') === 'none').length;
}

/** While the viewer's error covers the chart nothing of the previous frame is reported. */
export function radarStatus(v: RadarViewer): DG.IWidgetStatus & {values: Readings} {
  const error = v.renderError;
  const canvas = error === null ? v.chart.getDom().querySelector('canvas') : null;
  const hitAreas: {[name: string]: Box} = {};
  const values: Readings = {};

  values['message'] = error ?? MessageHandler._getMessage(v.root, WARNING_CLASS) ?? '';
  values['axes'] = canvas ? drawnAxes(v.chart).join(', ') : '';
  values['rows shown'] = canvas ? drawnRows(v.chart) : 0;

  if (canvas) {
    const box = canvas.getBoundingClientRect();
    hitAreas['view'] = {x: 0, y: 0, width: box.width, height: box.height};
  }

  return {
    parts: canvas ? {root: v.root, canvas} : {root: v.root},
    hitAreas, values, shortcuts: {}, events: [], description: null, error,
  };
}
