import * as DG from 'datagrok-api/dg';
import type {TimelinesViewer} from './timelines-viewer';
import {Box, Readings, categoryName, drawnItems} from '../../utils/utils';

/** The rows the plot draws: the lanes in the y axis' zoom window, whether or not an interval of theirs
 * falls inside the x window. */
function drawnLaneCount(chart: any): number {
  const axisModel = chart?.getModel?.()?.getComponent?.('yAxis', 0);
  const categories: unknown[] = axisModel?.getCategories?.() ?? [];
  const scale = axisModel?.axis?.scale;
  return (scale?.getTicks?.() ?? []).filter((tick: any) =>
    categories[scale.getRawOrdinalNumber(tick.value)] != null).length;
}

/** The category axis skips labels that would overlap, so only drawn lanes have one. AxisBuilder tags
 * a label `label_<ordinal>`, which is how its lane is found among the axis' categories. */
function laneLabels(chart: any): {name: string, box: Box}[] {
  const axisModel = chart?.getModel?.()?.getComponent?.('yAxis', 0);
  const categories: unknown[] = axisModel?.getCategories?.() ?? [];
  const view = axisModel == null ? null : chart.getViewOfComponentModel?.(axisModel);
  const labels: {name: string, box: Box}[] = [];
  view?.group?.traverse?.((el: any) => {
    const anid = typeof el.anid === 'string' && el.anid.startsWith('label_') ? Number(el.anid.slice(6)) : NaN;
    if (!Number.isInteger(anid) || categories[anid] == null || el.ignore || el.invisible)
      return;
    const rect = el.getBoundingRect().clone();
    const transform = el.getComputedTransform();
    if (transform)
      rect.applyTransform(transform);
    const box = {x: rect.x, y: rect.y, width: rect.width, height: rect.height};
    labels.push({name: categoryName(categories[anid]), box});
  });
  return labels;
}

/** While the viewer shows its message the chart is hidden, so nothing of it is reported. The series
 * data is what survived the legend filter in `render` and both dataZoom windows. */
export function timelinesStatus(v: TimelinesViewer): DG.IWidgetStatus & {values: Readings} {
  const error = v.renderError;
  const canvas = error === null ? v.chart.getDom().querySelector('canvas') : null;
  const model = (v.chart as any).getModel?.();
  const hitAreas: {[name: string]: Box} = {};
  const values: Readings = {};

  values['rows shown'] = v.filter.trueCount;
  values['lanes'] = canvas ? drawnLaneCount(v.chart) : 0;
  values['intervals'] = canvas ? drawnItems(model?.getSeriesByIndex?.(0)?.getData?.()).length : 0;

  if (canvas) {
    const grid = model?.getComponent?.('grid', 0)?.coordinateSystem?.getRect?.();
    if (grid)
      hitAreas['view'] = {x: grid.x, y: grid.y, width: grid.width, height: grid.height};
    for (const label of laneLabels(v.chart))
      hitAreas[`lane ${label.name}`] = label.box;
  }

  return {
    parts: canvas ? {root: v.root, canvas} : {root: v.root},
    hitAreas, values, shortcuts: {}, events: [], description: null, error,
  };
}
