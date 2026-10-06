import * as DG from 'datagrok-api/dg';
import type {SunburstViewer} from './sunburst-viewer';

type Box = {x: number, y: number, width: number, height: number};

type Readings = {[name: string]: number | string | boolean};

/** The side of the square a segment's hit area reports, centred on a point the sector is known to
 * contain. A ring sector's bounding box is not a hit area — for a wide sector its centre is the
 * hole, not the ring. */
const HIT = 6;

type Segment = {path: string, value: number, point: {x: number, y: number}};

/** Every sector the last layout drew, addressed by its tree path the same way a click is
 * (`treePathInfo` minus the root, joined with ` | `), with a point the sector itself confirms it
 * contains. Reads echarts internals, so every step is optional: a sector we cannot name or cannot
 * place reports nothing rather than a guess. */
function laidOutSegments(chart: any): Segment[] {
  const data = chart?.getModel?.()?.getSeriesByIndex?.(0)?.getData?.();
  if (!data)
    return [];
  const segments: Segment[] = [];
  for (let i = 0; i < data.count(); i++) {
    const el = data.getItemGraphicEl(i);
    const shape = el?.shape;
    if (shape?.r == null)
      continue;

    const path = segmentPath(data, i);
    if (path == null || path === '')
      continue;

    const point = containedPoint(el, shape);
    if (point == null)
      continue;

    segments.push({path, value: Number(data.getValues?.(i)?.[0]), point});
  }
  return segments;
}

/** A segment's path from the root, root excluded — `Cancer | Male` — which is what the click
 * handler builds from `treePathInfo` and what the selection is keyed by. Falls back to the node's
 * own name when the tree is not walkable. */
function segmentPath(data: any, index: number): string | null {
  let node = data.tree?.getNodeByDataIndex?.(index);
  if (node == null)
    return data.getName?.(index) ?? null;
  const names: string[] = [];
  for (; node != null && node.parentNode != null; node = node.parentNode)
    names.unshift(String(node.name ?? ''));
  return names.join(' | ');
}

/** The mid-radius, mid-angle point of a sector, in the chart's own pixel space. zrender builds a
 * sector with `cy + r*sin(angle)`, so that candidate is offered first and is the one taken when the
 * element has no `contain` to ask; the other is the fallback, and a sector that accepts neither is
 * left unplaced. */
function containedPoint(el: any, shape: any): {x: number, y: number} | null {
  const r = ((shape.r0 ?? 0) + shape.r) / 2;
  const angle = (shape.startAngle + shape.endAngle) / 2;
  const x = shape.cx + r * Math.cos(angle);
  for (const y of [shape.cy + r * Math.sin(angle), shape.cy - r * Math.sin(angle)]) {
    if (el.contain?.(x, y) !== false)
      return {x, y};
  }
  return null;
}

/** How many sectors the layout has placed. `setOption` returns before the layout runs, so a canvas
 * is not evidence of a sunburst — this count is. */
export function laidOutSegmentCount(chart: any): number {
  return laidOutSegments(chart).length;
}

/** What the sunburst shows, for automation: one hit area per drawn sector, and the readings a test
 * compares. Replaces the radius-by-direction sweep a spec had to do to land on a ring. */
export function sunburstStatus(v: SunburstViewer): DG.IWidgetStatus & {values: Readings} {
  const error = v.renderError;
  const canvas = error === null ? v.root.querySelector('canvas') as HTMLCanvasElement | null : null;
  const hitAreas: {[name: string]: Box} = {};
  const values: Readings = {};

  values['hierarchy columns'] = (v.eligibleHierarchyNames ?? []).join(', ');
  values['on click'] = v.onClick;
  values['include nulls'] = v.includeNulls;
  values['rows shown'] = v.filter.trueCount;

  if (canvas) {
    const box = canvas.getBoundingClientRect();
    hitAreas['view'] = {x: 0, y: 0, width: box.width, height: box.height};
    const segments = laidOutSegments(v.chart);
    for (const segment of segments) {
      hitAreas[`segment ${segment.path}`] = {
        x: segment.point.x - HIT / 2, y: segment.point.y - HIT / 2, width: HIT, height: HIT,
      };
      values[`rows of segment ${segment.path}`] = segment.value;
    }
    values['segments'] = segments.length;
    values['segment names'] = segments.map((s) => s.path).join(', ');
  }

  return {
    parts: canvas ? {root: v.root, canvas} : {root: v.root},
    hitAreas, values, shortcuts: {}, events: [], description: null, error,
  };
}
