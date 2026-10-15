import * as DG from 'datagrok-api/dg';
import type {SankeyViewer} from './sankey';
import {Box, Readings, categoryName, svgHits, svgPoint} from '../../utils/utils';

/** The side of the square a link's hit area reports, centred on a point of its path. */
const HIT = 4;

/** Where along a link, by its length, a click target is looked for — the middle first. */
const ALONG = [0.5, 0.3, 0.7, 0.15, 0.85];

const linkName = (d: any): string => `${categoryName(d?.source?.name)} -> ${categoryName(d?.target?.name)}`;

/** A point of [path] that a click lands on — the path itself or another link of the same pair, which
 * selects the same rows. Links cross, so the middle of one can be under another. */
function linkPoint(svg: SVGSVGElement, path: SVGPathElement): {x: number, y: number} | null {
  const name = linkName((path as any).__data__);
  const length = path.getTotalLength();
  for (const t of ALONG) {
    const local = path.getPointAtLength(length * t);
    const p = svgPoint(path, local.x, local.y);
    if (svgHits(svg, p, (hit) => hit.classList.contains('link') && linkName((hit as any).__data__) === name))
      return p;
  }
  return null;
}

/** Everything is read from the svg the last render built, so a render that drew nothing — no rows, a
 * message — reports nothing. */
export function sankeyStatus(v: SankeyViewer): DG.IWidgetStatus & {values: Readings} {
  const error = v.renderError;
  const svg = error === null ? v.root.querySelector('svg') : null;
  const hitAreas: {[name: string]: Box} = {};
  const values: Readings = {};

  const origin = svg?.getBoundingClientRect();
  const nodes: string[] = [];
  for (const rect of Array.from(svg?.querySelectorAll('g.node > rect') ?? [])) {
    const box = rect.getBoundingClientRect();
    if (box.width === 0 || box.height === 0)
      continue;
    const name = categoryName((rect as any).__data__?.name);
    nodes.push(name);
    hitAreas[`node ${name}`] = {
      x: box.left - origin!.left, y: box.top - origin!.top, width: box.width, height: box.height,
    };
  }

  const links = Array.from(svg?.querySelectorAll<SVGPathElement>('path.link') ?? []);
  // every row draws its own link; a pair is placed by its topmost one, the last drawn
  const tried = new Set<string>();
  for (const path of links.slice().reverse()) {
    const name = `link ${linkName((path as any).__data__)}`;
    if (tried.has(name))
      continue;
    tried.add(name);
    const p = linkPoint(svg!, path);
    if (p !== null)
      hitAreas[name] = {x: p.x - HIT / 2, y: p.y - HIT / 2, width: HIT, height: HIT};
  }

  values['nodes'] = nodes.length;
  values['node names'] = nodes.join(', ');
  values['links'] = links.length;
  values['rows shown'] = v.filter.trueCount;

  return {
    parts: {root: svg ?? v.root},
    hitAreas, values, shortcuts: {}, events: [], description: null, error,
  };
}
