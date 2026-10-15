import * as DG from 'datagrok-api/dg';
import type {ChordViewer} from './chord-viewer';
import {Box, Readings, categoryName, svgHits, svgPoint} from '../../utils/utils';

/** The side of the square a category's or chord's hit area reports, centred on a point inside it. */
const HIT = 4;

/** How far in from the circle, as a share of the inner radius, a chord's end is probed. */
const INWARD = [0.97, 0.9, 0.8, 0.6];

/** A circos block: its arc spans [start, end] radians, clockwise from twelve o'clock, and its chords
 * end on its inner radius. */
type Block = {start: number, end: number, len: number, inner: number};

/** Where the angle [a] at radius [r] is in the layout's user space — d3's arc convention. */
const polar = (a: number, r: number) => ({x: r * Math.sin(a), y: -r * Math.cos(a)});

/** The radii the arc of a block fills along the angle [a], found by walking out from the centre: the
 * layout's radii are on a configuration object every chord viewer shares, the path is this one's. */
function arcRadii(path: SVGPathElement, a: number, limit: number): {inner: number, outer: number} | null {
  let inner = -1;
  let outer = -1;
  for (let r = 0; r <= limit; r++) {
    const {x, y} = polar(a, r);
    if (path.isPointInFill(new DOMPoint(x, y))) {
      if (inner < 0)
        inner = r;
      outer = r;
    }
  }
  return inner < 0 ? null : {inner, outer};
}

/** The first point of [candidates] (user space of [el]) that is inside [el] and that a click lands on
 * as [el] — or as something [same] takes for it. */
function clickPoint(svg: SVGSVGElement, el: SVGGeometryElement, candidates: {x: number, y: number}[],
  same: (hit: Element) => boolean): {x: number, y: number} | null {
  for (const c of candidates) {
    if (!el.isPointInFill(new DOMPoint(c.x, c.y)))
      continue;
    const p = svgPoint(el, c.x, c.y);
    if (svgHits(svg, p, (hit) => hit === el || same(hit)))
      return p;
  }
  return null;
}

/** Everything is read from the svg the last render built. */
export function chordStatus(v: ChordViewer): DG.IWidgetStatus & {values: Readings} {
  const error = v.renderError;
  const svg = error === null ? v.root.querySelector('svg') : null;
  const hitAreas: {[name: string]: Box} = {};
  const values: Readings = {};
  const box = (p: {x: number, y: number}): Box => ({x: p.x - HIT / 2, y: p.y - HIT / 2, width: HIT, height: HIT});

  const arcs = Array.from(svg?.querySelectorAll<SVGPathElement>('.cs-layout path') ?? []);
  const blocks = new Map<string, Block>();
  const limit = svg === null ? 0 : Math.max(svg.width.baseVal.value, svg.height.baseVal.value) / 2;
  // every arc of one layout shares its radii, so they are walked for once
  let radii: {inner: number, outer: number} | null = null;
  for (const path of arcs) {
    const d = (path as any).__data__;
    const span = d.end - d.start;
    const mid = d.start + span / 2;
    radii ??= arcRadii(path, mid, limit);
    if (radii === null)
      continue;
    blocks.set(d.id, {start: d.start, end: d.end, len: d.len, inner: radii.inner});
    const r = (radii.inner + radii.outer) / 2;
    const p = clickPoint(svg!, path, [mid, d.start + span / 4, d.end - span / 4].map((a) => polar(a, r)), () => false);
    if (p !== null)
      hitAreas[`category ${categoryName(d.label)}`] = box(p);
  }

  const chords = Array.from(svg?.querySelectorAll<SVGPathElement>('path.chord') ?? []);
  for (const path of chords) {
    const d = (path as any).__data__;
    const name = `chord ${categoryName(d.source.label)} -> ${categoryName(d.target.label)}`;
    const candidates: {x: number, y: number}[] = [];
    for (const end of [d.source, d.target]) {
      const block = blocks.get(end.id);
      if (block === undefined)
        continue;
      // the same arithmetic circos places a chord's end with
      const span = block.end - block.start;
      const a = block.start + (end.start + end.end) / 2 / block.len * span;
      for (const k of INWARD)
        candidates.push(polar(a, block.inner * k));
    }
    const p = clickPoint(svg!, path, candidates, (hit) => {
      const h = (hit as any).__data__;
      return hit.classList.contains('chord') &&
        h?.source?.label === d.source.label && h?.target?.label === d.target.label;
    });
    if (p !== null)
      hitAreas[name] = box(p);
  }

  values['categories'] = arcs.length;
  values['chords'] = chords.length;
  values['rows shown'] = v.filter.trueCount;

  return {
    parts: {root: svg ?? v.root},
    hitAreas, values, shortcuts: {}, events: [], description: null, error,
  };
}
