/* eslint-disable max-len */
/**
 * Connections of a sequence drawn in its grid cell: the splitted sequence's graph info (BILN cyclizations, HELM bonds)
 * as brackets over the monomers they join, instead of the `C(1,3)` marks. Layout and drawing only; {@link MonomerPlacer}
 * decides when (a setting, a row with connections, room in the cell) and where.
 *
 * The drawing follows how disulfides are drawn over one-letter sequences: a square bracket with a leg on each monomer,
 * shorter brackets nearer the letters, overlapping ones on separate lanes; a chain's backbone ring under the chain;
 * a monomer bonded to three or more others (a scaffold on three Cys) as one comb.
 */
import {ISeqConnection, ISeqGraphInfoHelm, ISeqSplitted} from './macromolecule/types';

/** A connection between two flat positions of a row: positions across all chains, in the order they are drawn */
export interface IDrawnConnection {
  /** The flat positions, from < to */
  from: number;
  to: number;
  /** The R-groups at {@link from} and {@link to}: 1 the N side, 2 the C side, 3 and up a side chain */
  rFrom: number;
  rTo: number;
  /** The chains (0-based) of {@link from} and {@link to} */
  chainFrom: number;
  chainTo: number;
  /** One chain closed into a ring, its first monomer's R1 to its last monomer's R2: drawn under the chain */
  headToTail: boolean;
  /** A hydrogen bond (HELM `pair`): drawn dashed */
  hbond: boolean;
  /** Set by {@link assignLanes}: 1 above the monomers, -1 below */
  side: 1 | -1;
  /** Set by {@link assignLanes}: 0 is the lane nearest the monomers */
  lane: number;
  /** A hub's comb ({@link combineHubs}): one bar over the hub and its partners, standing for all of {@link conns}.
   * {@link from} and {@link to} then span the whole comb */
  comb?: { hub: number, members: number[], conns: IDrawnConnection[] };
}

/** Lane spacing and stroke of the brackets, in CSS px */
export interface IConnectionsGeometry {
  /** From the monomers to the first lane */
  tick0: number;
  /** Between lanes */
  laneGap: number;
  lineWidth: number;
}

/** Vertical metrics of the monomers' font, as offsets from the y a fillText with textBaseline 'top' gets */
export interface IConnectionsFontMetrics {
  /** The top of the capitals */
  capTopOff: number;
  /** The bottom of the capitals: the baseline */
  capBotOff: number;
  /** The em height (the font's ascent + descent) */
  emH: number;
  /** The part of the em under the baseline */
  descBelow: number;
}

/** How a row's connections fit its cell */
export interface IConnectionsLayout {
  /** The connections to draw, a hub's connections folded into one comb */
  conns: IDrawnConnection[];
  /** Lanes above and below the monomers (in a stacked layout: above the first line and below the last) */
  above: number;
  below: number;
  geo: IConnectionsGeometry;
  /** Stacked: the chains' first flat positions, one line each; empty when the row is drawn on one line */
  chains: number[];
}

/** Lines of bonds whose two ends are not the same monomer, and of rings */
export const CONNECTION_GREY = 'rgb(105,105,105)';

/** Colours of the bond numbers that stand for bonds crossing a line break (multiline) */
export const CONNECTION_NUMBER_COLORS = ['#4e79a7', '#e15759', '#f28e2b', '#59a14f', '#b07aa1', '#76b7b2', '#9c755f', '#ff9da7'];

/** The row's connections in flat positions, in graph info order; empty when it has none.
 * {@link ISeqConnection.monomerIndex1} counts within its chain: flat = chain start + monomer index. */
export function getDrawnConnections(seqSS: ISeqSplitted): IDrawnConnection[] {
  const gi = seqSS.graphInfo;
  if (!gi) return [];
  const hbonds = (gi as ISeqGraphInfoHelm).hbonds;
  if (!gi.connections?.length && !hbonds?.length) return [];
  const starts = gi.disjointSeqStarts?.length ? gi.disjointSeqStarts : [0];
  const len = seqSS.length;
  const chainLast = (c: number): number => (c + 1 < starts.length ? starts[c + 1] : len) - 1;
  /** the chain a flat position is in (the connection's own chain indices are not trusted) */
  const chainOf = (p: number): number => {
    let c = starts.length - 1;
    while (c > 0 && starts[c] > p) c--;
    return c;
  };
  const res: IDrawnConnection[] = [];
  const add = (c: ISeqConnection, hbond: boolean): void => {
    const p1 = (starts[c.seqIndex1] ?? NaN) + c.monomerIndex1;
    const p2 = (starts[c.seqIndex2] ?? NaN) + c.monomerIndex2;
    if (!(p1 >= 0 && p2 >= 0 && p1 < len && p2 < len) || p1 === p2) return;
    const fwd = p1 < p2;
    const chainFrom = chainOf(fwd ? p1 : p2);
    const chainTo = chainOf(fwd ? p2 : p1);
    const from = fwd ? p1 : p2;
    const to = fwd ? p2 : p1;
    const rFrom = fwd ? c.rGroup1 : c.rGroup2;
    const rTo = fwd ? c.rGroup2 : c.rGroup1;
    const headToTail = !hbond && chainFrom === chainTo && rFrom === 1 && rTo === 2 &&
      from === starts[chainFrom] && to === chainLast(chainFrom);
    res.push({from, to, rFrom, rTo, chainFrom, chainTo, headToTail, hbond, side: 1, lane: 0});
  };
  for (const c of gi.connections ?? []) add(c, false);
  for (const c of hbonds ?? []) add(c, true);
  return res;
}

/** A monomer bonded to three or more others through side bonds (a TATA/THA scaffold on three Cys, a branching Lys)
 * becomes one comb: a bar over it and its partners, a leg to each and a dot on the hub. One lane instead of one per
 * bond, and the picture of one object holding several. */
export function combineHubs(conns: IDrawnConnection[]): IDrawnConnection[] {
  const byPos = new Map<number, IDrawnConnection[]>();
  for (const c of conns) {
    if (c.headToTail || c.hbond) continue;
    for (const p of [c.from, c.to]) {
      let l = byPos.get(p);
      if (!l) byPos.set(p, l = []);
      l.push(c);
    }
  }
  const used = new Set<IDrawnConnection>();
  const res: IDrawnConnection[] = [];
  for (const [hub, cs] of byPos) {
    const free = cs.filter((c) => !used.has(c));
    if (free.length < 3) continue;
    const members = free.map((c) => c.from === hub ? c.to : c.from);
    for (const c of free) used.add(c);
    res.push({...free[0], from: Math.min(hub, ...members), to: Math.max(hub, ...members), comb: {hub, members, conns: free}});
  }
  if (used.size === 0) return conns;
  for (const c of conns) if (!used.has(c)) res.push(c);
  return res;
}

const overlaps = (a: IDrawnConnection, b: IDrawnConnection): boolean => a.from <= b.to && b.from <= a.to;

/** Outer lanes first: {@link drawConnections} draws in the given order, so an inner bracket stays on top where they cross */
const byLaneOuterFirst = (a: IDrawnConnection, b: IDrawnConnection): number => b.lane - a.lane;

const SOLID: number[] = [];
const DASHED = [2, 2];

/** Sides and lanes: rings under their chain and the rest above it, unless {@link sideOf} says otherwise; shorter
 * connections nearer the monomers; two connections that overlap or share a monomer never share a lane, so a bracket
 * that holds another is always higher.
 * @return the numbers of lanes above and below */
export function assignLanes(conns: IDrawnConnection[], sideOf?: (c: IDrawnConnection) => 1 | -1): { above: number, below: number } {
  const order = conns.slice().sort((a, b) => (a.to - a.from) - (b.to - b.from) || a.from - b.from);
  const placed: IDrawnConnection[] = [];
  let above = 0;
  let below = 0;
  for (const c of order) {
    c.side = sideOf ? sideOf(c) : c.headToTail ? -1 : 1;
    let lane = 0;
    while (placed.some((p) => p.side === c.side && p.lane === lane && overlaps(p, c))) lane++;
    c.lane = lane;
    placed.push(c);
    if (c.side === 1) above = Math.max(above, lane + 1);
    else below = Math.max(below, lane + 1);
  }
  return {above, below};
}

/** Lanes 0.3 em apart with a stroke of 1/12 em, or packed closer when the cell is short */
export function connectionsGeometry(fontSize: number, compressed: boolean): IConnectionsGeometry {
  return compressed ?
    {tick0: 2, laneGap: 2.5, lineWidth: Math.max(1, fontSize / 14)} :
    {tick0: Math.max(3, Math.round(fontSize * 0.3)), laneGap: Math.max(3, Math.round(fontSize * 0.3)), lineWidth: Math.max(1, fontSize / 12)};
}

/** The height the lanes of one side take */
export function sideHeight(lanes: number, geo: IConnectionsGeometry): number {
  return lanes > 0 ? geo.tick0 + (lanes - 1) * geo.laneGap + geo.lineWidth + 1 : 0;
}

/** The height of one line of monomers with its lanes: the capitals, the lanes above, and below either the lanes or the
 * descenders */
export function lineHeightWithLanes(above: number, below: number, geo: IConnectionsGeometry, met: IConnectionsFontMetrics): number {
  return sideHeight(above, geo) + (met.capBotOff - met.capTopOff) + Math.max(sideHeight(below, geo) + 2, met.descBelow);
}

/** The least gap between two stacked lines: room for a clearly visible connection between them (21 px at 12 px font,
 * so two chains go on two lines in rows of about 45 px and more, three in about 80 px) */
export function stackedLineGap(fontSize: number): number {
  return Math.round(fontSize * 1.75);
}

const fontMetricsCache = new Map<string, IConnectionsFontMetrics>();
let measureContext: CanvasRenderingContext2D | null = null;

/** The vertical metrics of {@link font}, measured once per font */
export function measureConnectionsFont(font: string, fontSize: number): IConnectionsFontMetrics {
  let res = fontMetricsCache.get(font);
  if (res) return res;
  try {
    measureContext ??= document.createElement('canvas').getContext('2d');
  } catch {
    // no DOM: the estimate below
  }
  if (measureContext) {
    measureContext.font = font;
    measureContext.textBaseline = 'top';
    const m = measureContext.measureText('M');
    const emH = m.fontBoundingBoxAscent + m.fontBoundingBoxDescent;
    res = {capTopOff: -m.actualBoundingBoxAscent, capBotOff: m.actualBoundingBoxDescent, emH,
      descBelow: emH - m.actualBoundingBoxDescent};
  } else
    res = {capTopOff: 0.25 * fontSize, capBotOff: 0.83 * fontSize, emH: 1.17 * fontSize, descBelow: 0.34 * fontSize};
  fontMetricsCache.set(font, res);
  return res;
}

/** One line per chain reads when every inter-chain bond joins neighbouring lines through a side chain (or pairs them),
 * and only the first and the last chain have bonds of their own (drawn above the first line and below the last);
 * every chain has monomers. Chains joined backbone to backbone (one chain's C-terminus to the next one's N-terminus)
 * are one backbone, which a single line shows better */
function isStackable(conns: IDrawnConnection[], chains: number[], len: number): boolean {
  if (chains.length < 2) return false;
  for (let k = 0; k < chains.length; k++)
    if ((k + 1 < chains.length ? chains[k + 1] : len) <= chains[k]) return false;
  const last = chains.length - 1;
  return conns.every((c) => c.chainFrom !== c.chainTo ?
    Math.abs(c.chainTo - c.chainFrom) === 1 && (c.hbond || c.rFrom > 2 || c.rTo > 2) :
    c.chainFrom === 0 || c.chainFrom === last);
}

/** How the row's connections fit a cell of height {@link h}, or null when it has none or they do not fit: the cell
 * then prints the monomers as it always did.
 * @param stackedAllowed one line per chain may be used (not in a scrolled sequence) */
export function layoutConnections(seqSS: ISeqSplitted, h: number, fontSize: number, met: IConnectionsFontMetrics,
  stackedAllowed: boolean,
): IConnectionsLayout | null {
  const all = getDrawnConnections(seqSS);
  if (all.length === 0) return null;
  const capH = met.capBotOff - met.capTopOff;
  const chains = seqSS.graphInfo?.disjointSeqStarts ?? [0];
  if (stackedAllowed && isStackable(all, chains, seqSS.length)) {
    const own = all.filter((c) => c.chainFrom === c.chainTo);
    const lanes = assignLanes(own, (c) => c.chainFrom === 0 ? 1 : -1);
    const geo = connectionsGeometry(fontSize, false);
    const need = sideHeight(lanes.above, geo) + chains.length * capH + (chains.length - 1) * stackedLineGap(fontSize) +
      Math.max(sideHeight(lanes.below, geo) + 2, met.descBelow) + 2;
    if (h >= need) {
      const inter = all.filter((c) => c.chainFrom !== c.chainTo);
      return {conns: own.sort(byLaneOuterFirst).concat(inter), ...lanes, geo, chains: chains.slice()};
    }
  }
  const conns = combineHubs(all);
  const lanes = assignLanes(conns);
  for (const compressed of [false, true]) {
    const geo = connectionsGeometry(fontSize, compressed);
    if (h >= lineHeightWithLanes(lanes.above, lanes.below, geo, met) + 2)
      return {conns: conns.slice().sort(byLaneOuterFirst), ...lanes, geo, chains: []};
  }
  return null;
}

/** Draws brackets (and combs) in the given order, outer lanes first ({@link layoutConnections} sorts them so).
 * {@link xOf} gives a position's x for an R-group: R1 at the monomer's left edge, R2 at its right edge, a side chain at
 * its centre; positions out of view lie past the cell's clip, so their brackets run off its edge. Brackets of the
 * hovered position are drawn thicker and the others fade. */
export function drawConnections(g: CanvasRenderingContext2D, conns: IDrawnConnection[], xOf: (pos: number, r?: number) => number,
  aboveY0: number, belowY0: number, geo: IConnectionsGeometry, colorOf: (c: IDrawnConnection) => string,
  hoverPos: number | null = null,
): void {
  const touches = (c: IDrawnConnection): boolean => hoverPos !== null && (c.comb ?
    c.comb.hub === hoverPos || c.comb.members.includes(hoverPos) : c.from === hoverPos || c.to === hoverPos);
  const anyHover = hoverPos !== null && conns.some(touches);
  g.save();
  try {
    g.lineCap = 'round';
    g.lineJoin = 'round';
    for (const c of conns) {
      const s = c.side;
      const y0 = s === 1 ? aboveY0 : belowY0;
      const laneY = y0 - s * (geo.tick0 + c.lane * geo.laneGap);
      const color = colorOf(c);
      const hovered = touches(c);
      g.strokeStyle = color;
      g.globalAlpha = anyHover && !hovered ? 0.3 : 1;
      g.lineWidth = hovered ? geo.lineWidth * 2 : geo.lineWidth;
      g.setLineDash(c.hbond ? DASHED : SOLID);
      g.beginPath();
      if (c.comb) {
        const xs = [c.comb.hub, ...c.comb.members].map((p) => xOf(p));
        g.moveTo(Math.min(...xs), laneY);
        g.lineTo(Math.max(...xs), laneY);
        for (const x of xs) {
          g.moveTo(x, laneY);
          g.lineTo(x, y0);
        }
        g.stroke();
        g.fillStyle = color;
        g.beginPath();
        g.arc(xOf(c.comb.hub), laneY, geo.lineWidth * 2, 0, Math.PI * 2);
        g.fill();
        continue;
      }
      const x1 = xOf(c.from, c.rFrom);
      const x2 = xOf(c.to, c.rTo);
      const r = Math.min(2.5, Math.abs(laneY - y0), Math.abs(x2 - x1) / 2);
      g.moveTo(x1, y0);
      g.lineTo(x1, laneY + s * r);
      g.arcTo(x1, laneY, x1 + r, laneY, r);
      g.lineTo(x2 - r, laneY);
      g.arcTo(x2, laneY, x2, laneY + s * r, r);
      g.lineTo(x2, y0);
      g.stroke();
    }
  } finally {
    g.restore();
  }
}

/** Draws small bond numbers centred over monomers, each number in its own colour: the bonds a line break cuts
 * (multiline). {@link items}: per monomer, its centre x, the y of the numbers' bottom and the numbers. */
export function drawBondNumbers(g: CanvasRenderingContext2D, items: { x: number, y: number, numbers: number[] }[],
  fontSize: number,
): void {
  g.save();
  try {
    g.font = `bold ${bondNumberFontSize(fontSize)}px sans-serif`;
    g.textBaseline = 'bottom';
    g.textAlign = 'left';
    for (const it of items) {
      const parts = it.numbers.map((n, i) => ({text: (i > 0 ? ',' : '') + n, n}));
      let x = it.x - parts.reduce((a, p) => a + g.measureText(p.text).width, 0) / 2;
      for (const p of parts) {
        g.fillStyle = CONNECTION_NUMBER_COLORS[(p.n - 1) % CONNECTION_NUMBER_COLORS.length];
        g.fillText(p.text, x, it.y);
        x += g.measureText(p.text).width;
      }
    }
  } finally {
    g.restore();
  }
}

/** The font size of {@link drawBondNumbers} */
export function bondNumberFontSize(fontSize: number): number {
  return Math.max(7, Math.round(fontSize * 0.62));
}

/** {@link css} darkened by {@link k} (0..1): a line in a monomer's colour, with more contrast than its letters */
export function darkerColor(css: string, k: number = 0.3): string {
  let r: number;
  let g: number;
  let b: number;
  const m = /^rgba?\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)/i.exec(css);
  if (m)
    [r, g, b] = [+m[1], +m[2], +m[3]];
  else if (/^#[0-9a-f]{6}$/i.test(css))
    [r, g, b] = [1, 3, 5].map((i) => parseInt(css.substring(i, i + 2), 16));
  else
    return css;
  const d = (v: number): number => Math.round(v * (1 - k));
  return `rgb(${d(r)},${d(g)},${d(b)})`;
}
