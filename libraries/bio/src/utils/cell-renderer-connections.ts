/* eslint-disable max-len */
/**
 * Connections of a sequence drawn in its grid cell: the splitted sequence's graph info (BILN cyclizations, HELM bonds)
 * as brackets over the monomers they join, instead of the `C(1,3)` marks. Layout and drawing only; {@link MonomerPlacer}
 * decides when (a setting, a row with connections, room in the cell) and where.
 *
 * The drawing follows how disulfides are drawn over one-letter sequences: a square bracket with a leg on each monomer,
 * shorter brackets nearer the letters, overlapping ones on separate lanes; a chain's backbone ring under the chain;
 * a monomer bonded to three or more others (a scaffold on three Cys) as one comb. Two chains bonded side by side (the
 * last monomer of one to the first of the next, as HELM writes a conjugate's CHEM parts) read as one line: a short
 * line between the two monomers in place of the chain break, its ends at the R-groups' heights.
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
  /** Two chains bonded side by side: {@link from} is the last monomer of its chain and {@link to} the first of the next
   * one. Drawn as a short line between them in place of the chain break ({@link drawInlineJoins}), not as a bracket;
   * the two chains are one line of a stacked layout. Not for a monomer bonded to three or more others (a comb) */
  inline: boolean;
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
  /** Stacked: the lines' first flat positions, a line being a chain or chains bonded side by side
   * ({@link IDrawnConnection.inline}); empty when the row is drawn on one line */
  chains: number[];
  /** The bonds drawn between two monomers side by side ({@link IDrawnConnection.inline}), not in {@link conns} */
  inline: IDrawnConnection[];
}

/** Lines of bonds whose two ends are not the same monomer, and of rings */
export const CONNECTION_GREY = 'rgb(105,105,105)';

/** Colours of the bond numbers that stand for bonds crossing a line break (multiline) */
export const CONNECTION_NUMBER_COLORS = ['#4e79a7', '#e15759', '#f28e2b', '#59a14f', '#b07aa1', '#76b7b2', '#9c755f', '#ff9da7'];

/** The index of the run that holds flat position {@link p}, of runs starting at {@link starts} (ascending) */
function runOf(starts: ArrayLike<number>, p: number): number {
  let k = starts.length - 1;
  while (k > 0 && starts[k] > p) k--;
  return k;
}

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
  const res: IDrawnConnection[] = [];
  const add = (c: ISeqConnection, hbond: boolean): void => {
    const p1 = (starts[c.seqIndex1] ?? NaN) + c.monomerIndex1;
    const p2 = (starts[c.seqIndex2] ?? NaN) + c.monomerIndex2;
    if (!(p1 >= 0 && p2 >= 0 && p1 < len && p2 < len) || p1 === p2) return;
    const fwd = p1 < p2;
    const from = fwd ? p1 : p2;
    const to = fwd ? p2 : p1;
    // the chain a flat position is in (the connection's own chain indices are not trusted)
    res.push({from, to, rFrom: fwd ? c.rGroup1 : c.rGroup2, rTo: fwd ? c.rGroup2 : c.rGroup1, chainFrom: runOf(starts, from),
      chainTo: runOf(starts, to), headToTail: false, hbond, inline: false, side: 1, lane: 0});
  };
  for (const c of gi.connections ?? []) add(c, false);
  for (const c of hbonds ?? []) add(c, true);
  // chains side by side: the last monomer of a chain bonded to the first of the next, one bond per chain break, unless
  // either monomer holds three or more bonds (a comb shows those better)
  let bonds: Map<number, number> | null = null;
  const joined = new Set<number>(); // the first positions of the chains joined to the chain before
  for (const c of res) {
    if (c.hbond || c.chainTo !== c.chainFrom + 1 || c.from !== chainLast(c.chainFrom) || c.to !== starts[c.chainTo] ||
      joined.has(c.to)) continue;
    if (!bonds) {
      bonds = new Map<number, number>();
      for (const b of res) if (!b.hbond) for (const p of [b.from, b.to]) bonds.set(p, (bonds.get(p) ?? 0) + 1);
    }
    if (bonds.get(c.from)! >= 3 || bonds.get(c.to)! >= 3) continue;
    c.inline = true;
    joined.add(c.to);
  }
  // a chain closed head to tail is a ring, its first monomer's R1 to its last monomer's R2, and so is a line of chains
  // side by side closed so
  const lines = joined.size > 0 ? starts.filter((s) => !joined.has(s)) : starts;
  const lineLast = (l: number): number => (l + 1 < lines.length ? lines[l + 1] : len) - 1;
  for (const c of res) {
    if (c.hbond || c.inline || c.rFrom !== 1 || c.rTo !== 2) continue;
    const l = runOf(lines, c.from);
    c.headToTail = c.chainFrom === c.chainTo && c.from === starts[c.chainFrom] && c.to === chainLast(c.chainFrom) ||
      joined.size > 0 && c.from === lines[l] && c.to === lineLast(l);
  }
  return res;
}

/** The first flat positions of the row's lines: its chains, where a chain bonded side by side to the one before it
 * ({@link IDrawnConnection.inline}) continues that one's line */
export function connectionLines(seqSS: ISeqSplitted, conns: IDrawnConnection[]): number[] {
  const starts = seqSS.graphInfo?.disjointSeqStarts?.length ? seqSS.graphInfo.disjointSeqStarts : [0];
  let joined: Set<number> | null = null;
  for (const c of conns) if (c.inline) (joined ??= new Set<number>()).add(c.to);
  return joined ? starts.filter((s) => !joined!.has(s)) : starts.slice();
}

/** The line ({@link connectionLines}) that holds flat position {@link p} */
export function lineOf(lines: number[], p: number): number {
  return runOf(lines, p);
}

/** A monomer bonded to three or more others through side bonds (a TATA/THA scaffold on three Cys, a branching Lys)
 * becomes one comb: a bar over it and its partners, a leg to each and a dot on the hub. One lane instead of one per
 * bond, and the picture of one object holding several. */
export function combineHubs(conns: IDrawnConnection[]): IDrawnConnection[] {
  const byPos = new Map<number, IDrawnConnection[]>();
  for (const c of conns) {
    if (c.headToTail || c.hbond || c.inline) continue;
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

/** One line per chain reads when every bond between lines joins neighbouring lines through a side chain (or pairs
 * them), and only the first and the last line have bonds of their own (drawn above the first line and below the last);
 * every line has monomers. A line is a chain, or chains bonded side by side ({@link connectionLines}): their bond is
 * drawn within the line. Chains joined backbone to backbone further along (one chain's C-terminus to the next one's
 * N-terminus) are one backbone, which a single line shows better */
function isStackable(conns: IDrawnConnection[], lines: number[], len: number): boolean {
  if (lines.length < 2) return false;
  for (let k = 0; k < lines.length; k++)
    if ((k + 1 < lines.length ? lines[k + 1] : len) <= lines[k]) return false;
  const last = lines.length - 1;
  return conns.every((c) => {
    if (c.inline) return true;
    const a = lineOf(lines, c.from);
    const b = lineOf(lines, c.to);
    return a !== b ? Math.abs(b - a) === 1 && (c.hbond || c.rFrom > 2 || c.rTo > 2) : a === 0 || a === last;
  });
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
  const inline = all.filter((c) => c.inline);
  const rest = inline.length > 0 ? all.filter((c) => !c.inline) : all;
  const lines = connectionLines(seqSS, all);
  if (stackedAllowed && isStackable(rest, lines, seqSS.length)) {
    const own = rest.filter((c) => lineOf(lines, c.from) === lineOf(lines, c.to));
    const lanes = assignLanes(own, (c) => lineOf(lines, c.from) === 0 ? 1 : -1);
    const geo = connectionsGeometry(fontSize, false);
    const need = sideHeight(lanes.above, geo) + lines.length * capH + (lines.length - 1) * stackedLineGap(fontSize) +
      Math.max(sideHeight(lanes.below, geo) + 2, met.descBelow) + 2;
    if (h >= need) {
      const between = rest.filter((c) => lineOf(lines, c.from) !== lineOf(lines, c.to));
      return {conns: own.sort(byLaneOuterFirst).concat(between), inline, ...lanes, geo, chains: lines};
    }
  }
  const conns = combineHubs(rest);
  const lanes = assignLanes(conns);
  for (const compressed of [false, true]) {
    const geo = connectionsGeometry(fontSize, compressed);
    if (h >= lineHeightWithLanes(lanes.above, lanes.below, geo, met) + 2)
      return {conns: conns.slice().sort(byLaneOuterFirst), inline, ...lanes, geo, chains: []};
  }
  return null;
}

/** Draws brackets (and combs) in the given order, outer lanes first ({@link layoutConnections} sorts them so).
 * {@link xOf} gives a position's x for an R-group: R1 at the monomer's left edge, R2 at its right edge, a side chain at
 * its centre; positions out of view lie past the cell's clip, so their brackets run off its edge. Brackets of the
 * hovered position are drawn thicker and the others fade. */
export function drawConnections(g: CanvasRenderingContext2D, conns: IDrawnConnection[], xOf: (pos: number, r?: number) => number,
  aboveY0: number, belowY0: number, geo: IConnectionsGeometry, colorOf: (c: IDrawnConnection) => string,
  hoverPos: number | null = null, fade?: boolean,
): void {
  const touches = (c: IDrawnConnection): boolean => touchesPos(c, hoverPos);
  // the others fade when the hovered monomer has a bond here, or the caller says it has one elsewhere in the cell
  const anyHover = hoverPos !== null && (fade ?? conns.some(touches));
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

/** Whether connection {@link c} (or comb) has its end at {@link pos} */
export function touchesPos(c: IDrawnConnection, pos: number | null): boolean {
  return pos !== null && (c.comb ? c.comb.hub === pos || c.comb.members.includes(pos) : c.from === pos || c.to === pos);
}

/** The height a bond's end takes beside a monomer, between the top of the capitals and their foot: R1 at the foot, R2
 * at the top, a side chain in the middle (as a bracket's legs stand at a monomer's left edge, right edge and middle) */
export function inlineJoinY(r: number, top: number, foot: number): number {
  return r === 1 ? foot : r === 2 ? top : (top + foot) / 2;
}

/** Draws the bonds of chains side by side ({@link IDrawnConnection.inline}) where the chain break would be: a line
 * from the left monomer's right edge to the right monomer's left edge, each end at its R-group's height
 * ({@link inlineJoinY}) and marked with a dot, so that a side chain to a side chain (both ends in the middle) does not
 * read as the separator of a chain. {@link ends} gives where the left label ends and the right one begins, and the
 * line's top and foot of the capitals; null skips a bond. The hovered monomer's bond is drawn thicker and, with
 * {@link fade}, the others fade, as the brackets do. */
export function drawInlineJoins(g: CanvasRenderingContext2D, joins: IDrawnConnection[],
  ends: (c: IDrawnConnection) => { x1: number, x2: number, top: number, foot: number } | null,
  lineWidth: number, colorOf: (c: IDrawnConnection) => string, hoverPos: number | null, fade: boolean,
): void {
  if (joins.length === 0) return;
  g.save();
  try {
    g.lineCap = 'round';
    g.setLineDash(SOLID);
    for (const c of joins) {
      const e = ends(c);
      if (!e) continue;
      const hovered = touchesPos(c, hoverPos);
      const color = colorOf(c);
      const w = hovered ? lineWidth * 2 : lineWidth;
      const y1 = inlineJoinY(c.rFrom, e.top, e.foot);
      const y2 = inlineJoinY(c.rTo, e.top, e.foot);
      g.strokeStyle = color;
      g.fillStyle = color;
      g.globalAlpha = fade && !hovered ? 0.3 : 1;
      g.lineWidth = w;
      g.beginPath();
      g.moveTo(e.x1, y1);
      g.lineTo(e.x2, y2);
      g.stroke();
      const r = Math.max(1.1, w * 1.1);
      g.beginPath();
      g.arc(e.x1, y1, r, 0, Math.PI * 2);
      g.moveTo(e.x2 + r, y2);
      g.arc(e.x2, y2, r, 0, Math.PI * 2);
      g.fill();
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
