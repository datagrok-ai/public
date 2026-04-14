/**
 * Pure utility functions for the interactive atom picker.
 *
 * Extracted from rdkit-cell-renderer.ts to keep the renderer class
 * focused on rendering while these handle geometry, SVG parsing,
 * and bond selection computations.
 */

/** Hit-tolerance in CSS pixels for atom proximity detection. */
const CLICK_TOLERANCE_PX = 12;

/** Returns the index of the atom whose center is nearest to the given
 *  point, but only if it falls within CLICK_TOLERANCE_PX (CSS pixels).
 *  Returns null if the point is too far from all atoms. */
export function findNearestAtom(
  positions: Map<number, {x: number; y: number}>,
  clickX: number,
  clickY: number,
): number | null {
  const tolSq = CLICK_TOLERANCE_PX * CLICK_TOLERANCE_PX;
  let nearestIdx: number | null = null;
  let nearestDistSq = tolSq;
  for (const [idx, p] of positions.entries()) {
    const dx = p.x - clickX;
    const dy = p.y - clickY;
    const dsq = dx * dx + dy * dy;
    if (dsq < nearestDistSq) {
      nearestDistSq = dsq;
      nearestIdx = idx;
    }
  }
  return nearestIdx;
}

/** Given a set of selected atom indices and a bondAtoms map, returns the
 *  bond indices where BOTH endpoint atoms are in the selection, and a
 *  matching highlightBondColors map with every selected bond painted the
 *  same color. */
export function computeSelectedBonds(
  atoms: Set<number>,
  bondAtoms: Map<number, [number, number]>,
  color: number[],
): {bondsArr: number[]; highlightBondColors: {[k: number]: number[]}} {
  const bondsArr: number[] = [];
  const highlightBondColors: {[k: number]: number[]} = {};
  for (const [bondIdx, [a1, a2]] of bondAtoms.entries()) {
    if (atoms.has(a1) && atoms.has(a2)) {
      bondsArr.push(bondIdx);
      highlightBondColors[bondIdx] = color;
    }
  }
  return {bondsArr, highlightBondColors};
}

/** Extracts atom positions and bond connectivity from an RDKit-generated
 *  SVG element. Atom positions come from two sources:
 *  - Text-labeled heteroatoms: center of the `<text class="atom-N">` bbox.
 *  - Carbons: average of bond-path endpoints touching that atom.
 *
 *  Bond connectivity is read from `<path class="bond-K atom-A atom-B">`
 *  class attributes. */
export function extractAtomPositionsFromSvg(
  svgEl: SVGSVGElement,
): {positions: Map<number, {x: number; y: number}>; bondAtoms: Map<number, [number, number]>} {
  const positions = new Map<number, {x: number; y: number}>();
  const bondAtoms = new Map<number, [number, number]>();

  // Text-labeled heteroatoms.
  const texts = svgEl.querySelectorAll('text[class*="atom-"]');
  for (let i = 0; i < texts.length; i++) {
    const t = texts[i];
    const cls = t.getAttribute('class') || '';
    const m = /(?:^|\s)atom-(\d+)/.exec(cls);
    if (!m) continue;
    const idx = parseInt(m[1], 10);
    try {
      const bb = (t as unknown as SVGGraphicsElement).getBBox();
      positions.set(idx, {x: bb.x + bb.width / 2, y: bb.y + bb.height / 2});
    } catch {
      /* ignore */
    }
  }

  // Carbons — average bond endpoints from <path> d attributes.
  const bondEnds = new Map<number, Array<{x: number; y: number}>>();
  const bondPaths = svgEl.querySelectorAll('path[class*="bond-"]');
  for (let i = 0; i < bondPaths.length; i++) {
    const p = bondPaths[i];
    const cls = p.getAttribute('class') || '';
    const atomIds: number[] = [];
    const atomMatches = cls.match(/atom-(\d+)/g) || [];
    for (const mm of atomMatches) {
      const n = parseInt(mm.slice(5), 10);
      if (!Number.isNaN(n))
        atomIds.push(n);
    }
    if (atomIds.length !== 2) continue;

    // Record the bond-index → atom-pair mapping.
    const bondMatch = /(?:^|\s)bond-(\d+)/.exec(cls);
    if (bondMatch) {
      const bondIdx = parseInt(bondMatch[1], 10);
      bondAtoms.set(bondIdx, [atomIds[0], atomIds[1]]);
    }

    const d = p.getAttribute('d') || '';
    const coords: Array<{x: number; y: number}> = [];
    const re = /[ML]\s*([\-\d.]+)[\s,]+([\-\d.]+)/g;
    let mm: RegExpExecArray | null;
    while ((mm = re.exec(d)) !== null)
      coords.push({x: parseFloat(mm[1]), y: parseFloat(mm[2])});
    if (coords.length < 2) continue;
    if (!bondEnds.has(atomIds[0]))
      bondEnds.set(atomIds[0], []);
    if (!bondEnds.has(atomIds[1]))
      bondEnds.set(atomIds[1], []);
    bondEnds.get(atomIds[0])!.push(coords[0]);
    bondEnds.get(atomIds[1])!.push(coords[coords.length - 1]);
  }
  for (const [idx, pts] of bondEnds.entries()) {
    if (positions.has(idx)) continue;
    const cx = pts.reduce((s, p) => s + p.x, 0) / pts.length;
    const cy = pts.reduce((s, p) => s + p.y, 0) / pts.length;
    positions.set(idx, {x: cx, y: cy});
  }

  return {positions, bondAtoms};
}
