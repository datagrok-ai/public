/**
 * Pairing attachment labels within one compound.
 *
 * `[*:n]` marks one end of one bond, so a compound is whole exactly when every number it carries is
 * carried by two of its pieces, or by one piece and capped by a position it leaves blank, which is
 * hydrogen. Both things the matrix needs follow from that count: whether a combination could exist,
 * and the order its pieces join in. Read per compound, so a column holding a terminal fragment on one
 * row and a connector on the next describes both instead of neither.
 */

/** The attachments a position's fragments fill, unioned over the column's non-blank values. A blank
 *  carries no number of its own, so the substituted rows are the only thing that says which point it
 *  erases — the one fact a single compound cannot supply about itself. */
export type PositionFills = {[position: string]: number[]};

/** Passes over the pieces, each mapping a position to the attachment it joins at. */
export type LinkStages = Array<{[position: string]: number}>;

export interface FragmentLinks {
  fills: PositionFills;
  /** Values RDKit could read. Text it could not is a label: it pairs with nothing, and joining onto
   *  one yields a query molecule the grid would draw as a real proposed compound. */
  structures: Set<string>;
}

export function attachmentNumbers(smiles: string): Set<number> {
  const numbers = new Set<number>();
  for (const match of smiles.matchAll(/\[\*:(\d+)\]/g))
    numbers.add(Number.parseInt(match[1], 10));
  return numbers;
}

/** The linker rewrites EVERY occurrence of `[*:n]` to the same ring-closure digit, so a piece carrying
 *  it twice closes on itself: core and fragment come back as separate rings that still parse. */
function repeatsAttachment(smiles: string): boolean {
  return [...smiles.matchAll(/\[\*:(\d+)\]/g)].length !== attachmentNumbers(smiles).size;
}

interface Piece {
  /** Empty for the core. */
  position: string;
  numbers: number[];
}

function piecesOf(core: string, values: {[position: string]: string}, positions: string[]):
  {pieces: Piece[], blanks: string[], carriers: Map<number, number[]>} {
  const pieces: Piece[] = [{position: '', numbers: [...attachmentNumbers(core)]}];
  const blanks: string[] = [];
  for (const position of positions) {
    const value = values[position] ?? '';
    if (value === '')
      blanks.push(position);
    else
      pieces.push({position, numbers: [...attachmentNumbers(value)]});
  }
  const carriers = new Map<number, number[]>();
  pieces.forEach((piece, i) => {
    for (const n of piece.numbers) {
      const ends = carriers.get(n);
      if (ends === undefined)
        carriers.set(n, [i]);
      else
        ends.push(i);
    }
  });
  return {pieces, blanks, carriers};
}

/** Whether the core reaches every piece over the bonds the pairing names. */
function connected(pieces: Piece[], carriers: Map<number, number[]>): boolean {
  const seen = new Uint8Array(pieces.length);
  seen[0] = 1;
  let found = 1;
  const queue = [0];
  while (queue.length > 0) {
    const i = queue.pop()!;
    for (const n of pieces[i].numbers) {
      const ends = carriers.get(n)!;
      if (ends.length !== 2)
        continue;
      const other = ends[0] === i ? ends[1] : ends[0];
      if (!seen[other]) {
        seen[other] = 1;
        found++;
        queue.push(other);
      }
    }
  }
  return found === pieces.length;
}

function readable(core: string, values: {[position: string]: string}, positions: string[],
  links: FragmentLinks): boolean {
  return links.structures.has(core) && !repeatsAttachment(core) && positions.every((p) => {
    const value = values[p] ?? '';
    return value === '' || (links.structures.has(value) && !repeatsAttachment(value));
  });
}

/**
 * Whether the decomposition can express this combination, judged on the attachment numbers alone.
 *
 * A limit of the R-group scheme, not of chemistry: a thiophene whose fragment carries no third point
 * may still be substitutable, it just has no R3 site here. Abstains where the numbers prove nothing —
 * a point no column fills anywhere, or text RDKit could not read — since a greyed cell is a claim.
 *
 * Must be answered before anything is built: `linkRGroupFragments` skips a fragment whose point the
 * piece it joins does not carry and hands back what it built, so an impossible cell comes back a
 * clean, drawable molecule of a different compound, usually one already measured in the same matrix.
 */
export function cellPossible(core: string, values: {[position: string]: string}, positions: string[],
  links: FragmentLinks): boolean {
  if (!readable(core, values, positions, links))
    return true;
  const {pieces, blanks, carriers} = piecesOf(core, values, positions);
  const covered = new Set(Object.values(links.fills).flat());
  for (const [n, ends] of carriers) {
    if (ends.length > 2)
      return false;
    if (ends.length === 1 && covered.has(n) && !blanks.some((p) => links.fills[p].includes(n)))
      return false;
  }
  return connected(pieces, carriers);
}

/**
 * The passes that assemble one compound from its own pieces, or null when they cannot be assembled.
 *
 * `reserved` is the site the matrix enumerates: a row key leaves it open, and a row fragment claiming
 * it would weld itself onto the very point the columns vary across. `closed` is what a whole cell
 * needs — a point left open there carries a stray dummy atom into a drawn structure. Null also when a
 * ring closes between fragments, since the linker joins by one bond at a time and a cycle has no walk.
 */
export function planLink(core: string, values: {[position: string]: string}, positions: string[],
  links: FragmentLinks, reserved: number[], closed: boolean): LinkStages | null {
  if (!readable(core, values, positions, links))
    return null;
  const {pieces, blanks, carriers} = piecesOf(core, values, positions);
  const held = new Set(reserved);
  for (const [n, ends] of carriers) {
    if (ends.length > 2 || (ends.length === 2 && held.has(n)))
      return null;
  }
  const stages: LinkStages = [];
  const attached = new Uint8Array(pieces.length);
  const spent = new Set<number>();
  attached[0] = 1;
  let frontier = [0];
  while (frontier.length > 0) {
    const stage: {[position: string]: number} = {};
    const next: number[] = [];
    for (const i of frontier) {
      for (const n of pieces[i].numbers) {
        if (spent.has(n))
          continue;
        spent.add(n);
        const ends = carriers.get(n)!;
        if (ends.length === 1) {
          // Capped here rather than in a pass of its own: the point may only have appeared when the
          // piece exposing it was attached.
          const blank = blanks.find((p) => links.fills[p].includes(n));
          if (blank !== undefined)
            stage[blank] = n;
          else if (closed && !held.has(n))
            return null;
          continue;
        }
        const other = ends[0] === i ? ends[1] : ends[0];
        if (attached[other])
          return null;
        attached[other] = 1;
        stage[pieces[other].position] = n;
        next.push(other);
      }
    }
    if (Object.keys(stage).length > 0)
      stages.push(stage);
    frontier = next;
  }
  return attached.some((a) => a === 0) ? null : stages;
}
