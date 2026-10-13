/**
 * Pairing attachment labels within one compound. `[*:n]` is one end of one bond, so a compound is whole
 * when every number is carried by two of its pieces, or by one piece and capped by a blank position.
 */

/** Attachment numbers per position, over the column's non-blank values. */
export type PositionFills = {[position: string]: number[]};

/** Join passes, each mapping a position to the attachment numbers it joins at. */
export type LinkStages = Array<{[position: string]: number[]}>;

export interface FragmentLinks {
  /** Every number a position's fragments carry: what a blank there caps. */
  fills: PositionFills;
  /** The numbers most of a position's fragments carry: what a row key leaves open for the axis. */
  sites: PositionFills;
  /** Values RDKit could read. Anything else is a label and never joins. */
  structures: Set<string>;
  covered: Set<number>;
  parsed: Map<string, {numbers: number[], repeats: boolean}>;
}

export function fragmentLinks(fills: PositionFills, sites: PositionFills, structures: Set<string>): FragmentLinks {
  return {fills, sites, structures, covered: new Set(Object.values(fills).flat()), parsed: new Map()};
}

export function attachmentNumbers(smiles: string): Set<number> {
  const numbers = new Set<number>();
  for (const match of smiles.matchAll(/\[\*:(\d+)\]/g))
    numbers.add(Number.parseInt(match[1], 10));
  return numbers;
}

function parsed(links: FragmentLinks, smiles: string): {numbers: number[], repeats: boolean} {
  let entry = links.parsed.get(smiles);
  if (entry === undefined) {
    const numbers = attachmentNumbers(smiles);
    // The linker turns every `[*:n]` into one ring digit, so a piece carrying n twice closes on itself.
    const repeats = [...smiles.matchAll(/\[\*:\d+\]/g)].length !== numbers.size;
    links.parsed.set(smiles, entry = {numbers: [...numbers], repeats});
  }
  return entry;
}

interface Piece {
  /** Empty for the core. */
  position: string;
  numbers: number[];
}

function piecesOf(core: string, values: {[position: string]: string}, positions: string[], links: FragmentLinks):
  {pieces: Piece[], blanks: string[], carriers: Map<number, number[]>} {
  const pieces: Piece[] = [{position: '', numbers: parsed(links, core).numbers}];
  const blanks: string[] = [];
  const bridges = new Set<string>();
  for (const position of positions) {
    const value = values[position] ?? '';
    if (value === '') {
      blanks.push(position);
      continue;
    }
    const numbers = parsed(links, value).numbers;
    // A fragment spanning two attachment points is written in each column it spans; it is one piece.
    if (numbers.length > 1 && bridges.has(value))
      continue;
    bridges.add(value);
    pieces.push({position, numbers});
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
  return links.structures.has(core) && !parsed(links, core).repeats && positions.every((p) => {
    const value = values[p] ?? '';
    return value === '' || (links.structures.has(value) && !parsed(links, value).repeats);
  });
}

/**
 * Whether the R-groups can form this combination, judged on the attachment numbers alone. Answers
 * true where the numbers prove nothing: a label, a structure with no numbered attachment point, or a
 * point no column fills.
 */
export function cellPossible(core: string, values: {[position: string]: string}, positions: string[],
  links: FragmentLinks): boolean {
  if (!readable(core, values, positions, links))
    return true;
  const {pieces, blanks, carriers} = piecesOf(core, values, positions, links);
  if (pieces.some((piece) => piece.numbers.length === 0))
    return true;
  for (const [n, ends] of carriers) {
    if (ends.length > 2)
      return false;
    if (ends.length === 1 && links.covered.has(n) && !blanks.some((p) => links.fills[p].includes(n)))
      return false;
  }
  return connected(pieces, carriers);
}

/**
 * The passes that join one compound from the core outwards, or null when they cannot. A row key
 * (`closed` false) leaves `reserved`, the axis's site, open; a whole cell must close every point.
 * A ring closing between two fragments is refused: the linker only joins onto what is already built.
 */
export function planLink(core: string, values: {[position: string]: string}, positions: string[],
  links: FragmentLinks, reserved: number[], closed: boolean): LinkStages | null {
  if (!readable(core, values, positions, links))
    return null;
  const {pieces, blanks, carriers} = piecesOf(core, values, positions, links);
  if (pieces.some((piece, i) => i > 0 && piece.numbers.length === 0))
    return null;
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
    const built = attached.slice();
    const stage: {[position: string]: number[]} = {};
    const add = (position: string, n: number): void => {
      (stage[position] ??= []).push(n);
    };
    const next: number[] = [];
    for (const i of frontier) {
      for (const n of pieces[i].numbers) {
        if (spent.has(n))
          continue;
        spent.add(n);
        const ends = carriers.get(n)!;
        if (ends.length === 1) {
          const blank = blanks.find((p) => links.fills[p].includes(n));
          if (blank !== undefined)
            add(blank, n);
          else if (closed && !held.has(n))
            return null;
          continue;
        }
        const other = ends[0] === i ? ends[1] : ends[0];
        if (attached[other])
          return null;
        attached[other] = 1;
        add(pieces[other].position, n);
        // A bridge forms every bond it shares with the built part in the same join.
        for (const m of pieces[other].numbers) {
          const partners = carriers.get(m)!;
          const partner = partners[0] === other ? partners[1] : partners[0];
          if (!spent.has(m) && partners.length === 2 && built[partner]) {
            spent.add(m);
            add(pieces[other].position, m);
          }
        }
        next.push(other);
      }
    }
    if (Object.keys(stage).length > 0)
      stages.push(stage);
    frontier = next;
  }
  return attached.some((a) => a === 0) ? null : stages;
}
