import * as DG from 'datagrok-api/dg';
import {KetSerializer, MolSerializer, Struct} from 'ketcher-core';

/* Ketcher's query features in a molblock RDKit reads (crux-sketch spike query-roundtrip, K1; the product owner, 2026-10-07:
 * "ketcher also has those query features ... lets try to fix"). Ketcher's own V2000 writer (ketcher-core) leaves out the
 * SMARTS-only atom properties (aromatic or aliphatic, implicit H count, ring membership, ring size, connectivity, a custom
 * atom query) and writes a custom bond query as "any"; Indigo's V3000 leaves the atom properties out too. Datagrok's
 * substructure filter reads the sketcher's molblock, so those queries filtered as if they were not there, and came back
 * without them. Here:
 * - `queryMolfile`: Ketcher's V2000 with each custom bond query MDL can hold written as its type and topology, and, for a
 *   drawing with a SMARTS-only atom property, the same drawing as V3000 with an RDKit SMARTSQ data group on every query
 *   atom: its SMARTS is RDKit's reading of the atom's MDL query (as RDKit reads Ketcher's V2000) and its SMARTS-only
 *   properties, in RDKit's SMARTS.
 * - `ketcherQuery`: the reverse, for a V3000 with SMARTSQ groups set on Ketcher (its own, or another sketcher's): each
 *   group's SMARTS on its atom as Ketcher's custom query, the groups removed.
 * Not held: an atom's chirality query (SMARTS `@`/`@@` depend on the neighbours' order), a custom bond query MDL's types
 * and topology cannot say (written "any", as before), and a drawing with S-groups, R-groups or attachment points beside a
 * SMARTS-only property (written V2000, as before). */

/** The SMARTS-only properties of a Ketcher atom, in RDKit's SMARTS: aromatic `a`, aliphatic `A`, implicit H count `h<n>`,
 * ring membership `R<n>`, ring size `r<n>`, connectivity `X<n>`. */
export interface QueryAtomFields {
  queryProperties?: {aromaticity?: string | null, ringMembership?: number | null, ringSize?: number | null,
    connectivity?: number | null, customQuery?: string | null} | null;
  implicitHCount?: number | null;
}

export function smartsOnly(atom: QueryAtomFields): string[] {
  const q = atom.queryProperties ?? {};
  const out: string[] = [];
  if (q.aromaticity === 'aromatic')
    out.push('a');
  else if (q.aromaticity === 'aliphatic')
    out.push('A');
  if (atom.implicitHCount != null)
    out.push(`h${atom.implicitHCount}`);
  if (q.ringMembership != null)
    out.push(`R${q.ringMembership}`);
  if (q.ringSize != null)
    out.push(`r${q.ringSize}`);
  if (q.connectivity != null)
    out.push(`X${q.connectivity}`);
  return out;
}

/** A custom query of an atom, or null. */
export function customAtomQuery(atom: QueryAtomFields): string | null {
  const c = atom.queryProperties?.customQuery?.trim();
  return c ? c : null;
}

/** Whether the atom has a property only SMARTS holds. */
export function needsSmarts(atom: QueryAtomFields): boolean {
  return customAtomQuery(atom) !== null || smartsOnly(atom).length > 0;
}

const BOND_ORDERS: {[k: string]: number} = {'-': 1, '=': 2, '#': 3, ':': 4};

/**
 * The MDL bond type (1 to 8) and topology (0 either, 1 ring, 2 chain) a custom bond query means, or null when MDL cannot
 * say it: an OR of single, double, triple and aromatic that MDL has a type for (`-`, `=`, `#`, `:`, `-,=`, `-,:`, `=,:`,
 * `~`), optionally ANDed with `@` or `!@` (by `;`, or by `&` after a single order), or `@` or `!@` alone (any bond).
 */
export function mdlBondQuery(custom: string): {type: number, topology: number} | null {
  const text = custom.replace(/\s+/g, '');
  // a lone topology: any bond, in a ring or in a chain
  const lone = /^(!?@)$/.exec(text);
  if (lone)
    return {type: 8, topology: lone[1] === '@' ? 1 : 2};
  const m = /^([-=#:~](?:,[-=#:~])*)(?:([;&])(!?@))?$/.exec(text);
  // `-,=&@` is `-,(=&@)` to a SMARTS parser: an OR of orders and a topology joined only by the low-precedence `;`
  if (!m || (m[2] === '&' && m[1].includes(',')))
    return null;
  const topology = m[3] === '@' ? 1 : m[3] === '!@' ? 2 : 0;
  const parts = Array.from(new Set(m[1].split(',')));
  if (parts.includes('~'))
    return parts.length === 1 ? {type: 8, topology} : null;
  const orders = parts.map((p) => BOND_ORDERS[p]).sort().join(',');
  const type = ({'1': 1, '2': 2, '3': 3, '4': 4, '1,2': 5, '1,4': 6, '2,4': 7} as {[k: string]: number})[orders];
  return type === undefined ? null : {type, topology};
}

const pad = (n: number | string, w: number) => String(n).padStart(w, ' ');

/** The V2000 lines a reader of this module handles: anything else (S-groups, R-groups, aliases …) leaves the molblock V2000. */
const HANDLED_PROPS = /^M {2}(CHG|RAD|ISO|ALS|RBC|SUB|UNS|END)/;

interface V2000Parsed {
  header: string[];
  atoms: {x: string, y: string, z: string, symbol: string, line: string}[];
  bonds: {a: number, b: number, type: number, stereo: number, topology: number}[];
  charges: Map<number, number>;
  radicals: Map<number, number>;
  isotopes: Map<number, number>;
  lists: Map<number, {not: boolean, symbols: string[]}>;
  others: string[];
}

export function parseV2000(text: string): V2000Parsed | null {
  const lines = text.split(/\r?\n/);
  const counts = lines[3];
  if (!counts || !counts.includes('V2000'))
    return null;
  const na = parseInt(counts.slice(0, 3));
  const nb = parseInt(counts.slice(3, 6));
  if (!Number.isFinite(na) || !Number.isFinite(nb))
    return null;
  const out: V2000Parsed = {header: lines.slice(0, 3), atoms: [], bonds: [], charges: new Map(), radicals: new Map(),
    isotopes: new Map(), lists: new Map(), others: []};
  for (let i = 0; i < na; i++) {
    const l = lines[4 + i] ?? '';
    out.atoms.push({x: l.slice(0, 10).trim(), y: l.slice(10, 20).trim(), z: l.slice(20, 30).trim(), symbol: l.slice(31, 34).trim(),
      line: l});
    const ccc = parseInt(l.slice(36, 39));
    if (ccc > 0 && ccc !== 4)
      out.charges.set(i + 1, 4 - ccc);
    if (ccc === 4)
      out.radicals.set(i + 1, 2);
  }
  for (let i = 0; i < nb; i++) {
    const l = lines[4 + na + i] ?? '';
    out.bonds.push({a: parseInt(l.slice(0, 3)), b: parseInt(l.slice(3, 6)), type: parseInt(l.slice(6, 9)),
      stereo: parseInt(l.slice(9, 12)) || 0, topology: parseInt(l.slice(15, 18)) || 0});
  }
  for (const l of lines.slice(4 + na + nb)) {
    if (!l.trim())
      continue;
    const key = l.slice(0, 6);
    const pairs = (into: Map<number, number>) => {
      const n = parseInt(l.slice(6, 9));
      for (let k = 0; k < n; k++)
        into.set(parseInt(l.slice(9 + 8 * k, 13 + 8 * k)), parseInt(l.slice(13 + 8 * k, 17 + 8 * k)));
    };
    if (key === 'M  CHG') {
      for (const [a] of out.charges)
        out.charges.delete(a);
      pairs(out.charges);
    } else if (key === 'M  RAD')
      pairs(out.radicals);
    else if (key === 'M  ISO')
      pairs(out.isotopes);
    else if (key === 'M  ALS') {
      const atom = parseInt(l.slice(7, 10));
      const n = parseInt(l.slice(10, 13));
      const not = l.slice(14, 15) === 'T';
      const symbols: string[] = [];
      for (let k = 0; k < n; k++)
        symbols.push(l.slice(16 + 4 * k, 20 + 4 * k).trim());
      out.lists.set(atom, {not, symbols});
    } else if (!HANDLED_PROPS.test(l))
      out.others.push(l);
  }
  return out;
}

/** The V2000 with each atom numbered by its atom-atom mapping field (1, 2 …): RDKit then writes each atom's expression in
 * SMARTS with that number. */
export function numbered(text: string): string {
  const lines = text.split(/\r?\n/);
  const na = parseInt(lines[3].slice(0, 3));
  for (let i = 0; i < na; i++) {
    const l = lines[4 + i].padEnd(69, ' ');
    lines[4 + i] = l.slice(0, 60) + pad(i + 1, 3) + l.slice(63);
  }
  return lines.join('\n');
}

/** Each atom's expression in a SMARTS whose atoms carry their numbers (`[#6&x0:2]`), by number, without the number. */
export function atomExpressions(smarts: string): Map<number, string> {
  const out = new Map<number, string>();
  for (let i = 0; i < smarts.length; i++) {
    if (smarts[i] !== '[')
      continue;
    let depth = 0;
    let j = i;
    for (; j < smarts.length; j++) {
      if (smarts[j] === '[')
        depth++;
      else if (smarts[j] === ']' && --depth === 0)
        break;
    }
    const inner = smarts.slice(i + 1, j);
    const m = /^(.*):(\d+)$/.exec(inner);
    if (m)
      out.set(parseInt(m[2]), m[1]);
    i = j;
  }
  return out;
}

/** A SMARTSQ data group's lines, its SMARTS split over continuation lines where it is long. */
function smartsqLines(index: number, atom: number, smarts: string): string[] {
  const out = [`M  V30 ${index} DAT 0 ATOMS=(1 ${atom}) QUERYTYPE=SMARTSQ QUERYOP== -`];
  let rest = `FIELDDATA="${smarts}"`;
  while (rest.length > 70) {
    out.push(`M  V30 ${rest.slice(0, 70)}-`);
    rest = rest.slice(70);
  }
  out.push(`M  V30 ${rest}`);
  return out;
}

/**
 * Ketcher's V2000 (`v2000`, ketcher-core's) of a drawing whose atoms and bonds are `atoms` and `bonds` in its order, made a
 * molblock RDKit reads as the query Ketcher shows: each custom bond query MDL can say written as its type and topology;
 * and, where an atom has a SMARTS-only property, the V3000 of it with a SMARTSQ group on every query atom.
 * `rdkitSmarts(molblock)` is RDKit's SMARTS of a query molblock (get_qmol). Returns `v2000` itself when there is nothing
 * to add.
 */
export function queryMolfile(v2000: string, atoms: QueryAtomFields[], bonds: {customQuery?: string | null}[],
  rdkitSmarts: (molblock: string) => string): string {
  const lines = v2000.split(/\r?\n/);
  const counts = lines[3] ?? '';
  const na = parseInt(counts.slice(0, 3));
  // custom bond queries as MDL's type and topology
  bonds.forEach((b, k) => {
    const custom = b.customQuery?.trim();
    const mdl = custom ? mdlBondQuery(custom) : null;
    const li = 4 + na + k;
    if (mdl === null || lines[li] === undefined)
      return;
    const l = lines[li].padEnd(21, ' ');
    lines[li] = l.slice(0, 6) + pad(mdl.type, 3) + l.slice(9, 15) + pad(mdl.topology, 3) + l.slice(18);
  });
  const patched = lines.join('\n');
  if (!atoms.some(needsSmarts))
    return patched;
  const parsed = parseV2000(patched);
  if (parsed === null || parsed.others.length > 0 || parsed.atoms.length !== atoms.length)
    return patched;
  const expressions = atomExpressions(rdkitSmarts(numbered(patched)));
  // every atom RDKit reads as more than an element (with its charge or isotope) gets its SMARTS: a list, a generic, an MDL
  // query field (as RDKit resolves it: a ring bond count "as drawn" is the count drawn, which RDKit's V3000 reader would
  // read as `x-2`), and the SMARTS-only properties
  const isQueryAtom = (k: number) => needsSmarts(atoms[k]) ||
    !/^(\d+)?#\d+(&[+-]\d*)?$/.test(expressions.get(k + 1) ?? '#0');
  const sgroups: string[] = [];
  let n = 0;
  const atomLines = parsed.atoms.map((a, k) => {
    const list = parsed.lists.get(k + 1);
    const symbol = !list ? a.symbol : list.not ? `"NOT [${list.symbols.join(',')}]"` : `[${list.symbols.join(',')}]`;
    const props = [
      parsed.charges.has(k + 1) && parsed.charges.get(k + 1) !== 0 ? `CHG=${parsed.charges.get(k + 1)}` : '',
      parsed.radicals.has(k + 1) ? `RAD=${parsed.radicals.get(k + 1)}` : '',
      parsed.isotopes.has(k + 1) ? `MASS=${parsed.isotopes.get(k + 1)}` : '',
    ].filter(Boolean).join(' ');
    if (isQueryAtom(k)) {
      const custom = customAtomQuery(atoms[k]);
      const base = expressions.get(k + 1);
      const extras = smartsOnly(atoms[k]);
      const smarts = custom !== null ? `[${custom.replace(/^\[(.*)\]$/, '$1')}]` :
        base !== undefined ? `[${base}${extras.length ? `;${extras.join('&')}` : ''}]` : null;
      if (smarts !== null)
        sgroups.push(...smartsqLines(++n, k + 1, smarts));
    }
    return `M  V30 ${k + 1} ${symbol} ${a.x} ${a.y} ${a.z} 0${props ? ' ' + props : ''}`;
  });
  const cfg = (b: {type: number, stereo: number}) => b.type === 1 ? ({1: 1, 4: 2, 6: 3} as {[s: number]: number})[b.stereo] ?? 0 :
    b.type === 2 && b.stereo === 3 ? 2 : 0;
  const bondLines = parsed.bonds.map((b, k) => {
    const c = cfg(b);
    return `M  V30 ${k + 1} ${b.type} ${b.a} ${b.b}${c ? ` CFG=${c}` : ''}${b.topology ? ` TOPO=${b.topology}` : ''}`;
  });
  return [...parsed.header, '  0  0  0     0  0            999 V3000', 'M  V30 BEGIN CTAB',
    `M  V30 COUNTS ${parsed.atoms.length} ${parsed.bonds.length} ${n} 0 0`, 'M  V30 BEGIN ATOM', ...atomLines,
    'M  V30 END ATOM', ...(bondLines.length ? ['M  V30 BEGIN BOND', ...bondLines, 'M  V30 END BOND'] : []),
    ...(n ? ['M  V30 BEGIN SGROUP', ...sgroups, 'M  V30 END SGROUP'] : []), 'M  V30 END CTAB', 'M  END', ''].join('\n');
}

/** The SMARTSQ groups of a V3000 molblock, by atom (1-based): `{atom: smarts}`; empty when it has none. */
export function smartsqGroups(v3000: string): Map<number, string> {
  const out = new Map<number, string>();
  if (!v3000.includes('V3000') || !v3000.includes('SMARTSQ'))
    return out;
  // join continuation lines, then read each DAT group's atoms and data
  const joined: string[] = [];
  let carry = '';
  for (const raw of v3000.split(/\r?\n/)) {
    if (!raw.startsWith('M  V30 ')) {
      joined.push(raw);
      continue;
    }
    const body = raw.slice(7);
    if (body.endsWith('-')) {
      carry += body.slice(0, -1);
      continue;
    }
    joined.push('M  V30 ' + carry + body);
    carry = '';
  }
  for (const l of joined) {
    const m = /^M {2}V30 \d+ DAT \d+ .*ATOMS=\(1 (\d+)\).*QUERYTYPE=SMARTSQ.*FIELDDATA="([^"]*)"/.exec(l);
    if (m)
      out.set(parseInt(m[1]), m[2]);
  }
  return out;
}

/** An atom's SMARTSQ as Ketcher's custom query: without its brackets and its atom-map number. */
export function customQueryOf(smarts: string): string {
  return smarts.trim().replace(/^\[(.*)\]$/, '$1').replace(/:\d+$/, '');
}

/** The V3000 without its SMARTSQ groups (the S-group block gone, its count 0), or null when it holds other S-groups too. */
export function withoutSmartsq(v3000: string): string | null {
  const lines = v3000.split(/\r?\n/);
  const begin = lines.findIndex((l) => l.startsWith('M  V30 BEGIN SGROUP'));
  const end = lines.findIndex((l) => l.startsWith('M  V30 END SGROUP'));
  if (begin < 0 || end < begin)
    return v3000;
  const block = lines.slice(begin + 1, end);
  // each group's first line starts with its index; a continuation does not
  const starts = block.filter((l, k) => /^M {2}V30 \d+ \w+ /.test(l) && !(k > 0 && block[k - 1].endsWith('-')));
  if (starts.some((l) => !/^M {2}V30 \d+ DAT /.test(l)) || smartsqGroups(v3000).size !== starts.length)
    return null;
  const out = [...lines.slice(0, begin), ...lines.slice(end + 1)];
  return out.map((l) => l.startsWith('M  V30 COUNTS ') ? l.replace(/^(M {2}V30 COUNTS \d+ \d+ )\d+/, '$10') : l).join('\n');
}


// ---------------------------------------------------------------- Ketcher's structures

/** Chem's RDKit module, as `DG.chem.convert` reaches it: Chem's function, called synchronously. */
let rdkit: any = null;

/** RDKit's SMARTS of a query molblock, read as a substructure search reads one (`get_qmol`). */
export function rdkitSmarts(molblock: string): string {
  if (rdkit === null) {
    const call = DG.Func.find({package: 'Chem', name: 'getRdKitModule'})[0].prepare();
    call.callSync();
    rdkit = call.getOutputParamValue();
  }
  const mol = rdkit.get_qmol(molblock);
  try {
    return mol.get_smarts();
  } finally {
    mol.delete();
  }
}

/** The structure each query molfile written here came from, as KET, by that molfile: set back on Ketcher (the
 * substructure filter's sketcher reopened, a switch back to Ketcher), it shows the drawing as it was drawn, its SMARTS-only
 * properties as fields. The last 32. */
const writtenFrom = new Map<string, string>();
const REMEMBERED = 32;

/** Ketcher's V2000 of `drawn` (`v2000`), made the query molblock RDKit reads as Ketcher shows it (`queryMolfile`). */
export function withQueries(drawn: Struct, v2000: string): string {
  const atoms = [...drawn.atoms.values()] as QueryAtomFields[];
  const bonds = [...drawn.bonds.values()] as {customQuery?: string | null}[];
  if (!atoms.some(needsSmarts) && !bonds.some((b) => b.customQuery?.trim()))
    return v2000;
  const text = queryMolfile(v2000, atoms, bonds, rdkitSmarts);
  if (text !== v2000) {
    writtenFrom.delete(text);
    writtenFrom.set(text, new KetSerializer().serialize(drawn));
    if (writtenFrom.size > REMEMBERED)
      writtenFrom.delete(writtenFrom.keys().next().value!);
  }
  return text;
}

/**
 * What Ketcher is given to show for a molecule string: for a query molfile written here, the structure it came from; for
 * another V3000 with SMARTSQ groups (Crux's), the same with each group's SMARTS on its atom as Ketcher's custom query;
 * anything else as it is.
 */
export function asKetcherQuery(text: string): string {
  const remembered = writtenFrom.get(text);
  if (remembered !== undefined)
    return remembered;
  const groups = smartsqGroups(text);
  if (groups.size === 0)
    return text;
  const stripped = withoutSmartsq(text);
  if (stripped === null)
    return text;
  const struct = new MolSerializer().deserialize(stripped);
  const ids = [...struct.atoms.keys()];
  for (const [atom, smarts] of groups) {
    const a: any = struct.atoms.get(ids[atom - 1]);
    if (a)
      a.queryProperties = {...(a.queryProperties ?? {}), customQuery: customQueryOf(smarts)};
  }
  return new KetSerializer().serialize(struct);
}
