/* eslint-disable max-len */
import {isMolBlock} from './chem-common';

type AliasTranslation = {
  symbol: string;                                       // RDKit-friendly atom symbol (A, Q, C, *)
  rbc?: number;                                         // M  RBC value: -1=acyclic, 1-3=exact, 4=>=4
  als?: {negate: boolean; atoms: string[]};             // M  ALS atom list
};

// Metals matched by Ketcher's `M` / `MH` generic-atom buttons. The strict
// CTfile spec defines 16f4 for M  ALS (16 atom-symbol slots per line) but
// buildAlsLine emits the full list regardless; RDKit's V2000 parser in
// practice handles longer ALS lines for our use here.
const METAL_ATOMS = [
  'Li', 'Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba',
  'Al', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
  'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
  'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'];
const HALOGEN_ATOMS = ['F', 'Cl', 'Br', 'I'];

// Query-atom aliases that sketchers emit but RDKit's V2000 parser silently
// downgrades to dummy atoms. Sources verified empirically:
//   - Ketcher "Extended Table" short codes — observed in molblocks emitted by
//     EPAM Ketcher (atom-symbol field, cols 32-34).
//   - Marvin JS alias-block text (e.g. `Heterocyclyl`) — observed in molblocks
//     emitted by ChemAxon Marvin (`A  <idx>\n<text>` after the bond block).
//
// "H-inclusive" Ketcher variants (CYH, CBH, CHH, CXH, ACH, ABH) collapse to
// their non-H counterparts because H atoms can't satisfy a ring / heteroatom
// constraint anyway. The truly broad H-inclusive variants (AH, QH, AHH) cannot
// be faithfully translated in V2000 and are listed in UNSUPPORTED_ALIASES.
const QUERY_ALIASES: Record<string, AliasTranslation> = {
  // -- Acyclic (any element) --
  // ACY: non-H, acyclic. ACH (= acyclic incl. H) is intentionally NOT here —
  // V2000 has no atom-level "include H" form (`*` parses as dummy [#0] in
  // RDKit), so it's listed in UNSUPPORTED_ALIASES and hidden in the sketcher
  // rather than approximated to ACY.
  'ACY': {symbol: 'A', rbc: -1},

  // -- Cyclic (any element) --
  // CYC: non-H, in ring. CYH: in ring OR H — encoded identically because H
  // can never be a ring member. The CYH button is hidden in the filter
  // sketcher (editor.css) since it's indistinguishable from CYC in practice,
  // but the entry is kept so externally-supplied molblocks still translate.
  // RBC=2 covers typical interior ring atoms; bridgeheads (≥3 ring bonds) miss.
  'CYC': {symbol: 'A', rbc: 2},
  'CYH': {symbol: 'A', rbc: 2},

  // -- Acyclic carbon --
  // ABC: just C. ABH (= C or H) is technically encoded as an atom list, but
  // RDKit drops the H entry during substructure search (explicit Hs in target
  // don't match), making it indistinguishable from ABC in practice. The ABH
  // button is hidden in editor.css; entry kept for externally-pasted molblocks.
  'ABC': {symbol: 'C', rbc: -1},
  'ABH': {symbol: 'L', rbc: -1, als: {negate: false, atoms: ['C', 'H']}},

  // -- Acyclic hetero --
  // AHC: heteroatom (not C, not H). AHH (= acyclic anything except C, incl. H)
  // intentionally NOT here — RDKit's V2000 parser ignores `M  ALS T` (exclude)
  // lists and treats the L atom as a plain wildcard, producing `[*;!R]` which
  // matches far too much. AHH is in UNSUPPORTED_ALIASES and hidden in the UI.
  'AHC': {symbol: 'Q', rbc: -1},

  // -- Cyclic carbon --
  // CBH ≡ CBC: H can't be ring-member. CBH button hidden in editor.css.
  'CBC': {symbol: 'C', rbc: 2},
  'CBH': {symbol: 'C', rbc: 2},

  // -- Cyclic hetero --
  // CHH ≡ CHC: H can't be ring-member. CHH button hidden in editor.css.
  'CHC': {symbol: 'Q', rbc: 2},
  'CHH': {symbol: 'Q', rbc: 2},

  // -- Cyclic non-carbon --
  // CXX in V2000 reduces to "cyclic heteroatom" (same as CHC); CXH collapses
  // to CXX since H can't be ring-member. CXH button hidden in editor.css.
  'CXX': {symbol: 'Q', rbc: 2},
  'CXH': {symbol: 'Q', rbc: 2},
  // Halogen / metal atom lists. Use 'L' (canonical atom-list placeholder) as
  // the atom-symbol field — with 'A' RDKit's V2000 parser treats the atom as
  // "any non-H" and ignores the M  ALS line.
  //
  // XH/MH technically include H in the list, but RDKit drops the H entry
  // during substructure matching — the +H semantics never surface in practice.
  // Both buttons are hidden in editor.css; entries kept so externally-pasted
  // molblocks still translate.
  'X':  {symbol: 'L', als: {negate: false, atoms: HALOGEN_ATOMS}},
  'XH': {symbol: 'L', als: {negate: false, atoms: [...HALOGEN_ATOMS, 'H']}},
  'M':  {symbol: 'L', als: {negate: false, atoms: METAL_ATOMS}},
  'MH': {symbol: 'L', als: {negate: false, atoms: [...METAL_ATOMS, 'H']}},
  // Marvin JS alias-block text observed in the wild
  'Heterocyclyl': {symbol: 'Q', rbc: 2},
};

// Ketcher labels we can recognise but cannot faithfully express in V2000.
// translateQueryAliasesV2000 logs a warning for each occurrence so the user
// understands why the substructure search may behave unexpectedly. Recommended
// fix when one of these is needed: build the query as SMARTS instead of via
// the sketcher, or convert the molblock to V3000 (some — like ARY, HAR — still
// require SMARTS even in V3000 since aromaticity isn't a standard query bit).
export const UNSUPPORTED_ALIASES: Record<string, string> = {
  'AH':  'any atom including H — V2000 has no clean way; * is parsed as dummy [#0] in RDKit',
  'ACH': 'acyclic atom incl. H — V2000 cannot express "any element including H" at the atom level',
  'AHH': 'acyclic anything except C — RDKit V2000 ignores M ALS T (exclude) lists; verified empirically (test returns wildcard-acyclic matches, not non-C-acyclic)',
  '*': 'wildcard (any atom incl. H per Ketcher) — RDKit V2000 parses * as dummy [#0], not a true wildcard',
  'QH':  'heteroatom or H (= not C) — needs L + M ALS exclude [C], rendered same way as AHH but without ring constraint; currently not encoded',
  'ALK': 'alkyl (sp3 acyclic carbon) — V2000 has no atom-level hybridization query',
  'ALH': 'alkyl incl. H — V2000 has no atom-level hybridization query',
  'AEL': 'alkenyl (sp2 acyclic carbon) — V2000 has no atom-level hybridization query',
  'AEH': 'alkenyl incl. H — V2000 has no atom-level hybridization query',
  'AYL': 'alkynyl (sp acyclic carbon) — V2000 has no atom-level hybridization query',
  'AYH': 'alkynyl incl. H — V2000 has no atom-level hybridization query',
  'AOX': 'Ketcher AOX (semantics undocumented; likely acyclic O-functional group)',
  'AOH': 'Ketcher AOH (semantics undocumented)',
  'ARY': 'aryl (aromatic ring carbon) — V2000 has no atom-level aromaticity query',
  'ARH': 'aryl incl. H — V2000 has no atom-level aromaticity query',
  'CAL': 'cycloalkyl (sp3 ring carbon) — V2000 has no atom-level hybridization query',
  'CAH': 'cycloalkyl incl. H — V2000 has no atom-level hybridization query',
  'CEL': 'cycloalkenyl (sp2 ring carbon) — V2000 has no atom-level hybridization query',
  'CEH': 'cycloalkenyl incl. H — V2000 has no atom-level hybridization query',
  'HAR': 'hetero aryl — V2000 has no atom-level aromaticity query',
  'HAH': 'hetero aryl incl. H — V2000 has no atom-level aromaticity query',
  'G':   'group generic — represents a functional group, not a single atom',
  'GH':  'group generic incl. H — represents a functional group',
  'G*':  'group generic — represents a functional group',
  'GH*': 'group generic — represents a functional group',
  'Pol': 'polymer attachment — SRU notation, not an atom-level query',
};

function pad3(n: number | string): string {
  return n.toString().padStart(3, ' ');
}

function buildRbcLines(entries: Array<[number, number]>): string[] {
  const out: string[] = [];
  for (let i = 0; i < entries.length; i += 8) {
    const chunk = entries.slice(i, i + 8);
    const pairs = chunk.map(([a, v]) => ` ${pad3(a)} ${pad3(v)}`).join('');
    out.push(`M  RBC${pad3(chunk.length)}${pairs}`);
  }
  return out;
}

function buildAlsLine(atom: number, negate: boolean, atoms: string[]): string {
  // Layout: `M  ALS aaa nnn e ssss ssss ...` — note the single space between
  // the "M  ALS" keyword and the atom-index field. The strict BIOVIA spec omits
  // it, but RDKit's V2000 parser reads the atom-index field starting at col 8
  // (matching Marvin/Indigo/OpenBabel output). Without the space, the parser
  // misreads the fields and the resulting qmol throws downstream on aromatic
  // perception.
  const sym = atoms.map((s) => s.padEnd(4, ' ')).join('');
  return `M  ALS ${pad3(atom)}${pad3(atoms.length)} ${negate ? 'T' : 'F'} ${sym}`;
}

/**
 * Rewrites V2000 query-atom aliases that Indigo / ChemAxon / Ketcher emit but
 * RDKit silently downgrades to dummy atoms. Returns the original molblock when
 * input is not a V2000 molblock or no known aliases are present.
 *
 * Recognized aliases (see QUERY_ALIASES): ACY/Acyclic, CYC/Cyclic, Alkyl, Aryl,
 * Heterocyclyl, Heteroaryl, Heteroatom, X/Halogen, M/Metal. Sources checked:
 *   - the atom-symbol field at columns 32-34 of each atom line, and
 *   - the CTfile alias block `A  <atomIdx>\n<aliasText>` after the bond block.
 *
 * Each rewrite replaces the symbol with A/Q/C/* and appends an `M  RBC` and/or
 * `M  ALS` line carrying the actual query constraint.
 */
export function translateQueryAliasesV2000(molblock: string): string {
  if (!molblock || !isMolBlock(molblock) || molblock.includes('V3000'))
    return molblock;

  const lines = molblock.split(/\r?\n/);
  if (lines.length < 5)
    return molblock;
  const counts = lines[3];
  const nAtoms = parseInt(counts.substring(0, 3), 10);
  const nBonds = parseInt(counts.substring(3, 6), 10);
  if (!Number.isFinite(nAtoms) || !Number.isFinite(nBonds))
    return molblock;
  const atomStart = 4;
  const bondStart = atomStart + nAtoms;
  const propsStart = bondStart + nBonds;

  const rewrites = new Map<number, AliasTranslation>();
  const warn = (atomIdx: number, alias: string) => console.warn(
    `translateQueryAliasesV2000: alias "${alias}" on atom ${atomIdx} cannot be ` +
    `represented in V2000 (${UNSUPPORTED_ALIASES[alias]}); leaving untranslated`);

  for (let i = 0; i < nAtoms; i++) {
    const line = lines[atomStart + i];
    if (!line || line.length < 34)
      continue;
    const sym = line.substring(31, 34).trim();
    const t = QUERY_ALIASES[sym];
    if (t)
      rewrites.set(i + 1, t);
    else if (UNSUPPORTED_ALIASES[sym])
      warn(i + 1, sym);
  }

  const keptProps: string[] = [];
  for (let i = propsStart; i < lines.length; i++) {
    const line = lines[i];
    if (line === undefined)
      continue;
    const aliasMatch = /^A\s+(\d+)\s*$/.exec(line);
    if (aliasMatch && i + 1 < lines.length) {
      const idx = parseInt(aliasMatch[1], 10);
      const text = lines[i + 1].trim();
      const t = QUERY_ALIASES[text];
      if (t) {
        rewrites.set(idx, t);
        i++;
        continue;
      }
      if (UNSUPPORTED_ALIASES[text])
        warn(idx, text);
      keptProps.push(line, lines[i + 1]);
      i++;
      continue;
    }
    keptProps.push(line);
  }

  if (rewrites.size === 0)
    return molblock;

  for (const [idx, t] of rewrites.entries()) {
    const li = atomStart + idx - 1;
    const line = lines[li];
    if (!line)
      continue;
    lines[li] = line.substring(0, 31) + t.symbol.padEnd(3, ' ') + line.substring(34);
  }

  const rbcEntries: Array<[number, number]> = [];
  const alsEntries: Array<{atom: number; negate: boolean; atoms: string[]}> = [];
  for (const [idx, t] of rewrites.entries()) {
    if (typeof t.rbc === 'number')
      rbcEntries.push([idx, t.rbc]);
    if (t.als)
      alsEntries.push({atom: idx, ...t.als});
  }

  const newProps = [...buildRbcLines(rbcEntries),
    ...alsEntries.map((e) => buildAlsLine(e.atom, e.negate, e.atoms))];

  // Drop both `M  END` (we re-add it) and any blank lines preserved from the
  // input (e.g. the trailing empty string left by Ketcher's trailing newline).
  // A blank line between the bond block and M-lines confuses the V2000 parser.
  const filteredProps = keptProps.filter((l) => l.trim() !== '' && !/^M\s+END\s*$/.test(l));
  const out = [...lines.slice(0, propsStart), ...filteredProps, ...newProps, 'M  END'];
  return out.join('\n');
}
